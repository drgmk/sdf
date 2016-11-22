import gzip
import re
import copy
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from . import filter
from . import utils
from .utils import SdfError
from . import config as cfg

c_micron = u.micron.to(u.Hz,equivalencies=u.spectral())

class Spectrum(object):
    """Basic spectral class (ha!).

    This class is used as the basis for ObsSpectrum and
    ModelSpectrum
    """
    
    def synthphot(self,filters):
        """Return synthetic photometry of spectrum for given filters.
            
        Here is where disambiguation of colours and indices (i.e.
        Stromgren) is done, so all high-level calls for synthetic
        photometry need to be done through Spectrum classes. A basic
        assumption is that colours/indices only used quoted values.

        First return value is fnu in Jy, second is colour correction.
        """
        
        fnus = np.array([],dtype=float)
        ccs = np.array([])  # dtype will be object unless all ccs defined

        # ensure we have a tuple/list
        if not isinstance(filters,(tuple,list)):
            filt = (filters,)
        else:
            filt = filters
        
        for f in filt:
        
            # a Filter object
            if isinstance(f,filter.Filter):
                fnu,cc = f.synthphot(self)
                fnus = np.append( fnus, fnu )
                ccs = np.append( ccs, cc)

            # a Colour object
            elif isinstance(f,filter.Colour):
                fnus = np.append(fnus, self.synthphot_mag_combo(f.filters,
                                                                f.weights ) )
                ccs = np.append(ccs, None)

            # a string (that might be a colour/index)
            elif isinstance(f,str):
                
                # colour/index
                if filter.iscolour(f):
                    col = filter.Colour.get(f)
                    fnus = np.append(fnus, self.synthphot_mag_combo(col.filters,
                                                                    col.weights ) )
                    ccs = np.append(ccs, None)

                else:
                    ftmp = filter.Filter.get(f)
                    fnu,cc = ftmp.synthphot(self)
                    fnus = np.append( fnus, fnu )
                    ccs = np.append( ccs, cc)
                        
            else:
                raise SdfError("synthphot requires either filter class or\
                               string, not {} or type {}".format(f,type(f)))
                    
        return fnus,ccs


    def synthphot_mag_combo(self,filters,weights):
        """Synthetic photometry for an arbitrary set of filters.
            
        Magnitudes are added using the weights given. E.g. the Stromgren
        M1 index (v-b)-(b-y), which is (v-2b+y), which would be done by
        passing the list of filters ['VS','BS','YS'] and weights 
        [1,-2,1].
        
        If mag is True, the synthetic photometry is converted to mag
        before the weights are applied.
        """
    
        if len(filters) != len(weights):
            raise SdfError("length of filters ({}) needs to be same as weights ({})".format(filters,weights))
        
        value = 0.0
        for i,f in enumerate(filters):
            filt = filter.Filter.get(f)
            fnu,_ = filt.synthphot(self)
            mag = filt.flux2mag(fnu)
            value += mag * weights[i]

        return value


    def sort(self,key):
        """Sort wavelength, frequency, and flux arrays by 'nu' or 'wave'.
            
        Duplicate wavelength values are also removed during the sorting.
        """
        if key == 'nu':
            _,srt = np.unique( self.nu_hz, return_index=True )
        elif key == 'wave':
            _,srt = np.unique( self.wavelength, return_index=True )
        else:
            raise SdfError("pass either 'nu' or 'wave' to sort by")

        # check we need to do something
        if len(self.nu_hz) == len(self.nu_hz[srt]):
            if np.all( np.equal( self.nu_hz, self.nu_hz[srt] ) ):
                return
            
        self.nu_hz = self.nu_hz[srt]
        self.wavelength = self.wavelength[srt]
        if hasattr(self,'fnujy'):
            self.fnujy = self.fnujy[srt]
        if hasattr(self,'e_fnujy'):
            if self.e_fnujy is not None:
                self.e_fnujy = self.e_fnujy[srt]
        if hasattr(self,'fnujy_sr'):
            self.fnujy_sr = self.fnujy_sr[srt]


    def fill_wave2hz(self):
        """Fill nu_hz array from wavelength
        """
        self.nu_hz = c_micron / self.wavelength
    

    def fill_hz2wave(self):
        """Fill wavelength array from nu_hz
        """
        self.wavelength = c_micron / self.nu_hz
    

    def copy(self):
        """Return a copy"""
        
        return copy.deepcopy(self)


    def extend_power_law(self,max_wav_micron=cfg.models['max_wav_micron'],
                         index=-2.):
        """Extrapolate a spectrum to a longer maximum wavelength.
            
        Extrapolate from the last few points using a power law, default
        is Rayleigh-Jeans (-2).
        """
    
        if max_wav_micron <= np.max(self.wavelength):
            raise ValueError("max_wav_micron {} should be longer than current\
                             max of {}".format(max_wav_micron,self.wavelength))
    
        self.sort('wave')
        self.wavelength = np.append(self.wavelength,max_wav_micron)
        self.nu_hz = np.append(self.nu_hz,c_micron/max_wav_micron)
        i = -5
        
        if hasattr(self,'fnujy'):
            meanfnu = np.mean(self.wavelength[i:-1]**(-1.*index) * self.fnujy[i+1:])
            self.fnujy = np.append(self.fnujy,meanfnu * max_wav_micron**index)
        
        if hasattr(self,'fnujy_sr'):
            meanfnu_sr = np.mean(self.wavelength[i:-1]**(-1.*index) * self.fnujy_sr[i+1:])
            self.fnujy_sr = np.append(self.fnujy_sr,meanfnu_sr * max_wav_micron**index)
        
        if hasattr(self,'e_fnujy'):
            raise SdfError("Cannot extrapolate spectrum with uncertainties!")


    def fill_gaps(self,npt=100,log=True):
        """Fill out a spectrum to be less sparse, interpolating fluxes.
            
        Points will be added until there are more than Npt per decade
        (if log=True), or their linear spacing is less than 1/Npt.
        """
    
        self.sort('wave')
        wave = self.wavelength
        spacing = 1./float(npt)
        if log:
            arr = np.log10(self.wavelength)
        else:
            arr = self.wavelength
        newarr = np.array([])
        while True:
            gaps = np.where( arr[1:] - arr[0:-1] > 1.1 * spacing )[0]
            if len(gaps) == 0:
                break
            # fill the first gap
            darr = arr[gaps[0]+1] - arr[gaps[0]]
            #print("filling {} to {} of size {}".format(arr[gaps[0]],arr[gaps[0]+1],darr))
            ninsert = np.ceil( darr/spacing )
            #print("with {} points".format(ninsert))
            insert = arr[gaps[0]] + (np.arange(ninsert-1)+1) * spacing
            #print("that are {}".format(insert))
            arr = np.insert(arr,gaps[0]+1,insert)
            newarr = np.append(newarr,insert)

        # fill these now, fnujy length will be checked by setter
        if log:
            arr = 10**newarr
        self.wavelength = np.append(self.wavelength,arr)
        self.nu_hz = np.append(self.nu_hz,c_micron/arr)

        # now add nu_hz and interpolate fnujy, if interpolation is in log space
        # remember that arr is already there
        if hasattr(self,'fnujy'):
            if log:
                fnujy_add = 10**np.interp(newarr,
                                          np.log10(wave),
                                          np.log10(self.fnujy))
            else:
                fnujy_add = np.interp(newarr,wave,self.fnujy)
            self.fnujy = np.append(self.fnujy,fnujy_add)

        if hasattr(self,'fnujy_sr'):
            if log:
                fnujy_sr_add = 10**np.interp(newarr,
                                             np.log10(wave),
                                             np.log10(self.fnujy_sr))
            else:
                fnujy_sr_add = np.interp(arr,wave,self.fnujy_sr)
            self.fnujy_sr = np.append(self.fnujy_sr,fnujy_sr_add)

        if hasattr(self,'e_fnujy'):
            raise SdfError("Cannot fill gaps in spectrum with uncertainties!")


    @property
    def wavelength(self):
        return self._wavelength
    @wavelength.setter
    def wavelength(self,value):
        self._wavelength = utils.validate_1d(value,None)

    @property
    def nu_hz(self):
        return self._nu_hz
    @nu_hz.setter
    def nu_hz(self,value):
        self._nu_hz = utils.validate_1d(value,None)


class ObsSpectrum(Spectrum):
    """Class used for observed spectra.

    Includes uncertainty and bibcode attributes.
    """
    def __init__(self,wavelength=None,nu_hz=None,fnujy=None,
                 e_fnujy=None,e_absolute=None,bibcode=None,
                 irradiance=None):
        self.wavelength = wavelength
        self.nu_hz = nu_hz
        self.fnujy = fnujy
        self.e_fnujy = e_fnujy
        self.e_absolute = e_absolute
        self.bibcode = bibcode
        self.irradiance = irradiance
    
        # fill wavelength/frequency if possible
        if self.wavelength is None and self.nu_hz is not None:
            self.fill_hz2wave()
        if self.nu_hz is None and self.wavelength is not None:
            self.fill_wave2hz()


    @classmethod
    def read_cassis(cls,file,module_split=False):
        """Read a file from the CASSIS database of IRS spectra.
        
        Check for nan values in uncertainty arrays.
        
        Optionally return a list of ObsSpectrum objects, split by
        IRS module (SL1, SL2, LL1, LL2), since these might need to
        be normalised individually.
            
        As the spectrum will be normalised as part of the fitting
        precedure, it's not clear which ucertainties to include.
        None are wavelength independent (and thus ignorable), so
        use the 'total' (rms + systematic).
        """
    
        # open the file and put it in a Table
        fh = fits.open(file)
        t = Table(fh[0].data)
    
        # get the column names
        coldefs = fh[0].header['COL*DEF']
        for i,name in enumerate(t.colnames):
            try:
                t[name].name = coldefs[i]
            except IndexError:
                pass
        # fix missing 16th column definition
        if t.colnames[-1] == 'col15' and fh[0].header['CASVE'] == '7':
            t[t.colnames[-1]].name = 'flag4'
        fh.close()

        # optionally split into modules
        modules = np.unique(t['module']).data
        if module_split and len(modules) > 1:
            s = [ObsSpectrum() for i in range(len(modules))]
            for i,mod in enumerate(modules):
                ok = (t['module'] == mod) & (np.isfinite(t['error (RMS+SYS)']))
                s[i].wavelength = t['wavelength'].data[ok]
                s[i].nu_hz = c_micron / t['wavelength'].data[ok]
                s[i].fnujy = t['flux'].data[ok] \
                             * u.Unit(fh[0].header['BUNIT']).to('Jy')
                s[i].e_fnujy = t['error (RMS+SYS)'].data[ok]
                s[i].bibcode = '2011ApJS..196....8L'
                s[i].sort('wave')
            return s

        else:
            self = cls()
            ok = np.isfinite(t['error (RMS+SYS)'])
            self.wavelength = t['wavelength'].data[ok]
            self.nu_hz = c_micron / t['wavelength'].data[ok]
            self.fnujy = t['flux'].data[ok] * u.Unit(fh[0].header['BUNIT']).to('Jy')
            self.e_fnujy = t['error (RMS+SYS)'].data[ok]
            self.bibcode = '2011ApJS..196....8L'
            self.sort('wave')
            return self


    @classmethod
    def read_file_of_type(cls,file,type,module_split=False):
        """Read a spectrum file, sorting out how to do it."""
        
        if type == 'irsstare':
            return ObsSpectrum.read_cassis(file,module_split=module_split)


    @classmethod
    def read_sdb_file(cls,file,module_split=False):
        """Read a sdb rawphot file and return a tuple of spectra.
        
        The file format is set by sdb_getphot.py, and is ascii.ipac.
        To get any photometry use photometry.get_sdb_file.
        
        There could be multiple spectra of a given type, so the
        keywords can have extra text in them to allow the keywords
        to be unique (e.g. 'irsstare0').
        """

        # get the data and keywords
        p = Table.read(file,format='ascii.ipac')
        kw = p.meta['keywords']
        
        # get spectra that match these names
        spnames = ['irsstare']
        s = []
        for sp in spnames:
            for key in kw.keys():
                if sp in key:
                    s.append( ObsSpectrum.read_file_of_type(kw[key]['value'],type=sp,
                                                            module_split=module_split) )
    
        # flatten the list in case we were appending lists to s
        if isinstance(s,ObsSpectrum):
            return (s,)
        elif len(s) > 0:
            sp = [item for sublist in s for item in sublist]
            return tuple(sp)
        else:
            return None
                

    @classmethod
    def vega_stis(cls,file='alpha_lyr_stis_008.fits'):
        """Return the CALSPEC Vega spectrum.
        
        The FITS file contains improper units, it was easier to just
        set these blank and fill them in below.
        """
    
        self = cls()
        t=Table.read(cfg.file['model_root']+file)
        t['WAVELENGTH'].unit = u.AA
        t['FLUX'].unit = u.Unit('erg/(s*cm*cm*AA)')

        self.wavelength = t['WAVELENGTH'].to(u.micron).value
        self.nu_hz = c_micron / self.wavelength
        self.fnujy = t['FLUX'].to(u.jansky,equivalencies=           \
                                  u.spectral_density(t['WAVELENGTH'])).value

        neg = self.fnujy <= 0.0
        self.fnujy[neg] = cfg.tiny

        self.fill_irradiance()
        return self


    @classmethod
    def vega_rieke(cls):
        """Return the "Vega" spectrum from Rieke et al. (2008)
            
        This spectrum has slightly lower fluxes than stated in the
        paper (7.14 vs. 7.15Jy @ 23.675um), and the adopted calibration
        has 7.17Jy, so to ensure that the calibration is as precise as
        possible multiply the spectrum up by 7.17/7.14. A possible
        reason for the discrepancy is the F_lambda to F_nu conversion,
        as there is also a small difference in the 10.6um values in
        Table 1 of that paper.
        """
        
        self = cls()
        t = Table.read(cfg.file['model_root']+'rieke08-vega-sun.txt',
                       format='ascii.cds')
        wave_um = t['Wave'].to(u.micron).value
        _,srt = np.unique(wave_um,return_index=True)
        srt = np.flipud(srt)
        self.wavelength = wave_um[srt]
        self.nu_hz = c_micron / wave_um[srt]
        fnu = t['Vega'].to(u.jansky,
                           equivalencies=u.spectral_density(t['Wave'])).value
        self.fnujy = fnu[srt]
        self.fnujy *= 7.17 / 7.1375583083242953
        
        neg = self.fnujy <= 0.0
        self.fnujy[neg] = cfg.tiny

        self.fill_irradiance()
        return self

    
    def fill_irradiance(self):
        """Compute the irradiance in W/m2."""
        
        self.irradiance = utils.sdf_int(self.fnujy,self.nu_hz) * 1e26 * u.W/u.m**2

        
    @property
    def e_fnujy(self):
        return self._e_fnujy
    @e_fnujy.setter
    def e_fnujy(self,value):
        if self.fnujy is None:
            expected_len = None
        else:
            expected_len = len(self.fnujy)
        self._e_fnujy = utils.validate_1d(value,expected_len)

        
    @property
    def bibcode(self):
        return self._bibcode
    @bibcode.setter
    def bibcode(self, value):
        self._bibcode = utils.validate_string(value)

    @property
    def e_absolute(self):
        return self._e_absolute
    @e_absolute.setter
    def e_absolute(self, value):
        self._e_absolute = utils.validate_float(value)

    @property
    def fnujy(self):
        return self._fnujy
    @fnujy.setter
    def fnujy(self,value):
        if self.nu_hz is None:
            expected_len = None
        else:
            expected_len = len(self.nu_hz)
        self._fnujy = utils.validate_1d(value,expected_len)

                
class ModelSpectrum(Spectrum):
    """Class used for model spectra

    Flux density is per steradian
    """
    
    def __init__(self,name=None,parameters=None,param_values=None,
                 wavelength=None,nu_hz=None,fnujy_sr=None):
        self.name = name
        self.parameters = parameters
        self.param_values = param_values
        self.wavelength = wavelength
        self.nu_hz = nu_hz
        self.fnujy_sr = fnujy_sr

        # fill wavelength/frequency if possible
        if self.wavelength is None and self.nu_hz is not None:
            self.fill_hz2wave()
        if self.nu_hz is None and self.wavelength is not None:
            self.fill_wave2hz()


    def resample(self,resolution=100,kernel=None):
        """Resample a model spectrum to a different resolution.
            
        Return the kernel that was used in case we want to use it again.
        """
    
        # get the new wavelength grid
        wav = np.power(10,np.arange(np.log10(np.min(self.wavelength)),
                                    np.log10(np.max(self.wavelength)),
                                    1.0/float(resolution)) )
    
        # get the convolution/resampling matrix (unless we got it)
        self.sort('wave')
        if kernel is None:
            kernel = utils.resample_matrix(self.wavelength,wav,
                                           new_R=resolution)

        self.wavelength = wav
        self.fill_wave2hz()
        self.fnujy_sr = self.fnujy_sr * kernel

        return kernel
    

    @classmethod
    def bnu_wave_micron(cls,wave_um,temp,lam0=None,beta=None):
        """Return a (modified) black body spectrum"""

        self = cls()
        
        self.name = 'bb'
        self.parameters = ['Temp']
        self.param_values = [temp]
        fnujysr = utils.bnu_wav_micron(wave_um,temp)
        srt = np.flipud(np.argsort(wave_um))
        self.wavelength = wave_um[srt]
        self.nu_hz = c_micron / wave_um[srt]
        self.fnujy_sr = fnujysr[srt]
        
        # modified if desired
        if lam0 is not None and beta is not None:
            lam0 = float(lam0)
            self.name = 'modbb'
            self.parameters = ['Temp','lam0','beta']
            self.param_values = [temp,lam0,beta]
            ok = self.wavelength > lam0
            self.fnujy_sr[ok] /= np.power(self.wavelength[ok]/lam0,beta)
        
        return self

    @classmethod
    def bnu_nu_hz(cls,nu_hz,temp,lam0=None,beta=None):
        """Return a black body spectrum
        """

        self = cls()
        wave_um = c_micron / nu_hz
        return self.bnu_wave_micron(wave_um,temp,lam0=None,beta=None)


    @classmethod
    def read_phoenix(cls,file):
        """Return the spectrum from a PHOENIX model
        
        Format is wavelength (A), flam (erg/s/cm^2/cm), bnu 
        (erg/s/cm^2/cm), and a bunch of other stuff we don't care 
        about. The flam column is log10(flam).

        A factor of pi is needed to get the PHOENIX models into Jy/sr
        (and get stellar radii correct). The docs
        https://phoenix.ens-lyon.fr/Grids/FORMAT suggest that
        (Rstar/d)^2 is the normalisation, so divide spectra by pi here.
        Comparing the results with Kurucz models shows this is correct.
        
        TODO: Can this go faster? BT-Settl files are epic

        """
        self = cls()

        self.name = 'phoenix'
        self.parameters = ['Teff','logg','MH']
        if 'BT-Settl' in file:
            par = re.search('lte([0-9]+)-([0-9.]+)([+-][0-9.]+)[a]([+-][0-9.]+).BT-Settl.7.bz2',file)
            teff = float(par.groups()[0])*100.0
            logg = float(par.groups()[1])
            mh = float(par.groups()[2])
        self.param_values = {'Teff':teff,'logg':logg,'MH':mh}
        
        c = lambda s: s.decode().replace("D", "E") # deal with fortran 1.0D+01
        w,f = np.loadtxt(file,dtype=float,usecols=(0,1),
                         converters={0:c,1:c},unpack=True)
        f = np.power(10,f)
        _,srt = np.unique( w, return_index=True )
        srt = np.flipud(srt)
        wav = u.Quantity( w[srt], u.Angstrom )
        flam = u.Quantity( f[srt]/np.pi, u.erg/u.s/u.cm**2/u.cm) # divide by pi
        self.wavelength = wav.to(u.micron).value
        self.nu_hz = wav.to( 'Hz', equivalencies=u.spectral() ).value
        self.fnujy_sr = flam.to(u.jansky,
                                equivalencies=u.spectral_density(wav) ).value
        self.extend_power_law()
        self.fill_gaps()
        self.sort('wave')
        return self


    def read_kurucz(file):
        """Read a grid of Castelli & Kurucz models from a specific file.
        Returns a tuple with (teff,logg,[M/H],models), where the first
        three are arrays with the parameters of each model, and models
        is a list of ModelSpectrum objects holding the spectra.
        
        TODO: not really clear whether this belongs here or in 
        SpecModel (or elsewhere). It's only here because it's similar
        to read_phoenix, which deals with single files.
        """

        fh = open(file,'r')
        lines = fh.readlines()

        cols = ((0,10),(10,20),(20,30),(30,40),(40,50),(50,60),(60,70),(70,80))

        # create list of models, and arrays of their parameters
        teff = np.array([],dtype=float)
        logg = np.array([],dtype=float)
        mh = np.array([],dtype=float)
        models = np.array([],dtype=object)
        wave = np.array([],dtype=float)
        model = []
        for i,l in enumerate(lines):
            l = l.rstrip()

            # wavelength grid comes first, starting on a line with only numbers
            if re.search('[a-zA-Z]',l) is None and len(teff) == 0 and 'TEFF' not in l:
                for c in cols:
                    if c[0] < len(l):
                        wave = np.append( wave, float( l[c[0]:c[1]] ) )
                continue

            # skip to next line if wave hasn't been started yet
            if len(wave) == 0:
                continue
            
            # new models have TEFF on first line
            if 'TEFF' in l:
                teffLogg = re.search('TEFF\s+(\d+\.\d*)\s+GRAVITY\s+(\d+\.\d*)',l)
                teff = np.append(teff,float(teffLogg.groups()[0]))
                logg = np.append(logg,float(teffLogg.groups()[1]))
                metalAlpha = re.search('\[([\s+-]?\d+\.\d+)([aAbB]?)\]?',l)
                mh = np.append(mh,float(metalAlpha.groups()[0]))
                if metalAlpha.groups()[1] != '':
                    raise SdfError("Alpha non-nothing ({}), wtf?".\
                                   format(metalAlpha.groups()))
            else:
                for c in cols:
                    if c[0] < len(l):
                        model.append( float( l[c[0]:c[1]] ) )
                
            # add complete model to list
            if len(model) > 0 and ('TEFF' in l or i == len(lines)-1):
                self = ModelSpectrum()
                self.name = 'kurucz'
                self.parameters = ['Teff','logg','MH']
                # if next model has been started we want
                # params from last model
                if 'TEFF' in l:
                    self.param_values = {'Teff':teff[-2],
                                         'logg':logg[-2],
                                         'MH':mh[-2]}
                else:
                    self.param_values = {'Teff':teff[-1],
                                         'logg':logg[-1],
                                         'MH':mh[-1]}
                
                self.wavelength = (wave * u.nm).to('micron').value
                self.nu_hz = c_micron / self.wavelength
                model = model[0:len(wave)]
                self.fnujy_sr = np.array(model,dtype=float) * 4e23
                # turn zeros into small numbers to avoid divide errors
                small = (self.fnujy_sr <= 0.0)
                self.fnujy_sr[small] = cfg.tiny
                # extend to long wavelengths, fill out, and sort again
                self.extend_power_law()
                self.fill_gaps()
                self.sort('wave')
                models = np.append(models,self)
                model = []

        return models,teff,logg,mh
    
    
    @property
    def fnujy_sr(self):
        return self._fnujy_sr
    @fnujy_sr.setter
    def fnujy_sr(self,value):
        if self.nu_hz is None:
            expected_len = None
        else:
            expected_len = len(self.nu_hz)
        self._fnujy_sr = utils.validate_1d(value,expected_len)
