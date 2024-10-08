"""Classes for spectra.
   
There are two classes for spectra, `ObsSpectrum` for real spectra, and
`ModelSpectrum` for models. Both are baed on the basic `Spectrum` class
so have a number of attributes and methods in common.

One of the most important methods is `synthphot`, which is how all
synthetic photometry calls are made.
"""

import os
import re
import copy
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from . import filter
from . import utils
from . import config as cfg

c_micron = u.micron.to(u.Hz, equivalencies=u.spectral())


class Spectrum(object):
    """Basic spectral class.

    Used as the basis for `ObsSpectrum` and `ModelSpectrum`.
    """
    
    def synthphot(self, filters):
        """Return synthetic photometry of spectrum for given filters.
            
        Here is where disambiguation of colours and indices (i.e.
        Stromgren) is done, so all high-level calls for synthetic
        photometry need to be done through Spectrum classes. A basic
        assumption is that colours/indices only used quoted values.

        First return value is fnu in Jy, second is colour correction.
        """
        
        fnus = np.array([], dtype=float)
        ccs = np.array([])  # dtype will be object unless all ccs defined

        # ensure we have a tuple/list
        if not isinstance(filters, (tuple, list)):
            filt = (filters, )
        else:
            filt = filters
        
        for f in filt:
        
            # a Filter object
            if isinstance(f, filter.Filter):
                fnu, cc = f.synthphot(self)
                fnus = np.append(fnus, fnu)
                ccs = np.append(ccs, cc)

            # a Colour object
            elif isinstance(f, filter.Colour):
                fnus = np.append(fnus, self.synthphot_mag_combo(f.filters,
                                                                f.weights))
                ccs = np.append(ccs, None)

            # a string (that might be a colour/index)
            elif isinstance(f, str):
                
                # colour/index
                if filter.iscolour(f):
                    col = filter.Colour.get(f)
                    fnus = np.append(fnus, self.synthphot_mag_combo(col.filters,
                                                                    col.weights))
                    ccs = np.append(ccs, None)

                else:
                    ftmp = filter.Filter.get(f)
                    fnu, cc = ftmp.synthphot(self)
                    fnus = np.append(fnus, fnu)
                    ccs = np.append(ccs, cc)
                        
            else:
                raise utils.SdfError("synthphot requires either filter class or\
                               string, not {} or type {}".format(f, type(f)))
                    
        return fnus, ccs

    def synthphot_mag_combo(self, filters, weights):
        """Synthetic photometry for an arbitrary set of filters.
            
        Magnitudes are added using the weights given. E.g. the Stromgren
        M1 index (v-b)-(b-y), which is (v-2b+y), which would be done by
        passing the list of filters ['VS', 'BS', 'YS'] and weights 
        [1, -2, 1].
        
        If mag is True, the synthetic photometry is converted to mag
        before the weights are applied.
        """
    
        if len(filters) != len(weights):
            raise utils.SdfError("length of filters ({}) needs to be same as weights ({})".format(filters, weights))
        
        value = 0.0
        for i, f in enumerate(filters):
            filt = filter.Filter.get(f)
            fnu, _ = filt.synthphot(self)
            mag = filt.flux2mag(fnu)
            value += mag * weights[i]

        return value

    def sort(self, key):
        """Sort wavelength, frequency, and flux arrays by 'nu' or 'wave'.
            
        Duplicate wavelength values are also removed during the sorting.
        """
        if key == 'nu':
            _, srt = np.unique(self.nu_hz, return_index=True)
        elif key == 'wave':
            _, srt = np.unique(self.wavelength, return_index=True)
        else:
            raise utils.SdfError("pass either 'nu' or 'wave' to sort by")

        # check we need to do something
        if len(self.nu_hz) == len(self.nu_hz[srt]):
            if np.all(np.equal(self.nu_hz, self.nu_hz[srt])):
                return
            
        self.nu_hz = self.nu_hz[srt]
        self.wavelength = self.wavelength[srt]
        if hasattr(self, 'fnujy'):
            self.fnujy = self.fnujy[srt]
        if hasattr(self, 'e_fnujy'):
            if self.e_fnujy is not None:
                self.e_fnujy = self.e_fnujy[srt]
        if hasattr(self, 'fnujy_sr'):
            self.fnujy_sr = self.fnujy_sr[srt]

    def fill_irradiance(self):
        """Compute the irradiance.
            
        Units are W/m2 for an ObsSpectrum, W/m2/sr for a ModelSpectrum.
        """
        
        self.sort('nu')
        if isinstance(self, ObsSpectrum):
            self.irradiance = 1e-26 * utils.sdf_int(self.fnujy, self.nu_hz)
        elif isinstance(self, ModelSpectrum):
            self.irradiance_sr = 1e-26 * utils.sdf_int(self.fnujy_sr, self.nu_hz)

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

    def extend_power_law(self, max_wav_micron=cfg.models['max_wav_micron'],
                         index=-2.):
        """Extrapolate a spectrum to a longer maximum wavelength.
            
        Extrapolate from the last few points using a power law, default
        is Rayleigh-Jeans (-2).
        """
    
        if max_wav_micron <= np.max(self.wavelength):
            raise ValueError("max_wav_micron {} should be longer than current\
                             max of {}".format(max_wav_micron, self.wavelength))
    
        self.sort('wave')
        self.wavelength = np.append(self.wavelength, max_wav_micron)
        self.nu_hz = np.append(self.nu_hz, c_micron/max_wav_micron)
        i = -5  # use average of -this many points at end of spectrum
        
        if hasattr(self, 'fnujy'):
            meanfnu = np.mean(self.wavelength[i:-1]**(-1.*index) * self.fnujy[i+1:])
            self.fnujy = np.append(self.fnujy, meanfnu * max_wav_micron**index)
        
        if hasattr(self, 'fnujy_sr'):
            meanfnu_sr = np.mean(self.wavelength[i:-1]**(-1.*index) * self.fnujy_sr[i+1:])
            self.fnujy_sr = np.append(self.fnujy_sr, meanfnu_sr * max_wav_micron**index)
        
        if hasattr(self, 'e_fnujy'):
            raise utils.SdfError("Cannot extrapolate spectrum with uncertainties!")

    def fill_gaps(self, npt=100, log=True):
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
            gaps = np.where(arr[1:] - arr[0:-1] > 1.1 * spacing)[0]
            if len(gaps) == 0:
                break
            # fill the first gap
            darr = arr[gaps[0]+1] - arr[gaps[0]]
            # print("filling {} to {} of size {}".format(arr[gaps[0]], arr[gaps[0]+1], darr))
            ninsert = np.ceil(darr/spacing)
            # print("with {} points".format(ninsert))
            insert = arr[gaps[0]] + (np.arange(ninsert-1)+1) * spacing
            # print("that are {}".format(insert))
            arr = np.insert(arr, gaps[0]+1, insert)
            newarr = np.append(newarr, insert)

        # fill these now, fnujy length will be checked by setter
        if log:
            arr = 10**newarr
        self.wavelength = np.append(self.wavelength, arr)
        self.nu_hz = np.append(self.nu_hz, c_micron/arr)

        # now add nu_hz and interpolate fnujy, if interpolation is in log space
        # remember that arr is already there
        if hasattr(self, 'fnujy'):
            if log:
                fnujy_add = 10**np.interp(newarr,
                                          np.log10(wave),
                                          np.log10(self.fnujy))
            else:
                fnujy_add = np.interp(newarr, wave, self.fnujy)
            self.fnujy = np.append(self.fnujy, fnujy_add)

        if hasattr(self, 'fnujy_sr'):
            if log:
                fnujy_sr_add = 10**np.interp(newarr,
                                             np.log10(wave),
                                             np.log10(self.fnujy_sr))
            else:
                fnujy_sr_add = np.interp(arr, wave, self.fnujy_sr)
            self.fnujy_sr = np.append(self.fnujy_sr, fnujy_sr_add)

        if hasattr(self, 'e_fnujy'):
            raise utils.SdfError("Cannot fill gaps in spectrum with uncertainties!")

    def resample(self, wav, kernel=None):
        """Resample a model spectrum to a different resolution.
            
        Parameters
        ----------
        wav : numpy.ndarray
            Wavelength grid to resample to.
        kernel : 2d numpy.ndarray, optional
            Pre-computed kernel for convolution.
            
        Returns
        -------
        kernel : 2d numpy ndarray
            The kernel that was generated for this convolution.

        See Also
        --------
        sdf.utils.resample_matrix : Where the kernel is computed
        """
    
        # get the convolution/resampling matrix (unless we got it)
        self.sort('wave')
        if kernel is None:
            kernel = utils.resample_matrix(self.wavelength, wav)

        self.wavelength = wav
        self.fill_wave2hz()
        if hasattr(self, 'fnujy'):
            self.fnujy = self.fnujy * kernel
        if hasattr(self, 'e_fnujy'):
            if self.e_fnujy is not None:
                self.e_fnujy = self.e_fnujy * kernel
        if hasattr(self, 'fnujy_sr'):
            self.fnujy_sr = self.fnujy_sr * kernel

        return kernel

    @property
    def wavelength(self):
        return self._wavelength

    @wavelength.setter
    def wavelength(self, value):
        self._wavelength = utils.validate_1d(value, None)

    @property
    def nu_hz(self):
        return self._nu_hz

    @nu_hz.setter
    def nu_hz(self, value):
        self._nu_hz = utils.validate_1d(value, None)


class ObsSpectrum(Spectrum):
    """Class used for observed spectra.

    Includes uncertainty and bibcode attributes.
    """
    def __init__(self, wavelength=None, nu_hz=None, fnujy=None,
                 e_fnujy=None, e_absolute=None, bibcode=None,
                 instrument=None, irradiance=None):
        self.wavelength = wavelength
        self.nu_hz = nu_hz
        self.fnujy = fnujy
        self.e_fnujy = e_fnujy
        self.e_absolute = e_absolute
        self.bibcode = bibcode
        self.instrument = instrument
        self.irradiance = irradiance
    
        # fill wavelength/frequency if possible
        if self.wavelength is None and self.nu_hz is not None:
            self.fill_hz2wave()
        if self.nu_hz is None and self.wavelength is not None:
            self.fill_wave2hz()
            
        if self.wavelength is not None and self.fnujy is not None:
            if len(self.wavelength) != len(fnujy):
                raise utils.SdfError('ObsSpectrum: wave/freq and flux arrays different length.')

    @classmethod
    def read_cassis_lores(cls, file, module_split=False):
        """Read a low res file from the CASSIS database of IRS spectra.
        
        Returns a tuple of ObsSpectrum objects, optionally split by
        IRS module (SL1, SL2, LL1, LL2), since these might need to
        be normalised individually.

        Check for nan values in uncertainty arrays.

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
        for i, name in enumerate(t.colnames):
            try:
                t[name].name = coldefs[i]
            except IndexError:
                pass
        # fix missing 16th column definition
        if t.colnames[-1] == 'col15' and fh[0].header['CASVE'] == '7':
            t[t.colnames[-1]].name = 'flag4'
        fh.close()

        # optionally split into modules and sort by wavelength
        modules = np.unique(t['module']).data
        names = np.array(['SL1', 'SL2', 'LL1', 'LL2'])
        module_names = names[modules.astype(int)]
        if module_split and len(modules) > 1:
            mod_sort = []
            s_list = []
            for i, mod in enumerate(modules):
                ok = ((t['module'] == mod) &
                      (np.isfinite(t['flux'])) &
                      (np.isfinite(t['error (RMS+SYS)'])))

                # there might be all nan errors for a module
                if np.any(ok) > 0:
                    s = ObsSpectrum()
                    s.instrument = 'Spitzer IRS '+module_names[i]
                    s.wavelength = t['wavelength'].data[ok]
                    s.nu_hz = c_micron / t['wavelength'].data[ok]
                    s.fnujy = t['flux'].data[ok] \
                        * u.Unit(fh[0].header['BUNIT']).to('Jy')
                    s.e_fnujy = t['error (RMS+SYS)'].data[ok]
                    s.bibcode = '2011ApJS..196....8L'
                    s.sort('wave')
                    s_list.append(s)
                    mod_sort.append(np.min(s.wavelength))

            if len(s_list) == 0:
                raise utils.SdfError('No usable data in {}'.format(file))

            s_list = [s_list[i] for i in np.argsort(mod_sort)]
            return tuple(s_list)

        else:
            self = cls()
            ok = np.isfinite(t['error (RMS+SYS)'])
            if np.any(ok) == 0:
                raise utils.SdfError('No usable data in {}'.format(file))
            self.wavelength = t['wavelength'].data[ok]
            self.nu_hz = c_micron / t['wavelength'].data[ok]
            self.fnujy = t['flux'].data[ok] * u.Unit(fh[0].header['BUNIT']).to('Jy')
            self.e_fnujy = t['error (RMS+SYS)'].data[ok]
            self.bibcode = '2011ApJS..196....8L'
            self.instrument = 'Spitzer IRS'
            self.sort('wave')
            return self,

    @classmethod
    def read_cassis_hires(cls, file, module_split=False):
        """Read a high res file from the CASSIS database of IRS spectra.
            
        Returns a tuple of ObsSpectrum objects, optionally split by IRS
        module (SH1, SH2), since these might need to be normalised
        individually.
        
        Check for nan values in uncertainty arrays.
        """
        
        # open the file and put it in a Table
        fh = fits.open(file)
        t = Table(fh[0].data)
        
        # get the column names
        coldefs = fh[0].header['COL*DEF']
        for i, name in enumerate(t.colnames):
            try:
                t[name].name = coldefs[i]
            except IndexError:
                pass
    
        fh.close()
        
        # optionally split into modules and sort by wavelength
        # module no. not given, so using wavelength of 19.475um
        module = ['SH', 'LH']
        wave_split = t['wavelength'] - 19.475
        split_sign = [-1, 1]
        if module_split:
            s_list = []
            for mod, sign in zip(module, split_sign):

                ok = ((wave_split * sign > 0) &
                      (np.isfinite(t['flux'])) &
                      (np.isfinite(t['flux_error'])))
                no = np.logical_or(t['flux_error'][ok] <= 0,
                                   t['flux'][ok] <= 0)
                ok[ok] = np.invert(no)
            
                # there might be all nan errors or no data
                if np.any(ok) > 0:
                    s = ObsSpectrum()
                    s.instrument = 'Spitzer IRS '+mod
                    s.wavelength = t['wavelength'].data[ok]
                    s.nu_hz = c_micron / t['wavelength'].data[ok]
                    s.fnujy = t['flux'].data[ok]
                    s.e_fnujy = t['flux_error'].data[ok]
                    s.bibcode = '2015ApJS..218...21L'
                    s.sort('wave')
                    s_list.append(s)

            if len(s_list) == 0:
                raise utils.SdfError('No usable data in {}'.format(file))

            return tuple(s_list)
        
        else:
            self = cls()
            ok = np.isfinite(t['flux_error'])
            no = np.logical_or(t['flux_error'][ok] <= 0,
                               t['flux'][ok] <= 0)
            ok[ok] = np.invert(no)
            if np.any(ok) == 0:
                raise utils.SdfError('No usable data in {}'.format(file))
            self.wavelength = t['wavelength'].data[ok]
            self.nu_hz = c_micron / t['wavelength'].data[ok]
            self.fnujy = t['flux'].data[ok]
            self.e_fnujy = t['flux_error'].data[ok]
            self.bibcode = '2015ApJS..218...21L'
            self.instrument = 'Spitzer IRS'
            self.sort('wave')
            return self,

    @classmethod
    def read_csv(cls, file):
        """Read a spectrum in a generic csv file.
        
        Assumes columns are wavelength (um), flux density (Jy), and
        uncertainty (Jy).
        """
        self = cls()
        t = np.genfromtxt(file, delimiter=',', comments='#')
        if t.shape[1] == 2:
            w, f = t.T
            e = 0.1 * np.abs(f)
            print('sdf.spectrum.ObsSpectrum.read_csv assuming 10% flux '
                  'error on spectrum {}'.format(file))
        elif t.shape[1] == 3:
            w, f, e = t.T
        elif t.shape[1] > 3:
            print(f'sdf.spectrum.ObsSpectrum.read_csv using first 3 cols as w, f, e in {file}')
            w, f, e = t.T[0], t.T[1], t.T[2]
        else:
            raise utils.SdfError('file {} has wrong number of columns'.format(file))

        # check for nans that may arise from comment columns
        ok = np.isfinite(w) & np.isfinite(f) & np.isfinite(e)

        self.wavelength = w[ok]
        self.nu_hz = c_micron / w[ok]
        self.fnujy = f[ok]
        self.e_fnujy = e[ok]
        self.bibcode = file
        self.instrument = '?'
        self.sort('wave')
        return self,

    @classmethod
    def read_spex(cls, file):
        """Read a spectrum in a fits file from SpeXTool.
        
        Assumes columns are wavelength, flux density, and uncertainty.
        """
        self = cls()
        
        fh = fits.open(file)
        t = Table(fh[0].data.T)
        t = t[np.isfinite(t['col1'])]
        
        w = t['col0'] * u.Unit(fh[0].header['XUNITS'])
        
        funit = fh[0].header['YUNITS'].replace('Wm', 'W m').replace('um', ' um')
        flam = t['col1'] * u.Unit(funit)
        fnu = flam.to('Jy', equivalencies=u.spectral_density(w))
        elam = t['col2'] * u.Unit(funit)
        enu = elam.to('Jy', equivalencies=u.spectral_density(w))

        self.wavelength = w.value
        self.nu_hz = c_micron / w.value
        self.fnujy = fnu.value
        self.e_fnujy = enu.value
        self.bibcode = file
        self.instrument = 'IRTF/SpeX'
        self.sort('wave')
        return self,

    @classmethod
    def read_file_of_type(cls, file, type_, module_split=False):
        """Read a spectrum file, sorting out how to do it.
            
        The calls should all return tuples of ObsSpectrum objects.
        
        Since we might have low or high resolution spectra, try low res
        first (which will fail).
        """
        
        if type_ == 'irsstare':
            try:
                return ObsSpectrum.read_cassis_lores(file, module_split=module_split)
            except:
                pass

            return ObsSpectrum.read_cassis_hires(file, module_split=module_split)
                
        elif type_ == 'spex':
            return ObsSpectrum.read_spex(file)
                
        elif type_ == 'csv':
            return ObsSpectrum.read_csv(file)

    def read_sdb_file(file, module_split=False, nspec=None):
        """Read a sdb rawphot file and return a tuple of spectra.
        
        The file format is set by sdb_getphot.py, and is ascii.ipac.
        To get any photometry use photometry.get_sdb_file, and for
        keywords use utils.get_sdb_keywords.
        
        The nspec keyword sets now many individual spectra of a given
        type are returned. These can still be split by module. If more
        than one spectrum of a given type exists then just the first is
        returned.
        
        There could be multiple spectra of a given type, so the
        keywords can have extra text in them to allow the keywords
        to be unique (e.g. 'irsstare0').
        """

        # get the data and keywords
        p = Table.read(file, format='ascii.ipac')
        kw = p.meta['keywords']

        if nspec is None:
            nspec = {'irsstare': 1,
                     'spex': 1,
                     'visir': 1,
                     'csv': 10}

        n = {}
        for k in nspec.keys():
            n[k] = 0

        # get spectra that match these names
        s = ()
        for key in kw.keys():
            # shortcut
            if np.sum([i in key for i in ['irsstare', 'spex', 'visir', 'csv']]) == 0:
                continue

            path = ''
            if kw[key]['value'][0] != '/':
                path = cfg.file['spectra']
            
            if 'irsstare' in key:
                if n['irsstare'] + 1 > nspec['irsstare']:
                    continue
                s += ObsSpectrum.read_file_of_type(
                                path+kw[key]['value'], type_='irsstare',
                                module_split=module_split
                                )
                n['irsstare'] += 1
            elif 'spex' in key:
                if n['spex'] + 1 > nspec['spex']:
                    continue
                s += ObsSpectrum.read_file_of_type(
                                path+kw[key]['value'], type_='spex')
                n['spex'] += 1
            elif 'visir' in key:
                if n['visir'] + 1 > nspec['visir']:
                    continue
                s += ObsSpectrum.read_file_of_type(
                                path+kw[key]['value'], type_='csv')
                n['visir'] += 1
            elif 'csv' in key:
                if n['csv'] + 1 > nspec['csv']:
                    continue
                s += ObsSpectrum.read_file_of_type(
                                path+kw[key]['value'], type_='csv')
                n['csv'] += 1

        if len(s) > 0:
            return s
        else:
            return None

    @classmethod
    def vega_stis(cls, file='alpha_lyr_stis_008.fits', file_loc=None):
        """Return the CALSPEC Vega spectrum.
        
        The FITS file contains improper units, it was easier to just
        set these blank and fill them in below.
        
        Parameters
        ----------
        file : str, optional
            Name of the Vega spectrum file.
        file_loc : str, optional
            Location of the file.
        """
    
        self = cls()
        if file_loc is None:
            file_loc = os.path.dirname(os.path.abspath(__file__)) + '/data/calibration/'
    
        t = Table.read(file_loc+file)
        t['WAVELENGTH'].unit = u.AA
        t['FLUX'].unit = u.Unit('erg/(s*cm*cm*AA)')

        self.wavelength = t['WAVELENGTH'].to(u.micron).value
        self.nu_hz = c_micron / self.wavelength
        self.fnujy = t['FLUX'].to(u.jansky,
                                  equivalencies=u.spectral_density(self.wavelength*u.micron)).value

        neg = self.fnujy <= 0.0
        self.fnujy[neg] = cfg.tiny

        self.fill_irradiance()
        return self

    @classmethod
    def vega_rieke(cls, file='rieke08-vega-sun.txt', file_loc=None):
        """Return the "Vega" spectrum from Rieke et al. (2008)
            
        This spectrum has slightly lower fluxes than stated in the
        paper (7.14 vs. 7.15Jy @ 23.675um), and the adopted calibration
        has 7.17Jy, so to ensure that the calibration is as precise as
        possible multiply the spectrum up by 7.17/7.14. A possible
        reason for the discrepancy is the F_lambda to F_nu conversion,
        as there is also a small difference in the 10.6um values in
        Table 1 of that paper.

        Parameters
        ----------
        file : str, optional
            Name of the Vega spectrum file.
        file_loc : str, optional
            Location of the file.
        """
        
        self = cls()
        if file_loc is None:
            file_loc = os.path.dirname(os.path.abspath(__file__)) + '/data/calibration/'
    
        t = Table.read(file_loc+file, format='ascii.cds')
        wave_um = t['Wave'].to(u.micron).value
        _, srt = np.unique(wave_um, return_index=True)
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

    @property
    def e_fnujy(self):
        return self._e_fnujy

    @e_fnujy.setter
    def e_fnujy(self, value):
        if self.fnujy is None:
            expected_len = None
        else:
            expected_len = len(self.fnujy)
        self._e_fnujy = utils.validate_1d(value, expected_len)

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
    def fnujy(self, value):
        if self.nu_hz is None:
            expected_len = None
        else:
            expected_len = len(self.nu_hz)
        self._fnujy = utils.validate_1d(value, expected_len)

                
class ModelSpectrum(Spectrum):
    """Class used for model spectra

    Flux density is per steradian
    """
    
    def __init__(self, name=None, parameters=None, param_values=None,
                 wavelength=None, nu_hz=None, fnujy_sr=None,
                 irradiance_sr=None):
        self.name = name
        self.parameters = parameters
        self.param_values = param_values
        self.wavelength = wavelength
        self.nu_hz = nu_hz
        self.fnujy_sr = fnujy_sr
        self.irradiance_sr = irradiance_sr

        # fill wavelength/frequency if possible
        if self.wavelength is None and self.nu_hz is not None:
            self.fill_hz2wave()
        if self.nu_hz is None and self.wavelength is not None:
            self.fill_wave2hz()

    @classmethod
    def bnu_wave_micron(cls, wave_um, temp, lam0=None, beta=None):
        """Return a (modified) black body spectrum"""

        self = cls()
        
        self.name = 'bb'
        self.parameters = ['log_Temp']
        self.param_values = [temp]
        fnujysr = utils.bnu_wav_micron(wave_um, temp)
        srt = np.flipud(np.argsort(wave_um))
        self.wavelength = wave_um[srt]
        self.nu_hz = c_micron / wave_um[srt]
        self.fnujy_sr = fnujysr[srt]
        
        # modified if desired
        if lam0 is not None and beta is not None:
            lam0 = float(lam0)
            self.name = 'modbb'
            self.parameters = ['log_Temp', 'log_lam0', 'beta']
            self.param_values = [temp, lam0, beta]
            ok = self.wavelength > lam0
            self.fnujy_sr[ok] /= np.power(self.wavelength[ok]/lam0, beta)
        
        return self

    @classmethod
    def bnu_nu_hz(cls, nu_hz, temp, lam0=None, beta=None):
        """Return a black body spectrum
        """

        self = cls()
        wave_um = c_micron / nu_hz
        return self.bnu_wave_micron(wave_um, temp, lam0=lam0, beta=beta)

    @classmethod
    def read_phoenix(cls, file):
        """Return the spectrum from a PHOENIX model
        
        Format is wavelength (A), flam (erg/s/cm^2/cm), bnu 
        (erg/s/cm^2/cm), and a bunch of other stuff we don't care 
        about. The flam column is log10(flam). The shortest wavelength
        is 10A (0.001um), so short of what we'll need here (~0.1um). 

        A factor of pi is needed to get the PHOENIX models into Jy/sr
        (and get stellar radii correct). The docs
        https://phoenix.ens-lyon.fr/Grids/FORMAT suggest that
        (Rstar/d)^2 is the normalisation, so divide spectra by pi here.
        Comparing the results with Kurucz models shows this is correct.

        .. todo:: Can read go faster? BT-Settl files are epic
        
        """

        self = cls()

        self.name = 'phoenix'
        self.parameters = ['Teff', 'logg', 'MH']
        if 'BT-Settl' in file:
            # first for 'normal' phoenix files, second for compatibility
            # with files without alpha enh and decimal temperatures
            par = re.search("lte([0-9]+)-([0-9.]+)([+-][0-9.]+)"
                            "[a]([+-][0-9.]+).BT-Settl.7.bz2", file)
            if par is None:
                par = re.search("lte([0-9.]+)-([0-9.]+)([+-][0-9.]+)"
                                ".BT-Settl.7.bz2", file)
            if par is None:
                raise utils.SdfError("couldn't re file {}".format(file))
            teff = float(par.groups()[0])*100.0
            logg = float(par.groups()[1])
            mh = float(par.groups()[2])
        self.param_values = {'Teff': teff, 'logg': logg, 'MH': mh}

        # attempt to load numpy save from last time, otherwise read in
        try:
            w, f = np.load(file+'.npy')
        except FileNotFoundError:
            c = lambda s: s.decode().replace("D", "E")  # deal with fortran 1.0D+01
            w, f = np.loadtxt(file, dtype=float, usecols=(0, 1),
                              converters={0: c, 1: c}, unpack=True)
            np.save(file, [w, f])
        
        f = np.power(10, f)
        _, srt = np.unique(w, return_index=True)
        srt = np.flipud(srt)
        wav = wave_vac2air(w[srt])
        wav = u.Quantity(wav, u.Angstrom)
        flam = u.Quantity(f[srt]/np.pi, u.erg/u.s/u.cm**2/u.cm)  # divide by pi
        self.wavelength = wav.to(u.micron).value
        self.nu_hz = wav.to('Hz', equivalencies=u.spectral()).value
        self.fnujy_sr = flam.to(u.jansky,
                                equivalencies=u.spectral_density(wav)).value
        self.extend_power_law()
        # make npt large, as we may resample with R>1000
        self.fill_gaps(npt=5000)
        self.sort('wave')
        return self

    @classmethod
    def read_koester(cls, file):
        """Return the spectrum from a Koester WD model

        Grid was obtained from VOSA, assume wavelengths already in air.
        Format is wavelength (A), flam (erg/s/cm^2/A). The shortest
        wavelength is about 900A (0.09um) and the longest 30, 000A (30um).
        """

        self = cls()

        self.name = 'koester_wd'
        self.parameters = ['Teff', 'logg']

        par = re.search("da([0-9]+)_([0-9]+).dk.dat.txt", file)
        if par is None:
            raise utils.SdfError("couldn't re file {}".format(file))
        teff = float(par.groups()[0])
        logg = float(par.groups()[1])/100
        self.param_values = {'Teff': teff, 'logg': logg}

        # read in
        w, f = np.loadtxt(file, unpack=True)

        # add short wavelength points
        w = np.append(w, [np.min(w)-np.median(np.diff(w)),
                          cfg.models['min_wav_micron']*1e4])
        f = np.append(f, [cfg.tiny, cfg.tiny])

        _, srt = np.unique(w, return_index=True)
        srt = np.flipud(srt)
        wav = u.Quantity(w[srt], u.Angstrom)
        flam = u.Quantity(f[srt]/np.pi, u.erg/u.s/u.cm**2/u.Angstrom)  # divide by pi
        self.wavelength = wav.to(u.micron).value
        self.nu_hz = wav.to('Hz', equivalencies=u.spectral()).value
        self.fnujy_sr = flam.to(u.jansky,
                                equivalencies=u.spectral_density(wav)).value
        self.extend_power_law()
        # make npt large, as we may resample with R>1000
        self.fill_gaps(npt=5000)
        self.sort('wave')
        return self

    def read_kurucz(file):
        """Read a grid of Castelli & Kurucz models from a specific file.
        Returns a tuple with (teff, logg, [M/H], models), where the first
        three are arrays with the parameters of each model, and models
        is a list of ModelSpectrum objects holding the spectra.
        
        .. todo:: not really clear whether this belongs here or in 
        SpecModel (or elsewhere). It's only here because it's similar
        to read_phoenix, which deals with single files.
        """

        fh = open(file, 'r')
        lines = fh.readlines()

        cols = ((0, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60), (60, 70), (70, 80))

        # create list of models, and arrays of their parameters
        teff = np.array([], dtype=float)
        logg = np.array([], dtype=float)
        mh = np.array([], dtype=float)
        models = np.array([], dtype=object)
        wave = np.array([], dtype=float)
        model = []
        for i, l in enumerate(lines):
            l = l.rstrip()

            # wavelength grid comes first, starting on a line with only numbers
            if re.search('[a-zA-Z]', l) is None and len(teff) == 0 and 'TEFF' not in l:
                for c in cols:
                    if c[0] < len(l):
                        wave = np.append(wave, float(l[c[0]:c[1]]))
                continue

            # skip to next line if wave hasn't been started yet
            if len(wave) == 0:
                continue
            
            # new models have TEFF on first line
            if 'TEFF' in l:
                teffLogg = re.search(r'TEFF\s+(\d+\.\d*)\s+GRAVITY\s+(\d+\.\d*)', l)
                teff = np.append(teff, float(teffLogg.groups()[0]))
                logg = np.append(logg, float(teffLogg.groups()[1]))
                metalAlpha = re.search(r'\[([\s+-]?\d+\.\d+)([aAbB]?)\]?', l)
                mh = np.append(mh, float(metalAlpha.groups()[0]))
                if metalAlpha.groups()[1] != '':
                    raise utils.SdfError("Alpha non-nothing ({}), wtf?".
                                         format(metalAlpha.groups()))
            else:
                for c in cols:
                    if c[0] < len(l):
                        model.append(float(l[c[0]:c[1]]))
                
            # add complete model to list
            if len(model) > 0 and ('TEFF' in l or i == len(lines)-1):
                self = ModelSpectrum()
                self.parameters = ['Teff', 'logg', 'MH']
                # if next model has been started we want
                # params from last model
                if 'TEFF' in l:
                    self.param_values = {'Teff': teff[-2],
                                         'logg': logg[-2],
                                         'MH': mh[-2]}
                else:
                    self.param_values = {'Teff': teff[-1],
                                         'logg': logg[-1],
                                         'MH': mh[-1]}
            
                # metalliticity
                mhstr = str(self.param_values['MH'])
                if file[1] == 'p':
                    mhstr = '+' + mhstr
                elif self.param_values['MH'] == 0.0:
                    mhstr = '-' + mhstr
                self.name = 'kurucz' + mhstr
                
                self.wavelength = (wave * u.nm).to('micron').value
                self.nu_hz = c_micron / self.wavelength
                model = model[0:len(wave)]
                self.fnujy_sr = np.array(model, dtype=float) * 4e23
                # turn zeros into small numbers to avoid divide errors
                small = (self.fnujy_sr <= 0.0)
                self.fnujy_sr[small] = cfg.tiny
                # extend to long wavelengths, fill out, and sort again
                self.extend_power_law()
                self.fill_gaps()
                self.sort('wave')
                models = np.append(models, self)
                model = []

        return models, teff, logg, mh

    @property
    def fnujy_sr(self):
        return self._fnujy_sr

    @fnujy_sr.setter
    def fnujy_sr(self, value):
        if self.nu_hz is None:
            expected_len = None
        else:
            expected_len = len(self.nu_hz)
        self._fnujy_sr = utils.validate_1d(value, expected_len)


def wave_vac2air(wavelength):
    """Convert wavelengths in vacuum to air."""
    div = 1 + 1e-6*nrefrac(wavelength)
    return wavelength / div


def nrefrac(wl, density=1.0):
    """Calculate refractive index of air from Cauchy formula.

    Note that Phoenix delivers synthetic spectra in the vaccum and
    that a line shift is necessary to adapt these synthetic spectra
    for comparisons to observations from the ground. For this, divide
    the vacuum wavelengths by (1+1.e-6*nrefrac) as returned from the
    function below to get the air wavelengths (or use the equation
    for AIR from it).
    https://osubdd.ens-lyon.fr/phoenix/doc/spectra.html

    Input: wavelength in Angstrom, density of air in amagat (relative to STP,
    e.g. ~10% decrease per 1000m above sea level).
    Returns N = (n-1) * 1.e6.
    """

    # The IAU standard for conversion from air to vacuum wavelengths is given
    # in Morton (1991, ApJS, 77, 119). For vacuum wavelengths (VAC) in
    # Angstroms, convert to air wavelength (AIR) via:

    #  AIR = VAC / (1.0 + 2.735182E-4 + 131.4182 / VAC^2 + 2.76249E8 / VAC^4)

    wl2inv = (1.e4/wl)**2
    refracstp = 272.643 + 1.2288 * wl2inv + 3.555e-2 * wl2inv**2
    return density * refracstp
