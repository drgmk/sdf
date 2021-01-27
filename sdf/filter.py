from functools import lru_cache
import warnings
import numpy as np
from astropy.utils.data import download_file
from astropy.io.votable import parse,exceptions
import astropy.units as u
from . import spectrum
from . import filter_info
from . import utils

# ignore warnings about SVO filter votables
warnings.simplefilter('ignore', exceptions.W42)

c_micron = u.micron.to(u.Hz,equivalencies=u.spectral())

def common_x(x1,y1,x2,y2):
    """Interpolate two y arrays to a single common x array"""
    if len(x1) == len(x2):
        if np.all( np.equal(x1,x2) ):
            _,srt = np.unique(x1,return_index=True)
            return (x1[srt],y1[srt],y2[srt])
    xall = np.append(x1,x2)
    xmin = np.max([np.min(x1),np.min(x2)])
    xmax = np.min([np.max(x1),np.max(x2)])
    keep = (xall >= xmin) & ( xall <= xmax)
    if np.any(keep) == False:
        raise utils.SdfError("x arrays don't overlap {} to {} and {} to {}".format(
                         np.min(x1),np.max(x1),np.min(x2),np.max(x2)))
    x = np.unique( xall[keep] )
    srt = np.argsort(x1)
    y1 = np.interp(x,x1[srt],y1[srt])
    srt = np.argsort(x2)
    y2 = np.interp(x,x2[srt],y2[srt])
    return (x,y1,y2)
    

class Filter(object):
    """Filter class to use for synthetic photometry

    A good primer on filter responses is Appendix A 2012PASP..124..140B.
    The key point is that photon-counting is equivalent to 
    energy-integration (A13), so synthetic photometry can always be done
    with the same equation if the filter responses are changed 
    accordingly. A quantum-efficiency-based response is S, and 
    S'=lambda.S is the energy-integrating equivalent, which can be
    integrated directly (F_lambda.S'.dlambda, A11). For 2MASS 
    multiplication by lambda was done already (2003AJ....126.1090C), but
    this is not always the case so must be checked for each filter, and
    photon-counting responses converted to energy-integrating.

    We also need to know the reference spectrum to do synthetic 
    photometry in systems where the flux density is given at a specific
    (e.g. mean) wavelength. This is commonly the case for IR photometry,
    but not in the optical (e.g. see appendix of 2011AJ....141..173B).
    It is also good to know how that reference wavelength was derived,
    but it can also be taken as a number. This correction converts the
    measurement (also the value that the reference spectrum would have
    at this wavelength) to the value the true spectrum must have at this
    wavelength in order to reproduce the observed signal. The most
    common reference spectrum is nu.F_nu=const, but others are also used
    (e.g. 10,000K blackbody for MIPS).
    
    Filter information comes from the filter_info.py file, all filters,
    colours, indices etc. need to be named as keys to the filters 
    dictionary that comes from compiliing that file, otherwise they will
    not be computed in the models.
    
    A zero point offset of zero is set by default.
    
    """

    # this needs to return all desired filters, colours, and indices
    all = list(filter_info.filters.keys())
    
    def __init__(self,name=None,system=None,nu_hz=None,fileloc=None,
                 response=None,response_type=None,
                 zero_point=None,zero_point_offset=None,
                 measurement_calibration=None,
                 magnitude_system=None,ref_wavelength=None,ref_nu_hz=None,
                 ref_spectrum=None,cc_denom=None,Af_Av=None):
        self.name = name
        self.nu_hz = nu_hz
        self.response = response
        self.response_type = response_type
        self.magnitude_system = magnitude_system
        self.zero_point = zero_point
        self.zero_point_offset = zero_point_offset
        self.measurement_calibration = measurement_calibration
        self.ref_wavelength = ref_wavelength
        self.ref_nu_hz = ref_nu_hz
        self.ref_spectrum = ref_spectrum
        self.cc_denom = cc_denom
        self.Af_Av = Af_Av


    @classmethod
    def svo_get(cls,name):
        """Get response and details for a filter from Spanish VO
            
        Other properties still need to be set "by hand", as this
        service doesn't provide either the reference spectrum or
        wavelength. It doesn't appear that their detector type is
        anything other than energy counting, so any conversion
        between the response the and "RSR" [lambda.R(lambda)] also
        needs to be figured out manually. All this is done in
        filter_info.py
        """
        
        self = cls()
        
        # get the filter, either from url or cache, if we want something
        # other than Vega use the PhotCalID in the name given
        if 'PhotCalID' in name:
            url = "http://svo2.cab.inta-csic.es/theory/fps3/fps.php?"+name
        else:
            url = "http://svo2.cab.inta-csic.es/theory/fps3/fps.php?ID="+name
        fileloc = download_file(url,cache=True)
        self.fileloc = fileloc
        
        # open the file and get the necessary info
        votable = parse(fileloc)
        
        # grab the filter name
        name_field = votable.get_field_by_id('filterID')
        self.name = name_field.value#.decode()
        
        # mean wavelength
        wmean = votable.get_field_by_id('WavelengthMean')
        self.ref_wavelength = wmean.value * wmean.unit
        
        # zero point in Jy, and offset in mag
        zp = votable.get_field_by_id('ZeroPoint')
        self.zero_point = (zp.value * zp.unit).to('Jy').value
        self.zero_point_offset = 0.0
        
        # magnitude system
        sys = votable.get_field_by_id('MagSys')
        self.magnitude_system = sys.value#.decode()
        
        # and the response, assuming there is no masking [hence filled()]
        vo_filt_table = votable.get_first_table()
        filt_table = vo_filt_table.to_table().filled()
        wav = filt_table['Wavelength']
        self.nu_hz = np.array( wav.to('Hz',equivalencies=u.spectral()) )
        self.response = np.array( filt_table['Transmission'] )
        return self


    @lru_cache(maxsize=128)
    def get_memo(name):
        """Get response and details for a filter (via filter_info.py)
            
        This function is memoized for speed, the maxsize just needs
        to be larger than the number of filters.
        """

        f = filter_info.filters

        # create (relatively narrow) generic filters
        if name[0:3] == 'WAV':
            self = Filter()
            self.name = name
            self.response_type = 'energy'
            cwav = float(name[3:])
            wave = np.arange( cwav*0.95, cwav*1.05, cwav/100. )
            self.nu_hz = c_micron / wave
            self.response = np.ones( len(wave) )
        elif name not in f:
            raise KeyError("Filter {} not in filter_info".format(name))
        else:
            if 'svo_name' in f[name]:
                self = Filter.svo_get(f[name]['svo_name'])
            else:
                self = Filter()
                if f[name]['wav_micron'] is not None:
                    self.nu_hz = c_micron / np.array(f[name]['wav_micron'])
                if f[name]['response'] is not None:
                    self.response = np.array(f[name]['response'])

            self.name = name
        
            # fill fields from filter_info, these override SVO
            if 'magnitude_system' in f[name]:
                if f[name]['magnitude_system'] is not None:
                    self.magnitude_system = f[name]['magnitude_system']
            if 'zero_point' in f[name]:
                if f[name]['zero_point'] is not None:
                    self.zero_point = f[name]['zero_point']
            if 'zero_point_offset' in f[name]:
                if f[name]['zero_point_offset'] is not None:
                    self.zero_point_offset = f[name]['zero_point_offset']
            if 'measurement_calibration' in f[name]:
                if f[name]['measurement_calibration'] is not None:
                    self.measurement_calibration = f[name]['measurement_calibration']
            if 'ref_wavelength' in f[name]:
                if f[name]['ref_wavelength'] is not None:
                    self.ref_wavelength = f[name]['ref_wavelength']
            if 'ref_spectrum' in f[name]:
                if f[name]['ref_spectrum'] is not None:
                    self.ref_spectrum = f[name]['ref_spectrum']
            if 'response_type' in f[name]:
                if f[name]['response_type'] is not None:
                    self.response_type = f[name]['response_type']
            
        # convert photon counting responses to energy
        if 'photon' in self.response_type:
            if self.nu_hz is not None and self.response is not None:
                self.response /= self.nu_hz
                self.response_type = 'photon -> energy'

        # sort, normalise, and fill
        if self.nu_hz is not None:
            self.sort()
            self.normalise_response()
            self.fill_mean_wavelength()
            self.fill_ref_nu_hz()
            self.fill_cc_denom()
        
        # compute zero-point in Vega system
        if self.magnitude_system == 'Vega':
            v = spectrum.ObsSpectrum.vega_stis()
            if (np.min(v.nu_hz) < np.min(self.nu_hz) and
                np.max(v.nu_hz) > np.min(self.nu_hz)):
                zp = self.synthphot(v)
                self.zero_point = zp[0]
            else:
                self.zero_point = None

        # or for AB
        elif self.magnitude_system == 'AB':
            self.zero_point = 3631.0

        return self


    @lru_cache(maxsize=128)
    def get_with_zpo(name,zpo_file='/Users/grant/astro/projects/sdf/sdf/calibration/zpos.txt'):
        """Get filter, using ZPO from file."""
        
        f = Filter.get_memo(name)
        t = np.genfromtxt(zpo_file,dtype=None)
        fs = np.array([filt.decode() for filt in t[0]])
        if name in fs:
            f.zero_point_offset = float( t[1][ fs == name ][0] )

        return f
    
    
    # decide what get() points to
    get = get_memo
    

    def mag2flux(self,mag):
        """Convert magnitudes to flux density.

        Use zero point associated with this filter to convert a given
        magnitude to flux density. Units returned are those associated
        with the zero point.
        """

        if self.zero_point is None or self.zero_point_offset is None:
            raise utils.SdfError("no zero point or offset for filter {})".
                           format(self.name))
        return self.zero_point * 10**(-0.4*(mag-self.zero_point_offset))
    
    
    def flux2mag(self,flux):
        """Convert flux density to magnitudes.
            
        Use zero point associated with this filter to convert a given
        flux density to magnitude. Units returned are those associated
        with the zero point.
        """
        
        if self.zero_point is None or self.zero_point_offset is None:
            raise utils.SdfError("no zero point or offset for filter {})".\
                           format(self.name))
        return self.zero_point_offset                               \
               - 2.5 * np.log10( flux / self.zero_point )


    def measflux2flux(self,flux):
        """Convert a measured flux to an actual flux.
        
        Use flux calibration function given for this filter.
        """

        if self.measurement_calibration is not None:
            return self.measurement_calibration(flux)
        else:
            return flux


    def sort(self):
        """Sort response in increasing order of frequency."""
        
        _,srt = np.unique(self.nu_hz,return_index=True)
        self.nu_hz = self.nu_hz[srt]
        self.response = self.response[srt]


    def normalise_response(self):
        """Normalise filter response so integral is one."""
        
        norm = utils.sdf_int(self.response,self.nu_hz)
        self.response /= norm

    def actual_flux(self,spectrum):
        """Return the spectrum's flux at the reference wavelength.
            
        Interpolate in log space.
        """
        
        spectrum.sort('nu')
        if hasattr(spectrum,'fnujy'):
            log_fnu = np.log10(spectrum.fnujy)
        elif hasattr(spectrum,'fnujy_sr'):
            log_fnu = np.log10(spectrum.fnujy_sr)
        
        log_nu = np.log10(spectrum.nu_hz)
        log_ref_fnu = np.interp(np.log10(self.ref_nu_hz),log_nu,log_fnu)
        return np.power(10,log_ref_fnu)


    def fill_ref_nu_hz(self):
        """Set the reference frequency of the filter."""
        
        if self.ref_wavelength is not None:
            self.ref_nu_hz = c_micron / self.ref_wavelength


    def fill_mean_wavelength(self):
        """Set the mean wavelength of the filter."""
        
        wave = c_micron / self.nu_hz
        self.mean_wavelength = utils.sdf_int(self.response,wave)            \
                               / utils.sdf_int(self.response/wave,wave)


    def pivot_wavelength(self):
        """Return the pivot wavelength.
            
        See Bessell & Murphy (2012) for definition. This wavelength is a
        property of the filter only, and results in exact conversion
        between mean F_lambda and mean F_nu.
        """
    
        wave = c_micron / self.nu_hz
        fp = utils.sdf_int( self.response, wave )                           \
             / utils.sdf_int( self.response / wave / wave, wave )
        return np.sqrt(fp)


    def fill_cc_denom(self):
        """Compute the denominator of the colour correction
        
        This calculation is specific to the filter response (i.e. 
        independent of observed spectrum) so only needs to be calculated
        once. See synthphot for details.
        """

        if self.ref_spectrum is not None:
            d = utils.sdf_int(self.response * self.ref_spectrum(self.nu_hz),
                              self.nu_hz)
            self.cc_denom = d / self.ref_spectrum(self.ref_nu_hz)


    def synthphot(self,spectrum):
        """Synthetic photometry of supplied spectrum object

        In the IR we have f_nu(quoted) (i.e. catalogue value) is:

        f_nu(quoted) = f_nu(lam_eff) K

        where the colour correction K is:

        K =        1        S_nu(lam_eff)   int( f_nu R dnu )
            ------------- -----------------
            f_nu(lam_eff) int( S_nu R dnu )

        For a correctly normalised bandpass R (energy counting), the
        integral in the numerator is "normal" synthetic photoemtry as
        done in the optical (e.g. appendix of 2011AJ....141..173B).

        Method is to separate the terms in synthetic photometry so that
        it can be done with a single set of equations. Separation is 
        into
          1 the "normal" integration as done for the optical
          2 the spectrum-independent part of the colour correction 
            (second term, which can be pre-computed)
        The latter is only done if necessary. f_nu(lam_eff) cancels so
        isn't needed to get f_nu(quoted), but is needed if we want to
        know what the colour correction was.
        
        Returns a tuple of (quoted flux, colour correction), where the
        colour correction will be None if there isn't one (i.e. no ref
        spectrum).
        
        For colours and indices the call needs to be done using the
        synthphot method for Spectrum objects.
        
        """

        if hasattr(spectrum,'fnujy'):
            fnu = spectrum.fnujy
        elif hasattr(spectrum,'fnujy_sr'):
            fnu = spectrum.fnujy_sr
        x,y1,y2 = common_x(self.nu_hz,self.response,spectrum.nu_hz,fnu)

        # this is the basic synthetic photometry
        cc_num = utils.sdf_int(y1 * y2, x)
        
        # if no reference spectrum we just return the integral because
        # we normalised the bandpass when it was loaded
        if self.ref_spectrum is None:
            return (cc_num,None)
        else:
            fnu_eff = self.actual_flux(spectrum)
            cc = cc_num / self.cc_denom / fnu_eff
            return (cc*fnu_eff,cc)


    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, value):
        self._name = utils.validate_string(value)

    @property
    def magnitude_system(self):
        return self._magnitude_system
    @magnitude_system.setter
    def magnitude_system(self, value):
        self._magnitude_system = utils.validate_string(value)
        
    @property
    def zero_point(self):
        return self._zero_point
    @zero_point.setter
    def zero_point(self, value):
        self._zero_point = utils.validate_float(value)

    @property
    def zero_point_offset(self):
        return self._zero_point_offset
    @zero_point_offset.setter
    def zero_point_offset(self, value):
        self._zero_point_offset = utils.validate_float(value)
    
    @property
    def measurement_calibration(self):
        return self._measurement_calibration
    @measurement_calibration.setter
    def measurement_calibration(self, value):
        self._measurement_calibration = utils.validate_function(value)

    @property
    def ref_spectrum(self):
        return self._ref_spectrum
    @ref_spectrum.setter
    def ref_spectrum(self, value):
        self._ref_spectrum = utils.validate_function(value)

    @property
    def response_type(self):
        return self._response_type
    @response_type.setter
    def response_type(self, value):
        self._response_type = utils.validate_string(value)

    @property
    def nu_hz(self):
        return self._nu_hz
    @nu_hz.setter
    def nu_hz(self, value):
        self._nu_hz = utils.validate_1d(value,None)

    @property
    def response(self):
        return self._response
    @response.setter
    def response(self,value):
        # always set nu_hz first, so no check for attribute nu_hz
        if self.nu_hz is None:
                expected_len = None
        else:
                expected_len = len(self.nu_hz)
        self._response = utils.validate_1d(value,expected_len)


def mean_wavelength(filternames):
    """Return mean_wavelengths given tuple/list of filter names."""
    
    wav = np.array([])
    for f in filternames:
        if iscolour(f):
            col = Colour.get(f)
            wav = np.append(wav,col.mean_wavelength)
        else:
            filt = Filter.get(f)
            wav = np.append(wav,filt.mean_wavelength)
    return wav


class Colour(object):
    """Class for colours/indices, largely for convenience."""

    def __init__(self,name=None,filters=None,weights=None,
                 mean_wavelength=None):
        self.name = name
        self.filters = filters
        self.weights = weights
        self.mean_wavelength = mean_wavelength


    @classmethod
    def get(cls,name):
        """Get a Colour object given a name."""
    
        if not iscolour(name):
            raise utils.SdfError("name given ({}) not a colour".format(name))

        self = cls()
        self.name = name
        self.fill_info()
        return self
        
        
    def fill_info(self):
        """Set filters, weights, and other info for this colour."""

        # colour
        if '_' in self.name:
            self.filters = np.array(self.name.split('_'),dtype=str)
            self.weights = np.array([1.,-1],dtype=float)
                
        # Stromgren M1, (v-b)-(b-y) or (v - 2b + y)
        elif self.name == 'STROMM1':
            self.filters = np.array(['VS','BS','YS'],dtype=str)
            self.weights = np.array([1.,-2.,1.],dtype=float)

        # Stromgren C1, (u-v)-(v-b) or (u - 2v + b)
        elif self.name == 'STROMC1':
            self.filters = np.array(['US','VS','BS'],dtype=str)
            self.weights = np.array([1.,-2.,1.],dtype=float)

        meanw = np.array([])
        for f in self.filters:
            filt = Filter.get(f)
            meanw = np.append(meanw,filt.mean_wavelength)
        self.mean_wavelength = np.mean(meanw)


@lru_cache(maxsize=128)
def iscolour(filter):
    """Return True if the given filter name is a colour
        
    In practise this means either the name has a '_' in it,
    e.g. BS_YS, indicating BS-YS, or that it's one of the
    Stromgren indices M1 or C1.
    
    If passed an array, list, or tuple, then return a list of
    booleans.
    """
    
    if isinstance(filter,(tuple,list,np.ndarray)):
        iscol = np.array([],dtype=bool)
        for f in filter:
            iscol = np.append(iscol,iscolour(f))
        return iscol
    else:
        if '_' in filter:
            fs = filter.split('_')
            if len(fs) != 2:
                raise utils.SdfError("filter name {} "
                               "has too many '_'s".format(filter))
            if fs[0] == fs[1]:
                raise utils.SdfError("filters in colour {} "
                               " are the same".format(filter))
            return True
        elif filter == 'STROMM1':
            return True
        elif filter == 'STROMC1':
            return True
        else:
            return False
