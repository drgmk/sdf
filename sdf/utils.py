from functools import lru_cache
from contextlib import contextmanager
import os
from hashlib import sha1

from scipy import sparse
import numpy as np
from scipy.integrate import simps
from scipy.ndimage import spline_filter
import astropy.units as u

from . import photometry
from . import spectrum
from . import config as cfg

c_micron = u.micron.to(u.Hz,equivalencies=u.spectral())


@contextmanager
def pushd(new_dir):
    """A context manager that implements the `pushd` command.
        
    This lets you run a block of commands while in a different 
    directory.
        
    From https://gist.github.com/theY4Kman/3583442
    """
    old_dir = os.getcwd()
    os.chdir(new_dir)
    try:
        yield old_dir
    finally:
        os.chdir(old_dir)


def bnu_wav_micron(wav_um,temp):
    """Return a Planck function, avoiding overflows
    """
    k1 = 3.9728949e19
    k2 = 14387.69
    fact1 = k1/(wav_um**3)
    fact2 = k2/(wav_um*temp)
    if isinstance(wav_um,np.ndarray):
        ofl = fact2 < 709
        bnu = np.zeros(len(wav_um)) + cfg.tiny
        if np.any(ofl) == False:
            return bnu
        else:
            bnu[ofl] = fact1[ofl]/(np.exp(fact2[ofl])-1.0)
            return bnu
    else:
        if fact2 > 709:
            return cfg.tiny
        else:
            return fact1/(np.exp(fact2)-1.0)

        
def bnu_nu_hz(nu_hz,temp):
    wav_um = c_micron / nu_hz
    return bnu_wav_micron(wav_um,temp)

def sdf_int(y,x):
    """Decide how we do integration in sdf"""
    return np.trapz(y,x)
#    return simps(y,x)

def validate_1d(value,expected_len,dtype=float):
    
    if type(value) in [list, tuple]:
        value = np.array(value,dtype=dtype)
        
    if value is None:
        value = value
    elif isinstance(value, np.ndarray) and value.ndim == 1:
        if expected_len is not None:
            if len(value) != expected_len:
                raise SdfError("incorrect length (expected {0} but found {1})".format(expected_len, len(value)))
        if value.dtype != dtype:
                value = np.array(value,dtype=dtype)
    else:
        raise TypeError("should be a 1-d sequence")
    
    return value

def validate_nd(value,expected_dim,dtype=float):
    
    if value is None:
        value = value
    elif isinstance(value, np.ndarray):
        if expected_dim is not None:
            if value.ndim != expected_dim:
                raise SdfError("incorrect dimension (expected {0} but found {1})".format(expected_dim, value.ndim))
        if value.dtype != dtype:
            value = np.array(value,dtype=dtype)
    else:
        raise TypeError("should be an n-d ndarray")
    
    return value

def validate_string(value):
    if value is None:
        return value
    elif isinstance(value, str):
        return value
    else:
        raise TypeError("should be a string")

def validate_dict(value):
    if value is None:
        return value
    elif isinstance(value, dict):
        return value
    else:
        raise TypeError("should be a dict")

def validate_float(value):
    if value is None:
        return value
    elif isinstance(value, float):
        return value
    else:
        raise TypeError("{} should be a float".format(value))

def validate_function(value):
    if value is None:
        return value
    elif callable(value):
        return value
    else:
        raise TypeError("should be a function")

def validate_quantity(value):
    if value is None:
        return value
    elif isinstance(value, u.Quantity):
        return value
    else:
        raise TypeError("should be an astropy units Quantity")

class SdfError(Exception):
    """Use this for non-standard errors
    """
    pass


def resample_matrix(wave_in,new_wave,new_R,old_R=np.inf,
                    kern_width=5):
    """Return a resampling/convolution kernel.
        
    Copied from Casey's sick, modified to ensure the range of indices
    included near each resampled wavelength are appropriate (rather
    than a single value for the whole spectrum).
    """

    # this will be the shape of the convolution matrix
    N, M = (new_wave.size, wave_in.size)

    # width of the kernel (dlambda=lambda/R) at each point
    fwhms = new_wave / float(new_R)
    if np.isfinite(old_R):
        assert old_R > new_R
        fwhms -= new_wave / float(old_R)

    # 2.355 ~= 2 * sqrt(2*log(2))
    sigmas = fwhms/2.3548200450309493

    # approx center indices for new wavs in old
    ioff = wave_in.searchsorted(new_wave)
    ioff = np.clip(ioff,0,M-1)

    # index ranges covered by convolution at a given point
    dlambda = np.diff(wave_in)
    dlambda = np.append(dlambda,dlambda[-1])
    dlambda = dlambda[ioff]
    nkern = np.ceil(kern_width * sigmas / dlambda).astype(int)

    # For +/- N_kernel_pixels at each point, calculate the kernel
    # and retain the indices.
    x_indices = np.array([],dtype=int)
    pdf = np.array([],dtype=int)
    for i in range(N):
        xtmp = np.arange(ioff[i] - nkern[i],ioff[i] + nkern[i],1)
        xtmp = np.clip(xtmp,0,M-1)
        x_indices = np.append(x_indices,xtmp)
        pdftmp = np.exp(-(wave_in[xtmp] - new_wave[i])**2/(2.*sigmas[i]**2))
        pdftmp /= pdftmp.sum()
        pdf = np.append(pdf,pdftmp)

    y_indices = np.repeat(np.arange(N), 2 * nkern)

    return sparse.coo_matrix((pdf, (x_indices, y_indices)), shape=(M, N))


def get_observations(file):
    """Get observations from a sdb file."""

    # first photometry
    p = photometry.Photometry.read_sdb_file(file) # Photometry object
    # then a tuple of ObsSpectrum objects
    s = spectrum.ObsSpectrum.read_sdb_file(file,module_split=True)
    if s is None:
        obs = (p,)
    else:
        obs = (p,) + s

    return obs


def uvby_convert(by,m1,c1):
    """Convert Stromgren photometry according to Bessell 2011.
        
    The coefficients are to convert synthetic photometry to observed,
    I_std = a_0 + a_1 * I_syn, so need to be inverted as we want to
    convert observed photometry to agree with synthetic.
    
    """
    
    # coeffecients from Table 2
    by_1 = [-0.007,0.997]
    by_2 = [ 0.004,0.979]
    m1_1 = [ 0.005,0.963]
    m1_2 = [ 0.011,0.951]
    c1_1 = [-0.016,0.994]
    c1_2 = [-0.003,1.018]

    if by < 0.5:
        by_out = (by - by_1[0])/by_1[1]
        m1_out = (m1 - m1_1[0])/m1_1[1]
        c1_out = (c1 - m1_1[0])/c1_1[1]
    else:
        by_out = (by - by_2[0])/by_2[1]
        m1_out = (m1 - m1_2[0])/m1_2[1]
        c1_out = (c1 - m1_2[0])/c1_2[1]

    return by_out,m1_out,c1_out


def linterp(newx,x,y):
    """Linear interpolation."""

    hi = x.searchsorted(newx)
    lo = hi - 1
    return y[lo] + (y[hi]-y[lo])*(newx-x[lo])/(x[hi]-x[lo])


@lru_cache(maxsize=16)
def spline_filter_mem(arr,order=None):
    """Filter array, memoizing result.
        
    Is passed an array using the hashable wrapper"""
    
    return spline_filter(arr.unwrap(),order=order)


class hashable(object):
    r'''Hashable wrapper for ndarray objects.
        
    Instances of ndarray are not hashable, meaning they cannot be added to
    sets, nor used as keys in dictionaries. This is by design - ndarray
    objects are mutable, and therefore cannot reliably implement the
    __hash__() method.
    
    The hashable class allows a way around this limitation. It implements
    the required methods for hashable objects in terms of an encapsulated
    ndarray object. This can be either a copied instance (which is safer)
    or the original object (which requires the user to be careful enough
    not to modify it).
        
    From https://gist.github.com/marquisthunder/9ef974477e1aef9dbb41
    '''
    def __init__(self, wrapped, tight=False):
        r'''Creates a new hashable object encapsulating an ndarray.
            wrapped
                The wrapped ndarray.
            tight
                Optional. If True, a copy of the input ndaray is created.
                Defaults to False.
        '''
        self.__tight = tight
        self.__wrapped = np.array(wrapped) if tight else wrapped
        self.__hash = int(sha1(wrapped.view(np.uint8)).hexdigest(), 16)

    def __eq__(self, other):
        return np.all(self.__wrapped == other.__wrapped)

    def __hash__(self):
        return self.__hash

    def unwrap(self):
        r'''Returns the encapsulated ndarray.
            If the wrapper is "tight", a copy of the encapsulated ndarray is
            returned. Otherwise, the encapsulated ndarray itself is returned.
        '''
        if self.__tight:
            return np.array(self.__wrapped)

        return self.__wrapped
