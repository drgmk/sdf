import glob
from os.path import exists,basename
from itertools import product

import numpy as np
import astropy.units as u

from . import convolve
from . import model
from . import spectrum
from . import filter
from .utils import SdfError
from . import config as cfg

c_micron = u.micron.to(u.Hz,equivalencies=u.spectral())

"""Helper functions to set up models"""

def setup_all():
    """Rederive all models."""
    setup_spec()
    setup_phot()

def setup_spec():
    """Rederive spectrum models."""
    bb_spectra()
    modbb_spectra()
    kurucz_spectra()
    phoenix_spectra()

def setup_phot():
    """Rederive convolved models.
        
    This needs to be run to propagate ZPOs into colours/indices.
    """
    specmodel2phot(mname='kurucz-0.0')
    specmodel2phot(mname='kurucz_m')
    specmodel2phot(mname='phoenix-0.0')
    specmodel2phot(mname='phoenix_m')
    specmodel2phot(mname='bb_disk')
    specmodel2phot(mname='bb_star')
    specmodel2phot(mname='modbb_disk')


def bb_spectra():
    """Generate SpecModel grid of blackbody models."""
    
    model.SpecModel.generate_bb_model(name='bb_disk',
                                      temperatures=cfg.models['bb_disk_temps'],
                                      write=True,overwrite=True)
    model.SpecModel.generate_bb_model(name='bb_star',
                                      temperatures=cfg.models['bb_star_temps'],
                                      write=True,overwrite=True)


def modbb_spectra():
    """Generate SpecModel grid of modified blackbody models."""
    
    model.SpecModel.generate_modbb_model(name='modbb_disk',
                                         temperatures=cfg.models['bb_disk_temps'],
                                         write=True,overwrite=True)


def specmodel2phot(mname='kurucz-0.0',overwrite=False):
    """Generate a PhotModel grid from SpecModel models."""
    
    convolve_specmodel(mname=mname,overwrite=overwrite)
    m = model.PhotModel.read_convolved_models(mname)
    m.write_model(m.name,overwrite=True)


def convolve_specmodel(mname='kurucz-0.0',overwrite=False):
    """Convolve a set of SpecModel models.
    
    Write files containing convolved fluxes for each filter. The source
    spectra are in the SpecModels as these are saved with high enough
    wavelength resolution.

    Range of parameters is the same as the spectra.
    
    """

    m = model.SpecModel.read_model(mname)

    # flat list of all indices (cartesian product)
    i_list = ()
    for i in range(len(m.parameters)):
        i_list += (list(range(m.param_shape()[i])),)
    i_list = [i for i in product(*i_list)]

    # grid with spectral dimension last
    fnujy_sr_roll = np.rollaxis(m.fnujy_sr,axis=0,start=len(m.fnujy_sr.shape))

    # loop through them and create the ConvolvedModels and the files
    filters = filter.Filter.all
    fnujy_sr = np.zeros(len(i_list))
    for fname in filters:
        
        outfile = cfg.model_loc[mname]+fname+'.fits'
        
        if exists(outfile) and overwrite == False:
            print("Skipping {}, file exists".format(fname))
        
        else:
            print("Convolving filter {}".format(fname))
            for i,ind in enumerate(i_list):
                s = spectrum.ModelSpectrum(nu_hz=c_micron/m.wavelength,
                                           fnujy_sr=fnujy_sr_roll[ind])
                conv,cc = s.synthphot(fname)
                fnujy_sr[i] = conv
        
            if len(m.parameters) > 1:
                fnujy_sr_reshaped = np.reshape(fnujy_sr,m.param_shape())
            else:
                fnujy_sr_reshaped = fnujy_sr

            cm = convolve.ConvolvedModel(name=mname,filter=fname,
                                         parameters=m.parameters,
                                         param_values=m.param_values,
                                         fnujy_sr=fnujy_sr_reshaped)

            cm.write_file(outfile,overwrite=overwrite)


def kurucz_spectra():
    """Generate a SpecModel grid of Castelli & Kurucz models.
        
    The models are called 'kurucz'.
    
    """
    
    # ranges, keep all [M/H]
    trange = [3499,26001]
    lrange = [2.9,5.1]
    
    # Solar metallicity
    f00 = cfg.file['kurucz_models']+'fp00k2odfnew.pck'
    m = model.SpecModel.read_kurucz(f00)
    m = model.crop(m,'Teff',trange)
    m = model.crop(m,'logg',lrange)
    m.write_model('kurucz-0.0',overwrite=True)

    # range of metallicities
    m = model.append_parameter(m,'MH',0.0)
    fs = glob.glob(cfg.file['kurucz_models']+'f*k2odfnew.pck')
    for f in fs:
        if f == f00:
            continue

        mx = model.SpecModel.read_kurucz(f)
        mx = model.crop(mx,'Teff',trange)
        mx = model.crop(mx,'logg',lrange)
        fn = basename(f)
        mhx = float( str(fn[2]) + '.' + str(fn[3]) )
        if fn[1] == 'm':
            mhx *= -1
        mx = model.append_parameter(mx,'MH',mhx)
        m = model.concat(m,mx)

    m.write_model('kurucz_m',overwrite=True)


def phoenix_spectra():
    """Convolve PHOENIX spectra to lower resolution and combine [M/H].
        
    TODO: allow for separate convolution and combination.
    """

    s00 = model.SpecModel.read_model('phoenix-0.0')
    s00 = model.append_parameter(s00,'MH',0.0)

    for m in ['+0.5',-0.5,-1.0,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0]:

        s = model.SpecModel.read_model('phoenix'+str(m))
        s = model.append_parameter(s,'MH',float(m))
        s00 = model.concat(s00,s)

    s00.write_model('phoenix_m',overwrite=True)


def phoenix_mh_spectra(resolution=500,mh=0.0,overwrite=False):
    """Generate SpecModel from phoenix spectra.

    Teff and logg range is hardcoded to 2600-29,000K and 2-4.5. This is
    a range where the grid is rectangular over all metallicities and
    covers a wide range.

    BT-Settl model files are large, so read and downsample files one at
    a time to avoid having them all in memory at once.
    
    TODO: this takes hours, parallelize...
    """

    # models from 2600-29000K, logg 2-4.5, at [M/H]=0.0,
    # sorted by temperature and then gravity
    if mh == 0.0:
        mhstr = '{:+}'.format(-np.abs(mh))
    else:
        mhstr = '{:+}'.format(mh)

    name = 'phoenix'+mhstr
    
    fs = glob.glob(cfg.file['phoenix_models']
                   +'lte[0-2][0-9][0-9]-[2-4].?'+mhstr
                   +'a?[0-9].[0-9].BT-Settl.7.bz2')
    fs.sort()

    # read in and resample, one at a time
    filters = list(filter.Filter.all)
    conv_fnujy_sr = np.zeros((len(filters),len(fs)))
    spec = []
    teff = []
    logg = []
    for i,f in enumerate(fs):
        
        s = spectrum.ModelSpectrum.read_phoenix(f)
        teff.append(s.param_values['Teff'])
        logg.append(s.param_values['logg'])
        print("Read teff:{}, logg:{}, [M/H]:{} ({} of {})".
              format(s.param_values['Teff'],
                     s.param_values['logg'],
                     mhstr,
                     i+1,len(fs)) )
        
        # convolve spectrum to much lower resolution
        kern = s.resample(resolution=resolution)
        spec.append(s)

    # sort spectra
    teffarr = np.unique(teff)
    loggarr = np.unique(logg)

    s = model.SpecModel()
    s.name = spec[0].name
    s.wavelength = spec[0].wavelength
    s.parameters = ['Teff','logg']
    s.param_values = {'Teff':teffarr,
                      'logg':loggarr}
        
    s.fnujy_sr = np.zeros((len(s.wavelength),
                          len(teffarr),
                          len(loggarr)),dtype=float)

    for i,sp in enumerate(spec):
        if not np.all( np.equal(s.wavelength,sp.wavelength) ):
            raise SdfError("wavelength grids not the same \
                            in files {} and {}".format(fs[0],fs[i]))
        j = np.where(teff[i] == teffarr)[0][0]
        k = np.where(logg[i] == loggarr)[0][0]
        s.fnujy_sr[:,j,k] = sp.fnujy_sr

    s.write_model(name,overwrite=overwrite)

    return s
