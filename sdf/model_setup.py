"""Functions to set up models.

Any model can in principle be used and fit to photometry and spectra
using sdf, the models just need to be set up correctly. For the default
models most of that is done here.

Basic structure is that there are SpecModels (spectra) and PhotModels
(photometry in filters). PhotModels are derived from SpecModels.

To start from scratch, it may be beneficial to create a fresh folder
in which models will go, since generation can take a long time. Do this
by changing the file:model_root setting in .sdf.conf when model_setup
is imported (but change it back once that is done).

Models with reddening can be created from those without, e.g. with
phoenix_reddened_spectra.

In general, run setup_phot to regenerate PhotModels when filters have
been added or modified.

For generic details of how models are structured and how they work, see
:py:mod:`sdf.model`.

Models
------

Aside from Phoenix, the functions below create only SpecModels, once
these are run then specmodel2phot must be run to generate PhotModels.

### Phoenix
To generate phoenix model grids, first create resampled spectra with
resample_phoenix_spectra, which will loop over metallicities to create
models with one resolution. Two are needed in the standard setup, r=100
for SpecModels, and r=1000 from which the PhotModels are generated.

Then setup_default_phoenix will put all of these together to make the
full grids.

### Kurucz
Run kurucz_spectra, which simply converts a model grid file into the
format used by sdf.

### Blackbodies
Run bb_spectra and/or modbb_spectra, which just call functions in
sdf.model.SpecModel

### Simplified size distribution
Run sdf.model.SpecModel.sd_spectra

### "Real grain"
Run real_grain_spectra, which reads in a file that was previously
created from IDL codes.

"""

import glob
from os.path import exists, basename
from itertools import product
from multiprocessing import Pool

from scipy.io import readsav
import numpy as np
import astropy.units as u
import extinction as ext

from . import convolve
from . import model
from . import spectrum
from . import filter
from . import utils
from . import config as cfg

c_micron = u.micron.to(u.Hz, equivalencies=u.spectral())


def setup_default_phoenix():
    """Setup default phoenix models - T_{eff}, logg, and [M/H] grid.
        
    PhotModel is derived from high resolution spectra, and SpecModel is
    lower resolution because high resolution isn't necessary, and saves
    memory while running.
    
    Does both a Solar (0.0) and full metallicity range model.
    """

    # create the high resolution SpecModel
    phoenix_spectra(in_name_postfix='-r1000', name='phoenix_sol', overwrite=True,
                    m_range=[])
    # compute convolved photometry in all filters and write PhotModel
    specmodel2phot('phoenix_sol', overwrite_filters=True, overwrite_model=True)
    # create the low resolution SpecModel
    phoenix_spectra(in_name_postfix='-r1000', name='phoenix_sol', overwrite=True,
                    m_range=[])

    # create the high resolution SpecModel
    phoenix_spectra(in_name_postfix='-r1000', name='phoenix_m', overwrite=True)
    # compute convolved photometry in all filters and write PhotModel
    specmodel2phot('phoenix_m', overwrite_filters=True, overwrite_model=True)
    # create the low resolution SpecModel
    phoenix_spectra(in_name_postfix='-r1000', name='phoenix_m', overwrite=True)


def setup_phot(overwrite_filters=False, overwrite_model=True,
               model_ignore=['phoenix-', 'phoenix+']):
    """Rederive convolved models and write combined PhotModel to disk.
    
    Parameters
    ----------
    overwrite_filters : bool, optional
        Set to True to overwrite files for each bandpass.
    overwrite_model : bool, optional
        Set to True to overwrite the PhotModel that was generated.
    model_ignore : list of str, optional
        Skip models that include these strings.

    See Also
    --------
    sdf.model_setup.specmodel2phot : Function called for each model.
    sdf.model_setup.convolve_specmodel : Function that does the heavy lifting.
    sdf.filter_info : Where filters are set up.

    Notes
    -----
    This function needs to be run to incorporate new filters and/or 
    propagate zero point offsets into colours/indices.
    """
    for name in cfg.models['names']:
        do_it = True
        for ig in model_ignore:
            if ig in name:
                print("skipping {}".format(name))
                do_it = False
        if do_it:
            print("doing {}".format(name))
            specmodel2phot(name, overwrite_filters=overwrite_filters,
                           overwrite_model=overwrite_model)


def specmodel2phot(mname, overwrite_filters=False, overwrite_model=False):
    """Generate a PhotModel grid from SpecModel models."""
    
    convolve_specmodel(mname, overwrite=overwrite_filters)
    m = model.PhotModel.read_convolved_models(mname)
    m.write_model(m.name, overwrite=overwrite_model)


def convolve_specmodel(mname, overwrite=False):
    """Convolve a SpecModel to generate a set of ConvolvedModels.

    Parameters
    ----------
    mname : string
        Name of the model to convolve. A SpecModel that is equal to or
        wider than the wavelength range required by the filters must
        exist in config['model_root'].
        
    overwrite : bool, optional
        Overwrite files for each filter. If they exist already they will
        simply be skipped. Thus, if set to False only convolved models
        that do not exist will be written, which would be the desired
        behaviour when new filters have been added (in sdf.filter_info).

    See Also
    --------
    sdf.convolve : ConvolvedModel class
    sdf.filter_info : Where new filters are added.
    sdf.config : Where paths to models are specified.
    
    Notes
    -----
    The range of parameters in the ConvolvedModels is the same as in the
    SpecModel used.
    """

    m = model.SpecModel.read_model(mname)

    # flat list of all indices (cartesian product)
    i_list = ()
    for i in range(len(m.parameters)):
        i_list += (list(range(m.param_shape()[i])), )
    i_list = [i for i in product(*i_list)]

    # grid with spectral dimension last
    fnujy_sr_roll = np.rollaxis(m.fnujy_sr, axis=0, start=len(m.fnujy_sr.shape))

    # loop through them and create the ConvolvedModels and the files
    filters = filter.Filter.all
    fnujy_sr = np.zeros(len(i_list))
    for fname in filters:
        
        outfile = cfg.model_loc[mname]+fname+'.fits'
        
        if exists(outfile) and overwrite is False:
            print("Skipping {}, file exists".format(fname))
        
        else:
            print("Convolving filter {}".format(fname))
            for i, ind in enumerate(i_list):
                s = spectrum.ModelSpectrum(nu_hz=c_micron/m.wavelength,
                                           fnujy_sr=fnujy_sr_roll[ind])
                conv, cc = s.synthphot(fname)
                fnujy_sr[i] = conv
        
            if len(m.parameters) > 1:
                fnujy_sr_reshaped = np.reshape(fnujy_sr, m.param_shape())
            else:
                fnujy_sr_reshaped = fnujy_sr

            cm = convolve.ConvolvedModel(name=mname,
                                         filter_=fname,
                                         parameters=m.parameters,
                                         param_values=m.param_values,
                                         fnujy_sr=fnujy_sr_reshaped)

            cm.write_file(outfile, overwrite=overwrite)


def kurucz_spectra():
    """Generate a SpecModel grid of Castelli & Kurucz models."""
    
    # ranges, keep all [M/H]
    trange = [3499, 26001]
    lrange = [2.9, 5.1]
    
    # Solar metallicity
    f00 = cfg.file['kurucz_models']+'fp00k2odfnew.pck'
    m = model.SpecModel.read_kurucz(f00)
    m = model.crop(m, 'Teff', trange)
    m = model.crop(m, 'logg', lrange)
    m.write_model('kurucz-0.0', overwrite=True)

    # range of metallicities
    m = model.append_parameter(m, 'MH', 0.0)
    fs = glob.glob(cfg.file['kurucz_models']+'f*k2odfnew.pck')
    for f in fs:
        if f == f00:
            continue

        mx = model.SpecModel.read_kurucz(f)
        mx = model.crop(mx, 'Teff', trange)
        mx = model.crop(mx, 'logg', lrange)
        fn = basename(f)
        mhx = float(str(fn[2]) + '.' + str(fn[3]))
        if fn[1] == 'm':
            mhx *= -1
        mx = model.append_parameter(mx, 'MH', mhx)
        m = model.concat(m, mx)

    m.write_model('kurucz_m', overwrite=True)


def phoenix_spectra(in_name_postfix='', name='phoenix_m', overwrite=True,
                    m_range=['+0.5', -0.5, -1.0, -1.5, -2.0, -2.5, -3.0, -3.5, -4.0]):
    """Combine PHOENIX spectra with a range of [M/H] and write to disk.

    Parameters
    ----------
    name: string, optional
        The name of the combined model.

    in_name_postfix: string, optional
        String that may appear on the end of the phoenix models we want
        to combine, in addition to "phoenix-X.X". These models must have
        already been created by phoenix_mh_spectra using name_postfix.
        
    overwrite: bool, optional
        Write the combined model to disk. This option would only need to
        be set to False for testing whether some models will actually
        combine OK.
        
    m_range: list, optional
        List of metallicities in addition to 0.0 to include in model
    """

    s00 = model.SpecModel.read_model('phoenix-0.0'+in_name_postfix)
    if len(m_range) > 0:
        s00 = model.append_parameter(s00, 'MH', 0.0)

    for m in m_range:

        s = model.SpecModel.read_model('phoenix'+str(m)+in_name_postfix)
        s = model.append_parameter(s, 'MH', float(m))
        s00 = model.concat(s00, s)

    s00.write_model(name, overwrite=overwrite)


def resample_phoenix_spectra(resolution=100, name_postfix=''):
    """Resample all phoenix spectra to common wavelength grid.
        
    This will take about 2.5h per metallicity for R=100 with 8 cores on
    a 5k iMac, or about 20h for R=1000.
    
    For reference, Spitzer's IRS instrument has low resolution and high
    resolution modules at R=60-130 and R=600. JSWT MIRI has low and
    medium resolution at R~100 and R~1550-3250. For IRS the low 
    resolution mode was by far the most common. Thus, for spectra that
    will be resampled for these instruments R~100 is most sensible.
    """

    for m in [0.5, 0.0, -0.5, -1.0, -1.5, -2.0, -2.5, -3.0, -3.5, -4.0]:
        phoenix_mh_spectra(resolution=resolution, mh=m,
                           name_postfix=name_postfix, overwrite=True)


def phoenix_mh_spectra_one(par):
    """Read in and convolve one phoenix spectrum.
    
    This is a helper function so that phoenix_mh/cool_spectra can
    process a set of phoenix spectra more quickly.
    
    Parameters
    ----------
    par : tuple
        Tuple of wavelengths to resample to, and phoenix file name.
    
    See Also
    --------
    phoenix_mh_spectra
    phoenix_cool_spectra
    """
    
    wave, f = par
    print("{}".format(f))
    s = spectrum.ModelSpectrum.read_phoenix(f)
    kern = s.resample(wave)
    return s


def phoenix_mh_spectra(resolution=100, mh=0.0, overwrite=False,
                       processes=cfg.calc['cpu'], name_postfix=''):
    """Generate a SpecModel at some metallicity from phoenix spectra.

    For the BT-Settl-AGS2009 models Teff and logg range is hardcoded to
    2600-29, 000K and 2-4.5. This is a range where the grid is
    rectangular over metallicities and covers a wide range (+0.5 to -4 
    at 0.5 steps) .
    
    For these models there is also a much smaller set of "cool" models,
    which go from 2600K down to 400K, with a restricted and non-square
    range of logg and at Solar metallicity (most metallicities are
    present above 2000K, but -1.0 is missing). Most realistic set of
    models is probably logg=3.5 (with a few non-existent 3.5 ones 
    replaced 4.0).
    
    Parameters
    ----------
    resolution : float, optional
        Resolution of generated models.
    mh : float
        Metallicity of desired models.
    overwrite : bool, optional
        Force overwrite of extant models.
    processes : int, optional
        Number of simultaneous processes for calculation.
    name_postfix : str, optional
        String to append to model names.
        
    See Also
    --------
    phoenix_cool_spectra
    """

    # models from 2600-29000K, logg 2-4.5, at [M/H]=0.0,
    # sorted by temperature and then gravity
    if mh == 0.0:
        mhstr = '{:+}'.format(-np.abs(mh))
    else:
        mhstr = '{:+}'.format(mh)

    name = 'phoenix'+mhstr+name_postfix
    
    # don't do the calculation if there will be a write error
    if overwrite is False:
        if name in cfg.model_loc.keys():
            if exists(cfg.model_loc[name]+name+'_SpecModel.fits'):
                raise utils.SdfError("{} exists, will not overwrite".
                                     format(cfg.model_loc[name]+name+'.fits'))

    # the files for the main set of models
    fs = glob.glob(cfg.file['phoenix_models']
                   + 'lte[0-2][0-9][0-9]-[2-4].?' + mhstr
                   + 'a?[0-9].[0-9].BT-Settl.7.bz2')
    fs.sort()

    # get the new wavelength grid, and ensure it goes to the max
    wave = np.power(10, np.arange(np.log10(cfg.models['min_wav_micron']),
                                  np.log10(cfg.models['max_wav_micron']),
                                  np.log10(1+1/float(resolution))))
    if np.max(wave) != cfg.models['max_wav_micron']:
        wave = np.append(wave, cfg.models['max_wav_micron'])

    # read in and resample, in parallel
    pool = Pool(processes=processes)
    par = zip([wave for i in range(len(fs))], fs)
    spec = pool.map(phoenix_mh_spectra_one, par)
    pool.close()

    # sort spectra
    teff = [s.param_values['Teff'] for s in spec]
    logg = [s.param_values['logg'] for s in spec]
    teffarr = np.unique(teff)
    loggarr = np.unique(logg)

    s = model.SpecModel()
    s.name = spec[0].name
    s.wavelength = spec[0].wavelength
    s.parameters = ['Teff', 'logg']
    s.param_values = {'Teff': teffarr,
                      'logg': loggarr}
        
    s.fnujy_sr = np.zeros((len(s.wavelength),
                          len(teffarr),
                          len(loggarr)), dtype=float)

    for i, sp in enumerate(spec):
        if not np.all(np.equal(s.wavelength, sp.wavelength)):
            raise utils.SdfError("wavelength grids not the same \
                                 in files {} and {}".format(fs[0], fs[i]))
        j = np.where(teff[i] == teffarr)[0][0]
        k = np.where(logg[i] == loggarr)[0][0]
        s.fnujy_sr[:, j, k] = sp.fnujy_sr

    s.write_model(name, overwrite=overwrite)

    return s


def reddened_spectra(name='phoenix_m', overwrite=True):
    """Generate a SpecModel of spectra, with reddening."""

    m = model.SpecModel.read_model(name)
    m0 = model.append_parameter(m, 'A_V', 0.0)

    for av in [0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0]:

        mx = model.append_parameter(m, 'A_V', av)

        red = ext.apply(ext.odonnell94(m.wavelength*1e4, av, 3.1), 1)

        dim_array = np.ones((1, mx.fnujy_sr.ndim), int).ravel()
        dim_array[0] = -1
        red_reshaped = red.reshape(dim_array)
        
        mx.fnujy_sr *= red_reshaped
        
        m0 = model.concat(m0, mx)

    m0.write_model(name+'_av', overwrite=overwrite)


def phoenix_cool_spectra(resolution=100, overwrite=False,
                         processes=cfg.calc['cpu'], name_postfix=''):
    """Generate a SpecModel of cool phoenix spectra.

    Use a subset of the BT-Settl-AGS2009 models for cool stars.

    .. todo:: yes, there alarming jumps in the spectra, these are present
    in the raw files and I'd guess where the opacities aren't available.

    Parameters
    ----------
    resolution : float, optional
        Resolution of generated models.
    overwrite : bool, optional
        Force overwrite of extant models.
    processes : int, optional
        Number of simultaneous processes for calculation.
    name_postfix : str, optional
        String to append to model names.
        
    See Also
    --------
    phoenix_mh_spectra
    """

    name = 'phoenix_cool'+name_postfix
    
    # don't do the calculation if there will be a write error
    if overwrite is False:
        if name in cfg.model_loc.keys():
            if exists(cfg.model_loc[name]+name+'_SpecModel.fits'):
                raise utils.SdfError("{} exists, will not overwrite".
                                     format(cfg.model_loc[name]+name+'.fits'))

    # the files for the main set of models
    fs = glob.glob(cfg.file['phoenix_cool_models']+'*.BT-Settl.7.bz2')
    fs.sort()

    # get the new wavelength grid, and ensure it goes to the max
    wave = np.power(10, np.arange(np.log10(cfg.models['min_wav_micron']),
                                  np.log10(cfg.models['max_wav_micron']),
                                  np.log10(1+1/float(resolution))))
    if np.max(wave) != cfg.models['max_wav_micron']:
        wave = np.append(wave, cfg.models['max_wav_micron'])

    # read in and resample, in parallel
    pool = Pool(processes=processes)
    par = zip([wave for i in range(len(fs))], fs)
    spec = pool.map(phoenix_mh_spectra_one, par)
    pool.close()

    # sort spectra
    teff = [s.param_values['Teff'] for s in spec]
    teffarr = np.unique(teff)

    s = model.SpecModel()
    s.name = spec[0].name
    s.wavelength = spec[0].wavelength
    s.parameters = ['Teff']
    s.param_values = {'Teff': teffarr}
        
    s.fnujy_sr = np.zeros((len(s.wavelength),
                          len(teffarr)), dtype=float)

    for i, sp in enumerate(spec):
        if not np.all(np.equal(s.wavelength, sp.wavelength)):
            raise utils.SdfError("wavelength grids not the same \
                                 in files {} and {}".format(fs[0], fs[i]))
        j = np.where(teff[i] == teffarr)[0][0]
        s.fnujy_sr[:, j] = sp.fnujy_sr

    s.write_model(name, overwrite=overwrite)

    return s


def koester_wd_spectra(resolution=100, overwrite=False,):
    """Generate a SpecModel for Koester white dwarf spectra.

    Parameters
    ----------
    resolution : float, optional
        Resolution of generated models.
    overwrite : bool, optional
        Force overwrite of extant models.
    """

    name = 'koester_wd'

    # don't do the calculation if there will be a write error
    if overwrite is False:
        if name in cfg.model_loc.keys():
            if exists(cfg.model_loc[name]+name+'_SpecModel.fits'):
                raise utils.SdfError("{} exists, will not overwrite".
                                     format(cfg.model_loc[name]+name+'.fits'))

    # the files for the main set of models
    fs = glob.glob(cfg.file['koester_models']+'da*.dk.dat.txt')
    fs.sort()

    # get the new wavelength grid, and ensure it goes to the max
    wave = np.power(10, np.arange(np.log10(cfg.models['min_wav_micron']),
                                  np.log10(cfg.models['max_wav_micron']),
                                  np.log10(1+1/float(resolution))))
    if np.max(wave) != cfg.models['max_wav_micron']:
        wave = np.append(wave, cfg.models['max_wav_micron'])

    # read in and resample
    spec = []
    for f in fs:
        s = spectrum.ModelSpectrum.read_koester(f)
        s.resample(wave)
        spec.append(s)

    # sort spectra
    teff = [s.param_values['Teff'] for s in spec]
    logg = [s.param_values['logg'] for s in spec]
    teffarr = np.unique(teff)
    loggarr = np.unique(logg)

    s = model.SpecModel()
    s.name = spec[0].name
    s.wavelength = spec[0].wavelength
    s.parameters = ['Teff', 'logg']
    s.param_values = {'Teff': teffarr,
                      'logg': loggarr}

    s.fnujy_sr = np.zeros((len(s.wavelength),
                          len(teffarr),
                          len(loggarr)), dtype=float)

    for i, sp in enumerate(spec):
        if not np.all(np.equal(s.wavelength, sp.wavelength)):
            raise utils.SdfError("wavelength grids not the same \
                                 in files {} and {}".format(fs[0], fs[i]))
        j = np.where(teff[i] == teffarr)[0][0]
        k = np.where(logg[i] == loggarr)[0][0]
        s.fnujy_sr[:, j, k] = sp.fnujy_sr

    s.write_model(name, overwrite=overwrite)

    return s


def bb_spectra():
    """Generate SpecModel grid of blackbody models."""
    
    model.SpecModel.bb_disk_r(name='bb_disk_r',
                              write=True, overwrite=True)
    model.SpecModel.bb_disk_r(name='bb_star',
                              temperatures=10**np.arange(2.7, 3.5, 0.1),
                              write=True, overwrite=True)


def modbb_spectra():
    """Generate SpecModel grid of modified blackbody models."""
    
    model.SpecModel.modbb_disk_r(name='modbb_disk_r',
                                 write=True, overwrite=True)
    model.SpecModel.bb_disk_r(name='modbb_l210b1_disk_r',
                              write=True, overwrite=True,
                              lam0=210.0, beta=1.0)
    model.SpecModel.modbb_disk_dr(name='modbb_disk_dr',
                                  write=True, overwrite=True)


def real_grain_spectra(file, overwrite=False):
    """Real dust grain models from IDL save file.
    
    Files are made by sdf/dust_spectra.pro, saving a grid of P(r)
    with dimensions [wav, temp, dmin, q]. When restored the dimensions are
    reversed.
    """

    pr = readsav(file)

    s = model.SpecModel()
    s.name = 'amsil_r'
    s.wavelength = pr['wavs']
    s.parameters = ['log_Temp', 'log_Dmin', 'q']
    s.param_values = {'log_Temp': np.log10(pr['tbbs']),
                      'log_Dmin': np.log10(pr['dmins']),
                      'q': pr['qs']}

    s.fnujy_sr = np.rollaxis(np.rollaxis(np.rollaxis(pr['pr'], 1), 2), 3)
    s.fnujy_sr[s.fnujy_sr < cfg.tiny] = cfg.tiny

    s.write_model(s.name, overwrite=overwrite)
