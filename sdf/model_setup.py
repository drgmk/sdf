import glob
from os.path import exists
import numpy as np
from . import convolve
from . import model
from . import spectrum
from . import filter
from .utils import SdfError
from . import config as cfg

"""Helper functions to set up models"""

def setup_all():
    """Rederive all models."""
    setup_spec()
    setup_phot()
    phoenix_phot_spec_models(overwrite=True)


def setup_spec():
    """Rederive spectrum models."""
    bb_spectra()
    modbb_spectra()
    kurucz_spectra()


def setup_phot():
    """Rederive convolved models."""
    bb_phot()
    modbb_phot()
    kurucz_phot()


def bb_phot():
    """Generate a PhotModel grid of blackbody models.
        
    This is logT spaced and the model called 'bb'.
    """
    model.PhotModel.generate_bb_model(write=True,overwrite=True)


def modbb_phot():
    """Generate a PhotModel grid of modified blackbody models.
        
    This is logT spaced and the model called 'modbb'.
    """
    model.PhotModel.generate_modbb_model(write=True,overwrite=True)


def bb_spectra():
    """Generate SpecModel grid of blackbody models.
    
    This is logT spaced and the model called 'bb'.
    """
    model.SpecModel.generate_bb_model(write=True,overwrite=True)


def modbb_spectra():
    """Generate SpecModel grid of modified blackbody models.
    
    This is logT spaced and the model called 'modbb'.
    """
    model.SpecModel.generate_modbb_model(write=True,overwrite=True)


def kurucz_phot():
    """Generate a PhotModel grid of Castelli & Kurucz models.
        
    The models are called 'kurucz'.
    """
    convolve_kurucz(overwrite=True)
    m = model.PhotModel.read_convolved_models('kurucz')
    m.write_model(m.name,overwrite=True)


def kurucz_spectra():
    """Generate a SpecModel grid of Castelli & Kurucz models.
        
    The models are called 'kurucz'.
    """
    m=model.SpecModel.read_kurucz(cfg.file['kurucz_models']+'fp00k2odfnew.pck')
    m.write_model('kurucz',overwrite=True)


def convolve_kurucz(file='fp00k2odfnew.pck',overwrite=False):
    """Convolve a set of Castelli & Kurucz models.
        
    These are at Solar metallicity (by default). Files are written for
    each filter.
    """

    # get the models and ensure they are sorted in order
    m,te,lg,mh = spectrum.ModelSpectrum.read_kurucz(cfg.file['kurucz_models']+file)
    trange = [3500,26000]
    lrange = [3,5]
    ok = (te >= trange[0]) & (te <= trange[1]) & (lg >= lrange[0]) & (lg <= lrange[1])
    par = np.zeros( np.sum(ok), dtype=[('Teff',float),('logg',float)] )
    par['Teff'] = te[ok]
    par['logg'] = lg[ok]
    models = m[ok]
    srt = np.argsort(par,order=('Teff','logg'))
    par = par[srt]
    models = models[srt]
    
    # get the parameter values and reshape the grid of spectra
    teff = np.unique(par['Teff'])
    logg = np.unique(par['logg'])
    grid = models.reshape(len(teff),len(logg))
    
    # loop through them and create the ConvolvedModels and the files
    filters = filter.Filter.all
    fnujy_sr = np.zeros(grid.shape)
    for fname in filters:
        outfile = cfg.model_loc['kurucz']+fname+'.fits'
        if exists(outfile) and overwrite == False:
            print("Skipping {}, file exists".format(fname))
        else:
            print("Convolving filter {}".format(fname))
            cm = convolve.ConvolvedModel(name='kurucz',
                                         filter=fname,parameters=['Teff','logg'],
                                         param_values={'Teff':teff,'logg':logg},
                                         fnujy_sr=fnujy_sr)
            for i in range(len(teff)):
                for j in range(len(logg)):
                    conv,cc = grid[i,j].synthphot(fname)
                    fnujy_sr[i,j] = conv
        
            cm.fnujy_sr = fnujy_sr
            cm.write_file(outfile,overwrite=overwrite)


def phoenix_phot_spec_models(overwrite=False):
    """Generate models from phoenix spectra.

    BT-Settl model files are large, so do both PhotModel and SpecModel
    processing at once, reading and downsampling files one at a time to
    avoid having them all in memory at once.
    """

    name = 'phoenix'
    
    # don't do the calculation if there will be a write error
    if overwrite == False:
        if ( exists(cfg.model_loc[name]+name+'_PhotModel.fits') or
            exists(cfg.model_loc[name]+name+'_SpecModel.fits') ):
            raise SdfError("model file(s) exist, will not overwrite")

    # models from 2600-29000K, logg 2-4.5, at [M/H]=0.0,
    # sorted by temperature and then gravity
    fs = glob.glob(cfg.file['phoenix_models']
                   +'lte[0-2][0-9][0-9]-[2-4]*BT-Settl.7.bz2')
    fs.sort()

    # read in, convolve, resample, one at a time
    filters = list(filter.Filter.all)
    conv_fnujy_sr = np.zeros((len(filters),len(fs)))
    spec = []
    teff = []
    logg = []
    for i,f in enumerate(fs):
        
        s = spectrum.ModelSpectrum.read_phoenix(f)
        teff.append(s.param_values['Teff'])
        logg.append(s.param_values['logg'])
        print("Read teff:{}, logg:{} ({} of {})".format(s.param_values['Teff'],
                                                        s.param_values['logg'],
                                                        i+1,len(fs)))
        
        # get convolved fluxes
        conv,_ = s.synthphot(filters)
        conv_fnujy_sr[:,i] = conv
        
        # checking for same wavelength grids (i.e. can we reuse kernel)
        if i > 0:
            if len(s.wavelength) != len(lastwav):
                kern = None
            elif not np.all( np.equal(s.wavelength,lastwav)):
                kern = None
            else:
                print("wave grids in {} and {} the same".format(fs[i-1],fs[i]))
        else:
            kern = None
        lastwav = s.wavelength

        # and spectrum (at much lower resolution)
        kern = s.resample(resolution=100,kernel=kern)
        spec.append(s)

    teffarr = np.unique(teff)
    loggarr = np.unique(logg)

    # sort photometry
    conv = conv_fnujy_sr.reshape( (len(filters),
                                   len(np.unique(teff)),
                                   len(np.unique(logg))) )
    cm = model.PhotModel(name=name,filters=filters,
                         parameters=['Teff','logg'],
                         param_values={'Teff': teffarr,
                                       'logg': loggarr},
                         fnujy_sr=conv)
    cm.write_model(name,overwrite=overwrite)

    # sort spectra
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
            raise SdfError("wavelength grids not the same in files {} and {}".format(fs[0],fs[i]))
        j = np.where(teff[i] == teffarr)[0][0]
        k = np.where(logg[i] == loggarr)[0][0]
        s.fnujy_sr[:,j,k] = sp.fnujy_sr

    s.write_model(name,overwrite=overwrite)
