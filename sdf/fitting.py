from functools import lru_cache
import os
import glob

import numpy as np
import pymultinest as pmn

from . import model
from . import photometry
from . import spectrum
from . import filter
from . import utils
from . import config as cfg

# these are for getting info to pymultinest
global_obs = ()
global_mod = ()
global_p_rng = ()


@lru_cache(maxsize=2)
def concat_obs(o):
    """Concatenate observations

    Concatenate the observations (filters and spectra), not
    shifting the normalisation of the spectra, but noting how
    many observations there are in each component (as the models
    will probably contain more due to colours/indices).

    """
    
    obs_wav = np.array([],dtype=float)
    obs_filt = np.array([],dtype=object)
    obs_fnu = np.array([],dtype=float)
    obs_e_fnu = np.array([],dtype=float)
    obs_uplim = np.array([],dtype=bool)
    obs_ignore = np.array([],dtype=bool)
    obs_bibcode = np.array([],dtype=str)
    obs_nel = np.array([],dtype=int)
    ispec = -2 # start one extra from end, we will put -1 there for phot
    obs_ispec = np.array([],dtype=int)
    for obs in o:
        if isinstance(obs,photometry.Photometry):
            obs_wav = np.append(obs_wav,obs.mean_wavelength())
            obs_filt = np.append(obs_fnu,obs.filters)
            obs_fnu = np.append(obs_fnu,obs.fnujy)
            obs_e_fnu = np.append(obs_e_fnu,obs.e_fnujy)
            obs_ispec = np.append(obs_ispec,np.repeat(-1,len(obs.filters)))
            obs_uplim = np.append(obs_uplim,obs.upperlim)
            obs_ignore = np.append(obs_ignore,obs.ignore)
            obs_bibcode = np.append(obs_bibcode,obs.bibcode)
        elif isinstance(obs,spectrum.ObsSpectrum):
            n = len(obs.wavelength)
            obs_wav = np.append(obs_wav,obs.wavelength)
            obs_filt = np.append(obs_filt,np.repeat(None,n))
            obs_fnu = np.append(obs_fnu,obs.fnujy)
            obs_e_fnu = np.append(obs_e_fnu,obs.e_fnujy)
            obs_ispec = np.append(obs_ispec,np.repeat(ispec,n))
            ispec -= 1
            obs_uplim = np.append(obs_uplim,np.zeros(n,dtype=bool))
            obs_ignore = np.append(obs_ignore,np.zeros(n,dtype=bool))
            obs_bibcode = np.append(obs_bibcode,np.repeat(obs.bibcode,n))
        obs_nel = np.append(obs_nel,len(obs.fnujy))

    return (obs_fnu,obs_e_fnu,obs_uplim,obs_ignore,obs_ispec,
           obs_nel,obs_wav,obs_filt,obs_bibcode)


def residual(param,*args):
    """Return residuals for a model compared to observation.
        
    The structure of the observation and models needs to
    match, and the filters and/or wavelengths present within
    each of the correspondong observation and model obects
    must match.
    
    For pure photometry the args would be:
    (Photometry,),(PhotModel,)
    
    Where there are spectra, the parameters for a given set of
    (PhotModel,SpecModel) are the same, and the args would be:
    (Photometry,ObsSpectrum),((PhotModel,SpecModel),)
    
    For combining multiple models the args would be:
    (Ph,ObsSp),((PhMod,SpMod),(PhMod,SpMod))
    
    The parameters are ordered similarly, with each tuple of
    models sharing the same parameters, and with the normalisations
    for observed spectra last, in the same order as the spectra.
    
    To deal with colours/indices, which need the model components
    to be added first, and then derived, PhotModels have extra
    elements with the colour/index base filters, which are used
    to derive the colours/indices as we go.
    
    """

    o,m = args # observations and models
    
    # concatenate observations
    obs_fnu,obs_e_fnu,obs_uplim,obs_ignore,obs_ispec,\
            obs_nel,obs_wav,obs_filt,_ = concat_obs(o)
    
    # multiply spectra by appropriate normalisation, ispec starts at
    # -2 so we can add 1.0 for photometry at the end of the params
    spec_norm = np.take(np.append(param,1.0),obs_ispec)
    obs_fnu = obs_fnu * spec_norm
    obs_e_fnu = obs_e_fnu * spec_norm

    # get model fluxes, including filling of colours/indices
    mod_fnu,_ = model.model_fluxes(m,param,obs_nel)

    # residuals in significance units, setting zero where (photometry)
    # is to be ignored, and for upper limits (but amended below)
    resid = np.zeros(len(obs_fnu))
    ok = np.invert( np.any([obs_uplim,obs_ignore],axis=0) )
    resid[ok] = (obs_fnu[ok] - mod_fnu[ok]) / obs_e_fnu[ok]

    # set residual if any 3sigma upper limits exceeded at the 1sigma
    # level, otherwise zero as set above
    for i,lim in enumerate(obs_uplim):
        if lim:
            if mod_fnu[i] > obs_fnu[i]/3.:
                resid[i] = -1. * mod_fnu[i] / (obs_fnu[i]/3.)

    return resid,obs_wav,obs_filt


def residual_phot(param,*args):
    """Return residuals for a model compared to observation.
    
    This is an unused version for only photometry.
    """

    p,m = args # photometry and model
    
    # get model fluxes
    if isinstance(m,(tuple,list)):
        model_fnujy = np.zeros(p.nused)
        i0 = 0
        for mod in m:
            nparam = len(mod.parameters)+1
            flux = mod.fnujy(param[i0:i0+nparam])
            model_fnujy = model_fnujy + flux
            i0 += nparam
    else:
        model_fnujy = m.fnujy(param)

    # residuals in significance units
    resid = np.zeros(p.nused)
    ok = np.invert(p.upperlim)
    resid[ok] = (p.fnujy[ok] - model_fnujy[ok]) / p.e_fnujy[ok]

    # set residual if any upper limits exceeded, otherwise zero
    # (e_fnujy may be zero so resid may be nan)
    if np.any(p.upperlim):
        lim = p.upperlim
        if model_fnujy[lim] > p.fnujy[lim]:
            resid[lim] = (p.fnujy[lim]-model_fnujy[lim])/(p.fnujy[lim]/3.)
        else:
            resid[lim] = np.zeros(np.sum(p.upperlim))

    print(resid)
    return resid


def chisq(param,*args):
    """Return sum of squared residuals."""
    
    res,_,_ = residual(param,*args)
    return np.sum( np.square( res ) )


def lnlike(param,*args):
    """Return log likelihood."""
    
    chi2 = chisq(param,*args)
    if np.isfinite(chi2):
        return -0.5 * chi2
    else:
        return -np.inf


def multinest_prior(cube,ndim,nparam):
    """Prior for pymultinest, turns values given in each element
    of cube from range 0-1 to the range for that parameter
    """
    
    # get parameter ranges
    global global_p_rng
    pars = global_p_rng

    for i in range(ndim):
        cube[i] = pars[i][0] + cube[i] * (pars[i][1]-pars[i][0])


def multinest_lnlike(cube,ndim,nparam):
    """Return log likelihood."""
    
    global global_obs,global_mod
    o,m = (global_obs,global_mod)
    param = np.array([])
    for i in range(ndim):
        param = np.append(param,cube[i])
    return lnlike(param,o,m)


def multinest(o,m,dir):
    """Run pymultinest to fit model(s) to photometry."""
    
    dir = dir.rstrip('/')
    
    m_info = model.models_info(m)
    pmn_out = dir+'/'+m_info['name']+cfg.fitting['pmn_model_suffix']
    global global_obs,global_mod,global_p_rng
    global_obs = o
    global_mod = m
    global_p_rng = m_info['p_rng']

    pmn.run(multinest_lnlike,multinest_prior,m_info['ndim'],
            n_live_points=cfg.fitting['n_live'],
            n_iter_before_update=cfg.fitting['n_update'],
            multimodal=True,sampling_efficiency=0.3,
            verbose=cfg.fitting['verb'],
            outputfiles_basename=pmn_out)


def pmn_models(dir):
    """Return the models fitted by multinest."""

    fs = glob.glob( dir + '/*' + cfg.fitting['pmn_model_suffix'] + '.txt' )
    models = ()
    mbase = []
    for f in fs:
        pmn_out = f.rstrip('.txt')
        mbase.append(pmn_out)
        tmp = pmn_out.rstrip(cfg.fitting['pmn_model_suffix'])
        tmp = os.path.basename(tmp)
        models += (tmp.split(cfg.fitting['model_join']),)

    return mbase,models


def sort_evidence(ev_in,ndim):
    """Return argsort for models using evidence.

    For a model with more parameters to be preferred, the log evidence
    must be more than ev_threshold higher (just higher if the number
    of dimensions is the same). Method is to initially sort models by
    evidence, and then go through and change things based on these
    criteria we are done.
    """

    if len(ev_in) != len(ndim):
        raise utils.SdfError("length of ev_in ({}) and ndim ({}) not equal".\
                       format(ev_in,ndim))
    if isinstance(ev_in,list):
        ev_in = np.array(ev_in)
    if isinstance(ndim,list):
        ndim = np.array(ndim)

#    print(ev_in,ndim)
    srt = np.argsort(ndim)
    ev_srt = ev_in[srt]
    dim_srt = ndim[srt]
    order = np.array(srt)
    last = np.zeros(len(order))
#    print(ev_srt,order)

    while not np.all( np.equal(order,last) ):
        last = np.array(order)
        for i in range(len(ndim)-1):
            if ev_in[order][i+1] > ev_in[order][i]+cfg.fitting['ev_threshold']\
               or ndim[order][i+1] == ndim[order][i]\
                  and ev_in[order][i+1] > ev_in[order][i]:
                tmp = order[i+1]
                order[i+1] = order[i]
                order[i] = tmp
#            print(i,ev_in[order])

#    print(order)
    return order
