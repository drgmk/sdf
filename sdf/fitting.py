from functools import lru_cache
import os
import glob

import binarytree as bt
import numpy as np
import emcee

# in case we don't have this module
try:
    import classifier.photometry
    import classifier.spectra
    classifier_module = True
except ImportError:
    classifier_module = False

from . import model
from . import photometry
from . import spectrum
from . import filter
from . import utils
from . import result
from . import db
from . import config as cfg

# these are for getting info to pymultinest
global_obs = ()
global_mod = ()
global_p_rng = ()


def fit_results(file,update_mn=False,update_an=False,
                update_json=False,update_thumb=False,
                sort=True,custom_sort=True,nospec=False):
    """Return a list of fitting results.
        
    Parameters
    ----------
    file : str
        The raw photometry file to use as input.
    update_mn : bool, optional
        Force update of multinest fitting.
    update_an : bool, optional
        Force update of post-multinest fitting analysis.
    sort : bool, optional
        Sort results by decreasing evidence.
    custom_sort: bool, optional
        Additonally sort results using per-target config.
    nospec : bool, optional
        Exclude observed specta from fitting (for speed).
    """

    print(" Fitting")
    results = []

    # fit any extra models that we'll append later
    extra = []
    if len(cfg.fitting['extra_models']) > 0:
        for m in cfg.fitting['extra_models']:

            print("  ",m)
            r = result.Result.get(
                            file,m,update_mn=update_mn,
                            update_an=update_an,update_json=update_json,
                            update_thumb=update_thumb,nospec=nospec
                                           )

            # check for files with no photometry
            if not hasattr(r,'obs'):
                print("  no photometry = no results")
                return None

            extra.append(r)

    # get results for models defined by conf, overrides default
    if len(cfg.fitting['models']) > 0:
        for m in cfg.fitting['models']:

            print("  ",m)
            r = result.Result.get(
                            file,m,update_mn=update_mn,
                            update_an=update_an,update_json=update_json,
                            update_thumb=update_thumb,nospec=nospec
                                           )

            # check for files with no photometry
            if not hasattr(r,'obs'):
                print("  no photometry = no results")
                return None

            results.append(r)

    else:
        # binary tree-based fitting
        t = model_director(file)
        try:
            print_model_tree(t)
        except:
            print_model_tree(t, cute=False)

        while t.left is not None and t.right is not None:

            print("  ",t.left.value,"vs.",t.right.value)

            r1 = result.Result.get(
                        file,t.left.value,update_mn=update_mn,
                        update_an=update_an,update_json=update_json,
                        update_thumb=update_thumb,nospec=nospec
                                            )
            r2 = result.Result.get(
                        file,t.right.value,update_mn=update_mn,
                        update_an=update_an,update_json=update_json,
                        update_thumb=update_thumb,nospec=nospec
                                            )

            # check for files with no photometry
            if not hasattr(r1,'obs'):
                print("  no photometry = no results")
                return None

            # append results, only append left result at start since where
            # on lower branches the left model has already been done
            if len(results) == 0:
                results.append(r1)
            results.append(r2)

            # move on down the tree
            if r2.evidence > r1.evidence + cfg.fitting['ev_threshold']:
                t = t.right
            else:
                t = t.left

    # sort list of results by evidence (required for subsequent custom sort)
    print(' Sorting')
    if sort or custom_sort:
        print('   sorting results by evidence')
        results = [results[i] for i in result.sort_results(results)]
    else:
        print('   no results sorting')

    # sort list of results by custom method
    if custom_sort:
        print('   applying db.custom_sort results sorting')
        srt = db.custom_sort(file, results)
        if srt is None:
            print("     couldn't get config")
        else:
            results = [results[i] for i in srt]

    # append extra results
    if len(extra) > 0:
        print(' Appending')
        [print('   {}'.format(m)) for m in cfg.fitting['extra_models']]
        results = results + extra

    # save a thumb of the best fit next to the input file, update every
    # time since we may have changed best fit (but not fitting itself)
    results[0].sed_thumbnail(file='{}/{}_thumb.png'.format(results[0].path, results[0].id), update=True)

    return results


def model_director(file,reddening=False,use_classifier=False):
    """Workflow for model fitting.

    Parameters
    ----------
    file : str
        Name of the photometry file we are fitting.
    reddening : bool, optional
        Use models with reddening.
    use_classifier : bool, optional
        Use classifier.
    """

    # default model tries star + up to two bb components
    if reddening:
        star = 'phoenix_m_av'
    else:
        star = 'phoenix_m'

    t_star = model_tree(top=(star,), extra='modbb_disk_r',n_extra=2)

    # cool star model
    if reddening:
        cool = 'phoenix_cool_av'
    else:
        cool = 'phoenix_cool'

    t_cool = model_tree(top=(cool,), extra='modbb_disk_r')

    # look for spectral type, LTY types get cool models, other types
    # default to star models, and M5-9 (or just M) get both
    tree = t_star
    cool = 0
    kw = utils.get_sdb_keywords(file)
    if 'sp_type' in kw.keys():
        if kw['sp_type'] is None:
            pass
        elif kw['sp_type'][:2] == 'DA' or kw['sp_type'][:2] == 'DB':
            tree = model_tree(top=('koester_wd',), extra='bb_disk_r')
        elif kw['sp_type'][0] in 'LTY':
            tree = t_cool
        elif kw['sp_type'][0] == 'M':
            if len(kw['sp_type']) > 1:
                if kw['sp_type'][1] in '56789':
                    cool = 1
        elif kw['sp_type'][0:2] == 'dM':
            if len(kw['sp_type']) > 2:
                if kw['sp_type'][2] in '56789':
                    cool = 1

    if cool == 1:
        tree = bt.Node(('top',))
        tree.left = t_cool
        tree.right = t_star

    tree.value = ('top',)
    return tree


def model_tree(top=('phoenix_m',),extra='modbb_disk_r',n_extra=1):
    """Return a binary tree for alternative models.
        
    Parameters
    ----------
    top : str, optional
        Name of the first model.
    extra : str, optional
        Name of the additional model component.
    n_extra : int, optional
        Include two extra component branch.
    """

    # the top node doesn't matter unless this tree becomes a branch
    # of a bigger tree
    t = bt.Node( top )
    t.left = bt.Node( top )
    t.right = bt.Node( top + (extra,) )
    if n_extra == 2:
        t.right.left = bt.Node( top + (extra,) )
        t.right.right = bt.Node( top + (extra, extra) )

    return t


def print_model_tree(t, cute=True):
    """Print a model tree, shortening the model names.
    
    Unicode symbols here https://unicode-table.com/en/
    
    Parameters
    ----------
    t : binary tree
        Tree to print.
    cute : bool, optional
        Print cute ascii symbols.
    """

    if cute:
        r = {'top':u'\u2602',
             'phoenix_m':u'\u2606',
             'phoenix_m_av':u'\u2605',
             'phoenix_cool':u'\u2733',
             'phoenix_cool_av':u'\u2739',
             'modbb_disk_r':u'\u29b8',
             'bb_disk_r':u'\u25cb'
             }
    else:
        r = {'top':'t',
             'phoenix_m':'p',
             'phoenix_m_av':'pr',
             'phoenix_cool':'c',
             'phoenix_cool_av':'cr',
             'modbb_disk_r':'mb',
             'bb_disk_r':'b'
             }

    l = bt.convert(t)

    for i in range(len(l)):
        x = ()
        if l[i] is not None:
            for j in range(len(l[i])):
                if l[i][j] in r.keys():
                    x += (r[l[i][j]],)
                else:
                    x += (l[i][j],)
            l[i] = x = ''.join(x)

    bt.convert(l).show()


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
    # level, ignored and otherwise zero as set above
    # see Johnson+2013, MNRAS 436, 2535 for some discussion
    for i,lim in enumerate(obs_uplim):
        if lim and not obs_ignore[i]:
            if mod_fnu[i] > obs_fnu[i]/3.:
                resid[i] = -1. * mod_fnu[i] / (obs_fnu[i]/3.)

    return resid,obs_wav,obs_filt


def residual_phot(param,*args):
    """Return residuals for a model compared to observation.
    
    This is a version for only photometry.
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
    
    import pymultinest as pmn
    
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


def pmn_pc(prob,samples,pcs,axis=0):
    """Return numpy-like percentile for multinest output.
        
    Accounts for probability of samples, which must therefore be given.
    Mostly a copy from analyzer.get_stats() in pymultinest.
    """

    # sort what to do depending on sample array dimension
    if np.ndim(samples) == 3:
    
        # loop over each column
        if axis == 2:
            out = np.zeros((len(pcs),samples.shape[0],samples.shape[1]))
            for i in range(samples.shape[0]):
                for j in range(samples.shape[1]):
                    out[:,i,j] = pmn_pc(prob,samples[i,j,:],pcs)
            
            return out
        else:
            raise utils.SdfError("axis must be 2 for 3d samples")

    if np.ndim(samples) == 2:
    
        if axis == 1:
            return pmn_pc(prob,samples.T,pcs) # do for the transpose
        
        # loop over each column
        elif axis == 0:
            out = np.zeros((len(pcs),samples.shape[1]))
            for i in range(samples.shape[1]):
                out[:,i] = pmn_pc(prob,samples[:,i],pcs)
            
            return out
        else:
            raise utils.SdfError("axis must be 0 or 1")

    elif np.ndim(samples) > 2:
        raise utils.SdfError("can't do more than 2d")
    elif np.ndim(samples) == 0:
        raise utils.SdfError("need an array/list/tuple of values")
    else:
        # organise and sort probabilities and samples
        if len(prob) != len(samples):
            raise utils.SdfError("prob and samples have different lengths"
                                 " {} and {}".format(len(prob),len(samples)))
        b = list(zip(prob,samples))
        b.sort(key=lambda x: x[1])
        b = np.array(b)
        b[:,0] = b[:,0].cumsum()
        
        # additional normalisation step, since expect to use sub-samples
        b[:,0] /= np.max(b[:,0])
        
        # interpolation function to get percentiles
        bi = lambda x: np.interp(x, b[:,0], b[:,1], left=b[0,1], right=b[-1,1])

        if isinstance(pcs,(int,float)):
            return bi(pcs/100.0)
        elif isinstance(pcs,(np.ndarray,list,tuple)):
            return np.array([bi(pc/100.0) for pc in pcs])
        else:
            raise utils.SdfError("wrong type {}".format(type(pcs)))


def weighted_dist(prob):
    """Select samples from distribution, accounting for weights.

    Returns boolean array with same length as input, with True
    indicating samples to keeep.

    Parameters
    ----------
    prob : np.array
        Array of probabilities
    """

    p_rand = np.random.uniform(high=np.max(prob),size=len(prob))
    return prob > p_rand


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


def emcee_prior(param,m):
    """Prior to keep parameters in range allowed by model."""

    m_info = model.models_info(m)
    p_rng = m_info['p_rng']
    for p,rng in zip(param,p_rng):
        if p < rng[0] or p > rng[1]:
            return -np.inf

    return 0.0


def run_emcee(r,nwalkers=8,nstep=100,start_pos=None):
    """Run emcee MCMC fitting."""

    if r.models == '':
        raise utils.SdfError('result.models empty, use sdf.model.get_models()')

    if start_pos is None:
        start_pos = [r.best_params + \
                     r.best_params_1sig*np.random.normal(size=r.model_info['ndim'])
                     for i in range(nwalkers)]

    def emcee_lnlike(param,*args):
        o,m = args
        return emcee_prior(param,m) + lnlike(param,*args)

    sampler = emcee.EnsembleSampler(nwalkers,r.model_info['ndim'],
                                    emcee_lnlike,args=(r.obs,r.models))
    pos,lnprob,rstate = sampler.run_mcmc(start_pos, nstep)

    return pos,sampler


def find_outliers(r, max_wav=3):
    '''Find probable outliers after fitting.
    
    Loops through observations shorter than some wavelength, ignoring
    each point, rescaling the model without this point, and computes
    the chi^2. Tests whether the lowest point is significantly
    different to the others.
    '''

    ok = (filter.mean_wavelength(r.filters) < max_wav) & \
          np.invert(r.filters_ignore)
    
    chi_ref = np.sum( ((r.obs_fnujy[ok] - r.model_fnujy[ok])/r.obs_e_fnujy[ok])**2 )

    dchi = np.zeros(len(r.obs_fnujy))
    for i in range(len(dchi)):
        if not ok[i]:
            continue
        
        ok_tmp = ok.copy()
        ok_tmp[i] = False
        model_tmp = r.model_fnujy * np.mean(r.obs_fnujy[ok_tmp]) / np.mean(r.model_fnujy[ok_tmp])
        dchi[i] = np.sum( ((r.obs_fnujy[ok_tmp] - model_tmp[ok_tmp])/r.obs_e_fnujy[ok_tmp])**2 )

    mini = np.argmin(dchi[ok])
    print('most likely outlier: {}: {}'.format(mini, r.filters[mini]))
    ok[mini] = False
    print('{:f} vs. range of {:f} to {:f}'.format(dchi[mini],
                                                  np.min(dchi[ok]),
                                                  np.max(dchi[ok])))

    return chi_ref, dchi
