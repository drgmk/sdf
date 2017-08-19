from functools import lru_cache
import os.path
import pickle
import glob
import time
import json

import numpy as np
from scipy.stats import truncnorm
from scipy.optimize import minimize
import pymultinest as pmn
import matplotlib.pyplot as plt
import corner
import astropy.units as u

from . import photometry
from . import spectrum
from . import model
from . import filter
from . import fitting
from . import plotting
from . import utils
from . import config as cfg


class BaseResult(object):
    """Basic class to compute and handle fitting results."""

    @lru_cache(maxsize=128)
    def __init__(self,rawphot,model_comps):
        """Basic instantiation of the Result object."""
        
        self.file_info(rawphot,model_comps)
    
    
    @lru_cache(maxsize=128)
    def file_info(self,rawphot,model_comps):
        """Basic file info."""
                  
        # component info
        self.model_comps = model_comps
        self.star_or_disk = ()
        for comp in model_comps:
            if comp in cfg.models['star']:
                self.star_or_disk += ('star',)
            elif comp in cfg.models['disk']:
                self.star_or_disk += ('disk',)
            else:
                raise utils.SdfError("couldn't assigm comp {} to star or disk "
                               "given lists in {} and {}".
                               format(comp,cfg.models['star'],
                                      cfg.models['disk']))
    
        self.n_comps = len(model_comps)
        
        # where the rawphot file is
        self.rawphot = rawphot
        self.path = os.path.dirname(rawphot)
        
        # id
        self.id = os.path.basename(rawphot).rstrip('-rawphot.txt')

        # where the multinest output is (or will be), create if needed
        self.pmn_dir = self.path + '/' + self.id            \
                       + cfg.fitting['pmn_dir_suffix']
        if not os.path.exists(self.pmn_dir):
            os.mkdir(self.pmn_dir)
        
        # the base name for multinest files
        self.pmn_base = self.pmn_dir + '/'                  \
                        + cfg.fitting['model_join'].join(self.model_comps) \
                        + cfg.fitting['pmn_model_suffix']

        # plot names, pickle, and json, files may not exist yet
        self.corner_plot = self.pmn_base+'corner.png'
        self.distributions_plot = self.pmn_base+'distributions.png'
        self.sed_thumb = self.pmn_base + 'sed_thumb.png'
        self.pickle = self.pmn_base + '.pkl'
        self.json = self.pmn_base + '.json'

        # rawphot modification time
        self.rawphot_time = os.path.getmtime(self.rawphot)


    def fill_data_models(self,nospec=False):
        """Get photometry/spectra and the models needed to fit these.
        
        Read in the observations and corresponding models. We must do
        this since models are not saved in the pickle to save space.

        Parameters
        ----------
        nospec : bool, optional
            Don't include any spectra when reading in the observations.
        """

        # observations; keywords, tuples of photometry and spectra. if
        # there is nothing in the photometry file then don't fill
        # anything else
        p = photometry.Photometry.read_sdb_file(self.rawphot)
        if p is None:
            return
        elif np.sum(p.ignore) == p.nphot:
            return

        self.obs = (p,)
        self.obs_keywords = utils.get_sdb_keywords(self.rawphot)
        self.exclude_spectra = nospec
        if not nospec:
            s = spectrum.ObsSpectrum.read_sdb_file(self.rawphot,
                                                   module_split=True,
                                                   nspec=1)
            if s is not None:
                self.obs = (p,) + s

        # models
        mod,plmod = model.get_models(self.obs,self.model_comps)
        self.models = mod
        self.pl_models = plmod
        self.model_info = model.models_info(self.models)


    def fill_observations(self):
        """Fill observed attributes.
        
        Returns obs_nel, needed to compute model fluxes.
        """
        
        # this is largely copied from fitting.residual
        tmp = fitting.concat_obs(self.obs)
        self.obs_fnujy,self.obs_e_fnujy,self.obs_upperlim,self.filters_ignore,\
            obs_ispec,obs_nel,self.wavelengths,self.filters,self.obs_bibcode = tmp
        spec_norm = np.take(self.best_params+[1.0],obs_ispec)
        self.obs_fnujy = self.obs_fnujy * spec_norm
        self.obs_e_fnujy = self.obs_e_fnujy * spec_norm

        return obs_nel


    def fill_best_fit_spectra(self):
        """Fill best fit model spectra.
        
        .. todo:: resample spectrum rather than interpolate model

        """

        # ObsSpectrum for each component, at original resolution
        wave = cfg.models['default_wave']
        star_spec = np.zeros(len(wave))
        disk_spec = np.zeros(len(wave))
        total_spec = np.zeros(len(wave))
        self.comp_spectra = ()
        for i,comp in enumerate(self.pl_models):
            for mtmp in comp:
                if not isinstance(mtmp,model.SpecModel):
                    continue
                
                m = mtmp.copy()

                s = spectrum.ObsSpectrum(wavelength=m.wavelength,
                                         fnujy=m.fnujy(self.comp_best_params[i]))
                s.fill_irradiance()
                self.comp_spectra += (s,)

                # TODO: resample spectrum rather than interpolate model
                m.interp_to_wavelengths(wave)
                total_spec += m.fnujy(self.comp_best_params[i])
                if self.star_or_disk[i] == 'star':
                    star_spec += m.fnujy(self.comp_best_params[i])
                elif self.star_or_disk[i] == 'disk':
                    disk_spec += m.fnujy(self.comp_best_params[i])

        # and total/star/disk spectra, at common wavelengths
        self.total_spec = spectrum.ObsSpectrum(wavelength=wave,fnujy=total_spec)
        if np.max(star_spec) > 0:
            self.star_spec = spectrum.ObsSpectrum(wavelength=wave,fnujy=star_spec)
        else:
            self.star_spec = None

        if np.max(disk_spec) > 0:
            self.disk_spec = spectrum.ObsSpectrum(wavelength=wave,fnujy=disk_spec)
        else:
            self.disk_spec = None


    def main_results_text(self):
        """Return nicely formatted tuple of text of results."""
    
        # the array sets the order, and the dict the conversion
        text_ord = ['Teff','lstar','rstar',
                    'Temp','rdisk_bb','ldisk_lstar',
                    'lam0','beta','Dmin','q']
                    
        text_sub = {'Teff': ['T<sub>star</sub>','K'],
                    'MH':   ['[M/H]',''],
                    'logg': ['logg',''],
                    'lstar':['L<sub>star</sub>','L<sub>Sun</sub>'],
                    'rstar':['R<sub>star</sub>','R<sub>Sun</sub>'],
                    'Temp': ['T<sub>dust</sub>','K'],
                    'lam0': ['&lambda;<sub>0</sub>','&mu;m'],
                    'beta': ['&beta;',''],
                    'Dmin': ['D<sub>min</sub>','&mu;m'],
                    'q':    ['q',''],
                    'ldisk_lstar':['L<sub>disk</sub>/L<sub>star</sub>',''],
                    'rdisk_bb':['R<sub>BB</sub>','au']}
    
        text = ()
        for res in self.main_results:

            string = ''
            i = 0
            for par in text_ord:
                if par in res.keys():
                    unc,meas = utils.rnd1sf([res['e_'+par],res[par]])
                    if i > 0:
                        string += ' , '
                    string += '{} = {:g} &plusmn; {:g} {}'.format(text_sub[par][0],meas,unc,text_sub[par][1])
                    i += 1

            text = text + (string,)
                
        return text
    
    
    def basic_results(self):
        """A dictionary of basic results, unreliant on classes."""

        r = {}

        # some info
        r['id'] = self.id
        r['model_comps'] = self.model_comps
        r['main_results'] = self.main_results
        r['parameters'] = self.parameters
        r['best_params'] = self.best_params
        r['best_params_1sig'] = self.best_params_1sig
        r['chisq'] = self.chisq

        # observed photometry
        r['phot_wavelength'] = []
        r['phot_fnujy'] = []
        r['phot_e_fnujy'] = []
        r['phot_upperlim'] = []
        r['phot_ignore'] = []
        for p in self.obs:
            if not isinstance(p,photometry.Photometry):
                continue
            r['phot_wavelength'].append(p.mean_wavelength().tolist())
            r['phot_fnujy'].append(p.fnujy.tolist())
            r['phot_e_fnujy'].append(p.e_fnujy.tolist())
            r['phot_upperlim'].append(p.upperlim.tolist())
            r['phot_ignore'].append(p.ignore.tolist())

        # observed spectra
        ispec = -1
        r['spectra'] = []
        for s in self.obs:
            if not isinstance(s,spectrum.ObsSpectrum):
                continue
            t = {}
            t['wavelength'] = s.wavelength.tolist()
            t['fnujy'] = (s.fnujy * self.best_params[ispec]).tolist()
            t['e_fnujy'] = (s.e_fnujy * self.best_params[ispec]).tolist()
            ispec -= 1
            r['spectra'].append(t)

        # spectra of each model component
        r['model_spectra'] = []
        for s in self.comp_spectra:
            r['model_spectra'].append({'wavelength': s.wavelength.tolist(),
                                      'fnujy': s.fnujy.tolist()})

        # total spectra for "star" and "disk" components
        if self.star_spec is not None:
            r['star_spec'] = {}
            r['star_spec']['wavelength'] = self.star_spec.wavelength.tolist()
            r['star_spec']['fnujy'] = self.star_spec.fnujy.tolist()
        if self.disk_spec is not None:
            r['disk_spec'] = {}
            r['disk_spec']['wavelength'] = self.disk_spec.wavelength.tolist()
            r['disk_spec']['fnujy'] = self.disk_spec.fnujy.tolist()

        return r


    def pickle_output(self):
        """Pickle the results for later."""

        # see if we need to write
        if hasattr(self,'pickle_time'):
            if self.analysis_time < self.pickle_time:
                return

        # save for later in a pickle, updating the pickle_time to now
        self.pickle_time = time.time()
        with open(self.pickle,'wb') as f:
            pickle.dump(self,f)


    def load(file):
        """Load a previously save result pickle."""

        with open(file,'rb') as f:
            r = pickle.load(f)

        # reload the necessary models (that weren't saved)
        r.models,r.pl_models = model.get_models(r.obs,r.model_comps)

        return r


    def sed_thumbnail(self,file=None,update=False):
        """Generate a thumbnail of the SED."""
    
        if file is None:
            file = self.sed_thumb
    
        if os.path.exists(file) and not update:
            if os.path.getmtime(file) > self.pickle_time:
                return

        plotting.hardcopy_sed(self,file=file,axis_labels=False,dpi=40)


    def write_json(self,update=False):
        """Write basic results as a json file."""

        if os.path.exists(self.json) and not update:
            if os.path.getmtime(self.json) > self.pickle_time:
                return

        with open(self.json,'w') as f:
            json.dump(self.basic_results(),f)


class FixedResult(BaseResult):
    """Class for unfitted (i.e. specified) models."""

    def __init__(self,rawphot,model_comps,parameters,
                 nospec=False):
        """Do everything at initialisation.
        
        Parameters
        ----------
        rawphot : string
            Name of rawphot file.
        model_comps : tuple of strings
            Models to use.
        parameters : tuple of lists
            Parameters for each model component, excluding normalisation.
        nospec : bool, optional
            Exclude spectra from fitting.
        """
        
        self.file_info(rawphot,model_comps)
        self.fill_data_models(nospec=nospec)
        
        # reorganise models for zero dimensions
        mod = ()
        for mod1,param in zip(self.models,parameters):
            mod_tmp = ()
            for mod_comp in mod1:
                mod_tmp += (model.reduce_zerod(mod_comp,param),)
            mod += (mod_tmp,)

        pl_mod = ()
        for pl1,param in zip(self.pl_models,parameters):
            pl_tmp = ()
            for pl_comp in pl1:
                pl_tmp += (model.reduce_zerod(pl_comp,param),)
            pl_mod += (pl_tmp,)

        self.models = mod
        self.pl_models = pl_mod
        self.model_info = model.models_info(self.models)

        # find best fit normalisation
        res = minimize(fitting.chisq,
                       np.ones(self.model_info['ndim']),
                       args=(self.obs,self.models))
        
        # fill model parameters
        self.parameters = self.model_info['parameters']
        self.best_params = res['x']
        self.comp_best_params = ()
        for param in self.best_params:
            self.comp_best_params += ([param],)

        # fill with post-processing steps
        obs_nel = self.fill_observations()
    
        # get model fluxes and spectra
        self.model_fnujy,model_comp_fnujy = \
            model.model_fluxes(self.models,self.best_params,obs_nel)
        self.residuals,_,_ = fitting.residual(self.best_params,
                                              self.obs,self.models)
        self.fill_best_fit_spectra()


class SampledResult(BaseResult):
    """Class for results from Monte-Carlo or other sampling methods."""

    def star_results(self):
        """Return tuple of dicts of star-specifics, if result has star."""

        star = ()
        distributions = ()
        for i,comp in enumerate(self.model_comps):
            if comp in cfg.models['star']:
                star_one, dist_one = self.star_results_one(i)
                star = star + (star_one,)
                distributions = distributions + (dist_one,)

        return star,distributions


    def star_results_one(self,i):
        """Return dict of star-specifics for ith model component.
        
        .. todo: make stellar parameters consistent, e.g. re-compute
                 lstar from rstar and teff if we have a distance
        """

        star = {}
        star['comp_no'] = i
        distributions = {}
        for j,par in enumerate(self.comp_parameters[i]):
            star[par] = self.best_params[j]
            star['e_'+par] = self.best_params_1sig[j]
        
        # rstar and lstar if star were at 1pc
        rstar_1pc_dist = np.zeros(cfg.fitting['n_samples'])
        lstar_1pc_dist = np.zeros(cfg.fitting['n_samples'])
        for j,par in enumerate(self.comp_param_samples[i]):
            rstar_1pc_dist[j] = np.sqrt(cfg.ssr * 10**par[-1]/np.pi) \
                                * u.pc.to(u.m)
            lstar_1pc_dist[j] = 4 * np.pi * rstar_1pc_dist[j]**2 \
                                * par[0]**4 * 5.670373e-08 / u.L_sun.to(u.W)
        
        distributions['lstar_1pc'] = lstar_1pc_dist
        self.distributions['lstar_1pc_tot'] += lstar_1pc_dist
        lo,star['lstar_1pc'],hi = fitting.pmn_pc(self.param_sample_probs,
                                                 lstar_1pc_dist,[16.0,50.0,84.0])
        star['e_lstar_1pc_lo'] = star['lstar_1pc'] - lo
        star['e_lstar_1pc_hi'] = hi - star['lstar_1pc']
        star['e_lstar_1pc'] = (star['e_lstar_1pc_lo']+star['e_lstar_1pc_hi'])/2.0
    
        # distance-dependent params
        if 'parallax' in self.distributions.keys():

            star['plx_arcsec'] = self.obs_keywords['plx_value'] / 1e3
            star['e_plx_arcsec'] = self.obs_keywords['plx_err'] / 1e3

            # combine lstar_1pc and plx distributions for lstar
            lstar_dist = lstar_1pc_dist / self.distributions['parallax']**2
            distributions['lstar'] = lstar_dist
            lo,star['lstar'],hi = fitting.pmn_pc(self.param_sample_probs,
                                                 lstar_dist,[16.0,50.0,84.0])
            star['e_lstar_lo'] = star['lstar'] - lo
            star['e_lstar_hi'] = hi - star['lstar']
            star['e_lstar'] = (star['e_lstar_lo']+star['e_lstar_hi'])/2.0
     
            rstar_dist = rstar_1pc_dist / self.distributions['parallax'] / u.R_sun.to(u.m)
            distributions['rstar'] = rstar_dist
            lo,star['rstar'],hi = fitting.pmn_pc(self.param_sample_probs,
                                                 rstar_dist,[16.0,50.0,84.0])
            star['e_rstar_lo'] = star['rstar'] - lo
            star['e_rstar_hi'] = hi - star['rstar']
            star['e_rstar'] = (star['e_rstar_lo']+star['e_rstar_hi'])/2.0
                
        return (star,distributions)


    def disk_r_results(self):
        """Return tuple of dicts of disk-specifics, if result has disk_r."""

        disk_r = ()
        distributions = ()
        for i,comp in enumerate(self.model_comps):
            if comp in cfg.models['disk_r']:
                disk_r_one,dist_one = self.disk_r_results_one(i)
                disk_r = disk_r + (disk_r_one,)
                distributions = distributions + (dist_one,)

        return disk_r,distributions
    
    
    def disk_r_results_one(self,i):
        """Return dict of disk_r-specifics for ith model component."""

        disk_r = {}
        disk_r['comp_no'] = i
        distributions = {}
        for j,par in enumerate(self.comp_parameters[i]):
            if 'log_' in par:
                par_in = par.replace('log_','')
                disk_r[par_in] = 10**self.comp_best_params[i][j]
                disk_r['e_'+par_in] = 10**self.comp_best_params_1sig[i][j]
            else:
                par_in = par
                disk_r[par_in] = self.comp_best_params[i][j]
                disk_r['e_'+par_in] = self.comp_best_params_1sig[i][j]
    
            # array of disk temperature samples
            if 'Temp' in par:
                distributions['tdisk'] = np.zeros(cfg.fitting['n_samples'])
                for k,sample in enumerate(self.comp_param_samples[i]):
                    if par == 'log_Temp':
                        distributions['tdisk'][k] = 10**sample[j]
                    elif par == 'Temp':
                        distributions['tdisk'][k] = sample[j]

        # disk and fractional luminosity
        ldisk_1pc_dist = np.zeros(cfg.fitting['n_samples'])
        for j,par in enumerate(self.comp_param_samples[i]):
            
            # there will only be one SpecModel in the ith component
            for m in self.pl_models[i]:
                if not isinstance(m,model.SpecModel):
                    continue
                s = spectrum.ObsSpectrum(wavelength=m.wavelength,
                                         fnujy=m.fnujy(par))
                s.fill_irradiance()

            ldisk_1pc_dist[j] = s.irradiance \
                        * 4 * np.pi * (u.pc.to(u.m))**2 / u.L_sun.to(u.W)
        
        distributions['ldisk_1pc'] = ldisk_1pc_dist
        lo,disk_r['ldisk_1pc'],hi = fitting.pmn_pc(self.param_sample_probs,
                                                   ldisk_1pc_dist,[16.0,50.0,84.0])
        disk_r['e_ldisk_1pc_lo'] = disk_r['ldisk_1pc'] - lo
        disk_r['e_ldisk_1pc_hi'] = hi - disk_r['ldisk_1pc']
        disk_r['e_ldisk_1pc'] = (disk_r['e_ldisk_1pc_lo']+disk_r['e_ldisk_1pc_hi'])/2.0
        
        # stellar luminosities (if >1 star) were summed already
        if np.sum(self.distributions['lstar_1pc_tot']) > 0.0:

            ldisk_lstar_dist = ldisk_1pc_dist / self.distributions['lstar_1pc_tot']
            distributions['ldisk_lstar'] = ldisk_lstar_dist
            lo,disk_r['ldisk_lstar'],hi = fitting.pmn_pc(self.param_sample_probs,
                                                         ldisk_lstar_dist,[16.0,50.0,84.0])
            disk_r['e_ldisk_lstar_lo'] = disk_r['ldisk_lstar'] - lo
            disk_r['e_ldisk_lstar_hi'] = hi - disk_r['ldisk_lstar']
            disk_r['e_ldisk_lstar'] = (disk_r['e_ldisk_lstar_lo']+
                                       disk_r['e_ldisk_lstar_hi'])/2.0

            # distance (and stellar L)-dependent params
            if 'parallax' in self.distributions.keys():
                lstar = self.distributions['lstar_1pc_tot'] /   \
                        self.distributions['parallax']**2
                rdisk_bb_dist = lstar**0.5 * (278.3/distributions['tdisk'])**2

                distributions['rdisk_bb'] = rdisk_bb_dist
                lo,disk_r['rdisk_bb'],hi = fitting.pmn_pc(self.param_sample_probs,
                                                          rdisk_bb_dist,[16.0,50.0,84.0])
                disk_r['e_rdisk_bb_lo'] = disk_r['rdisk_bb'] - lo
                disk_r['e_rdisk_bb_hi'] = hi - disk_r['rdisk_bb']
                disk_r['e_rdisk_bb'] = (disk_r['e_rdisk_bb_lo']+disk_r['e_rdisk_bb_hi'])/2.0

        return (disk_r,distributions)


class Result(SampledResult):
    """Class to compute and handle multinest results.
    
    .. todo:: change this to MultinestResult, doing so will cause issues
              with loaded pickles, which will still be Result.
    """

    @lru_cache(maxsize=128)
    def get(rawphot,model_comps,update_mn=False,
            update_an=False,update_json=False,update_thumb=False,
            nospec=False):
        """Take photometry file and model_name, and fill the rest.
        
        The process works on a heirarchy of update times, each of which
        must be more recent than the previous, or a redo will be forced.

        rawphot_time < mn_time < mn_a_time < analysis_time < pickle_time
        
        These must also be more recent than the configuration variables
        for the oldest allowed analysis and multinest times.

        Each step can be called, and based on these the code will be
        (re)run or not. fill_data_models will always be run because the
        models are not saved in the pickle to save space.
        
        Parameters
        ----------
        rawphot : str
            Rawphot file from sdb.
        model_comps : tuple of str
            Tuple of model names to fit.
        udpate_mn : bool, optional
            Force update of mutinest fitting.
        udpate_an : bool, optional
            Force update of post-multinest analysis.
        update_thumb : bool, optional
            Force update of thumbnail plot.
        update_json : bool, optional
            Force update of json file.
        nospec : bool, optional
            Exclude and spectra when observations are read in.
            
        See Also
        --------
        result.fill_data_models
        result.run_multinest
        result.run_analysis
        """

        self = Result(rawphot,model_comps)

        # see if we have a pickle of results already
        if os.path.exists(self.pickle):
            with open(self.pickle,'rb') as f:
                self = pickle.load(f)

        # update object with local file info since processing may have
        # been done elsewhere
        self.file_info(rawphot,model_comps)

        # see if we can skip everything except the json
        if hasattr(self,'rawphot_time') and hasattr(self,'mn_time') and \
            hasattr(self,'mn_a_time') and hasattr(self,'mn_time') and   \
            hasattr(self,'rawphot_time') and hasattr(self,'pickle_time'):
            if not update_an and not update_mn and \
                    self.exclude_spectra == nospec and \
                    self.mn_time > cfg.fitting['mn_oldest'] and \
                    self.mn_a_time > cfg.fitting['an_oldest'] and \
                    self.pickle_time > self.analysis_time > self.mn_a_time > \
                    self.mn_time > self.rawphot_time:

                self.sed_thumbnail(update=update_thumb)
                self.write_json(update=update_json)
                return self

        mn_up = update_mn
        if hasattr(self,'exclude_spectra'):
            mn_up = (self.exclude_spectra != nospec or mn_up)

        # run everything, these will check if anything needs to be done
        self.fill_data_models(nospec=nospec)
        # stop if no photometry
        if not hasattr(self,'obs'):
            return self
        self.run_multinest(update_mn=mn_up)
        self.run_analysis(update_an=update_an)

        # delete the models to save space, we don't need them again
        self.models = ''
        self.pl_models = ''

        self.pickle_output()
        self.sed_thumbnail(update=update_thumb)
        self.write_json(update=update_json)

        return self


    def run_multinest(self,update_mn=False):
        """Run multinest.
        
        Run multinest fitting, get the results as a pymultinest
        analyzer object, and make a corner plot of the points. The
        fitting will be run if it has never been run (no files or no
        mn_time attribute), is older than the read time of the 
        photometry file, or is older than the mn_oldest time in
        config.fitting.
        
        All multinest files are updated, even if a previous run was
        used, so the file modification times are not useful for testing
        for whether an update is needed.

        Parameters
        ----------
        update_mn : bool, optional
             Force update of multinest fitting.
        """

        # if we want to re-run multinest, delete previous output first
        run_mn = update_mn
        if os.path.exists(self.pmn_base+'phys_live.points'):
            if hasattr(self,'mn_time'):
                if self.rawphot_time > self.mn_time:
                    run_mn = True
                if cfg.fitting['mn_oldest'] > self.mn_time:
                    run_mn = True
            else:
                run_mn = True

        # multinest does checkpointing, so we can force a re-run by
        # deleting the files
        if run_mn:
            self.delete_multinest()
        
        # we must go there, multinest only takes 100 char paths
        with utils.pushd(self.pmn_dir):
            fitting.multinest( self.obs,self.models,'.' )
            self.mn_time = os.path.getmtime(self.pmn_base+'phys_live.points')

        # update the analyzer if necessary
        get_a = False
        if hasattr(self,'mn_a_time'):
            if self.mn_time > self.mn_a_time:
                get_a = True
        else:
            get_a = True

        if get_a:
            a = pmn.Analyzer(outputfiles_basename=self.pmn_base,
                             n_params=self.model_info['ndim'])
            self.analyzer = a
            self.mn_a_time = time.time()

        # parameter fitting corner plot
        plot = False
        if not os.path.exists(self.corner_plot):
            plot = True
        else:
            if os.path.getmtime(self.corner_plot) < self.mn_a_time:
                plot = True
            
        if plot:
            d = self.analyzer.get_data()
            fig = corner.corner(d[:,2:],labels=self.model_info['parameters'])
            fig.savefig(self.corner_plot)
            plt.close(fig) # not doing this causes an epic memory leak

        # parameter names and best fit
        self.evidence = self.analyzer.get_stats()['global evidence']
        self.parameters = self.model_info['parameters']

        self.best_params = []
        self.best_params_1sig = []
        for i in range(len(self.parameters)):
            self.best_params.append(self.analyzer.get_stats()\
                                    ['marginals'][i]['median'])
            self.best_params_1sig.append(self.analyzer.get_stats()\
                                         ['marginals'][i]['sigma'])
        
        # tuple of multinest samples to use for uncertainty estimation
        self.param_samples = ()
        self.param_sample_probs = []
        randi = []
        for i in np.random.randint(0,high=len(self.analyzer.data),
                                   size=cfg.fitting['n_samples']):
            randi.append(i)
            self.param_samples += (self.analyzer.data[i,2:],)
            self.param_sample_probs.append(self.analyzer.data[i,0])

        # split the parameters into components
        self.n_parameters = len(self.parameters)
        self.comp_best_params = ()
        self.comp_best_params_1sig = ()
        self.comp_parameters = ()
        self.comp_param_samples = ()
        i0 = 0
        for comp in self.models:
            nparam = len(comp[0].parameters)+1
            self.comp_parameters += (comp[0].parameters,)
            self.comp_best_params += (self.best_params[i0:i0+nparam],)
            self.comp_best_params_1sig += (self.best_params_1sig[i0:i0+nparam],)
            
            comp_i_samples = ()
            for i in randi:
                comp_i_samples += (self.analyzer.data[i,i0+2:i0+nparam+2],)
            self.comp_param_samples += (comp_i_samples,)

            i0 += nparam
        

    def run_analysis(self,update_an=False):
        """Run analysis of the multinest results.
        
        Run post-multinest fitting analysis. Will be run if forced to,
        or if the multinest analyzer time is more recent than the
        previous analysis time.
        
        Parameters
        ----------
        update_an : bool, optional
            Force update of analysis
        """

        # see if we have anything to do
        run_an = update_an
        if hasattr(self,'analysis_time'):
            if self.mn_a_time > self.analysis_time:
                run_an = True
        else:
            run_an = True
        if not run_an:
            return

        # fluxes and uncertainties etc. using parameter samples
        self.distributions = {}
        
        # we will want lstar at 1pc below
        self.distributions['lstar_1pc_tot'] = np.zeros(cfg.fitting['n_samples'])
        
        # generate a normal distribution of parallaxes, truncated to
        # contain no negative values, if there is an uncertainty
        if self.obs_keywords['plx_err'] is not None \
            and self.obs_keywords['plx_value'] is not None:
            if self.obs_keywords['plx_err'] > 0 \
                and self.obs_keywords['plx_value'] > 0:
                
                lo_cut = -1. * ( self.obs_keywords['plx_value'] /   \
                                 self.obs_keywords['plx_err'] )
                                 
                self.distributions['parallax'] = \
                    truncnorm.rvs(lo_cut,np.inf,
                                  loc=self.obs_keywords['plx_value']/1e3,
                                  scale=self.obs_keywords['plx_err']/1e3,
                                  size=cfg.fitting['n_samples'])
                    
        # observed fluxes
        obs_nel = self.fill_observations()

        # model fluxes, including colours/indices
        model_dist = np.zeros((len(self.filters),cfg.fitting['n_samples']))
        model_comp_dist = np.zeros((self.n_comps,len(self.filters),
                                    cfg.fitting['n_samples']))

        # all fluxes, including colours/indices
        p_all = photometry.Photometry(filters=filter.Filter.all)
        self.all_filters = p_all.filters
        p_all_mod,_ = model.get_models((p_all,),self.model_comps)
        all_dist = np.zeros((p_all.nphot,cfg.fitting['n_samples']))
        all_comp_dist = np.zeros((self.n_comps,p_all.nphot,
                                    cfg.fitting['n_samples']))

        for i,par in enumerate(self.param_samples):

            # TODO: compute once, grab model fluxes from all fluxes...
            model_fnujy,model_comp_fnujy = \
                model.model_fluxes(self.models,par,obs_nel)
            model_dist[:,i] = model_fnujy
            model_comp_dist[:,:,i] = model_comp_fnujy

            all_fnujy,all_comp_fnujy = \
                model.model_fluxes(p_all_mod,par,[p_all.nphot])
            all_dist[:,i] = all_fnujy
            all_comp_dist[:,:,i] = all_comp_fnujy

        # summed model fluxes
        self.distributions['model_fnujy'] = model_dist
        lo,self.model_fnujy,hi = fitting.pmn_pc(self.param_sample_probs,
                                                model_dist,[16.0,50.0,84.0],
                                                axis=1)
        self.model_fnujy_1sig_lo = self.model_fnujy - lo
        self.model_fnujy_1sig_hi = hi - self.model_fnujy

        # per-component model fluxes
        self.distributions['model_comp_fnujy'] = model_comp_dist
        lo,self.model_comp_fnujy,hi = fitting.pmn_pc(self.param_sample_probs,
                                                     model_comp_dist,[16.0,50.0,84.0],
                                                     axis=2)
        self.model_comp_fnujy_1sig_lo = self.model_comp_fnujy - lo
        self.model_comp_fnujy_1sig_hi = hi - self.model_comp_fnujy

        # residuals and fitting results
        self.residuals,_,_ = fitting.residual(self.best_params,
                                              self.obs,self.models)
        self.chisq = np.sum( np.square( self.residuals ) )
        self.dof = len(self.wavelengths)-len(self.parameters)-1

        # star and disk photometry for all filters
        star_phot_dist = np.zeros((p_all.nphot,cfg.fitting['n_samples']))
        disk_phot_dist = np.zeros((p_all.nphot,cfg.fitting['n_samples']))
        for i,comp in enumerate(self.model_comps):
            if self.star_or_disk[i] == 'star':
                star_phot_dist += all_comp_dist[i,:,:]
            elif self.star_or_disk[i] == 'disk':
                disk_phot_dist += all_comp_dist[i,:,:]

        # star photometry in all filters
        self.distributions['star_phot'] = star_phot_dist
        lo,self.all_star_phot,hi = fitting.pmn_pc(self.param_sample_probs,
                                              star_phot_dist,[16.0,50.0,84.0],
                                              axis=1)
        self.all_star_phot_1sig_lo = self.all_star_phot - lo
        self.all_star_phot_1sig_hi = hi - self.all_star_phot

        if np.sum(self.all_star_phot) == 0:
            self.all_star_phot = None

        # disk photometry in all filters
        self.distributions['disk_phot'] = disk_phot_dist
        lo,self.all_disk_phot,hi = fitting.pmn_pc(self.param_sample_probs,
                                              disk_phot_dist,[16.0,50.0,84.0],
                                              axis=1)
        self.all_disk_phot_1sig_lo = self.all_disk_phot - lo
        self.all_disk_phot_1sig_hi = hi - self.all_disk_phot

        if np.sum(self.all_disk_phot) == 0:
            self.all_disk_phot = None

        # component photometry in all filters
        self.distributions['all_comp_phot'] = all_comp_dist
        lo,self.all_comp_phot,hi = fitting.pmn_pc(self.param_sample_probs,
                                                  all_comp_dist,[16.0,50.0,84.0],
                                                  axis=2)
        self.all_comp_phot_1sig_lo = self.all_comp_phot - lo
        self.all_comp_phot_1sig_hi = hi - self.all_comp_phot

        # total photometry in all filters
        self.distributions['all_phot'] = all_dist
        lo,self.all_phot,hi = fitting.pmn_pc(self.param_sample_probs,
                                             all_dist,[16.0,50.0,84.0],
                                             axis=1)
        self.all_phot_1sig_lo = self.all_phot - lo
        self.all_phot_1sig_hi = hi - self.all_phot

        # model-specifics, also combined into a single tuple
        self.star,self.star_distributions = self.star_results()
        self.disk_r,self.disk_r_distributions = self.disk_r_results()
        self.main_results = self.star + self.disk_r
        
        # corner plot of distributions, join them together first
        if not os.path.exists(self.distributions_plot):
            samples = np.zeros(cfg.fitting['n_samples'])
            labels = []
            if 'parallax' in self.distributions.keys():
                samples = np.vstack((samples,self.distributions['parallax']))
                labels.append('parallax')
            for dist in self.star_distributions:
                for key in dist.keys():
                    samples = np.vstack((samples,dist[key]))
                    labels.append(key)
            for dist in self.disk_r_distributions:
                for key in dist.keys():
                    samples = np.vstack((samples,dist[key]))
                    labels.append(key)
            samples = samples[1:]

            fig = corner.corner(samples.transpose(),labels=labels)
            fig.savefig(self.distributions_plot)
            plt.close(fig)

        # set analysis finish time
        self.analysis_time = time.time()


    def delete_multinest(self):
        """Delete multinest output so it can be run again."""

        fs = glob.glob(self.pmn_base+'*')
        for f in fs:
            os.remove(f)


def sort_results(results):
    """Return indices to sort a list of Result objects by evidence."""

    ev = []
    ndim = []
    for r in results:

        ndim.append( r.model_info['ndim'] )
        ev.append( r.evidence )

    return fitting.sort_evidence(ev,ndim)
