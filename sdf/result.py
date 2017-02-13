from functools import lru_cache
import os.path
import pickle
import glob

import numpy as np
import pymultinest as pmn
import matplotlib.pyplot as plt
import corner
import astropy.units as u

from . import photometry
from . import spectrum
from . import model
from . import filter
from . import fitting
from . import utils
from .utils import SdfError
from . import config as cfg


class Result(object):
    """Class to compute and handle multinest results."""

    @lru_cache(maxsize=128)
    def __init__(self,rawphot,model_comps,update=False,nospec=False):
        """Take photometry file and model_name, and fill the rest."""
        
        self.rawphot = rawphot
        self.model_comps = model_comps
        self.star_or_disk = ()
        for comp in model_comps:
            if comp in cfg.models['star']:
                self.star_or_disk += ('star',)
            elif comp in cfg.models['disk']:
                self.star_or_disk += ('disk',)
            else:
                raise SdfError("couldn't assigm comp {} to star or disk "
                               "given lists in {} and {}".
                               format(comp,cfg.models['star'],
                                      cfg.models['disk']))
    
        self.n_comps = len(model_comps)
        
        # id
        self.id = os.path.basename(rawphot).rstrip('-rawphot.txt')

        # where the rawphot file is
        self.path = os.path.dirname(rawphot)
        
        # where the multinest output is (or will be), create if needed
        self.pmn_dir = self.path + '/' + self.id            \
                       + cfg.fitting['pmn_dir_suffix']
        if not os.path.exists(self.pmn_dir):
            os.mkdir(self.pmn_dir)
        
        # the base name for multinest files
        self.pmn_base = self.pmn_dir + '/'                  \
                        + '+'.join(self.model_comps)        \
                        + cfg.fitting['pmn_model_suffix']

        # observations; keywords, tuples of photometry and spectra. if
        # there is nothing in the photometry file then don't fill
        # anything else
        p = photometry.Photometry.read_sdb_file(self.rawphot)
        if p is None:
            return

        self.obs = (p,)
        self.obs_keywords = utils.get_sdb_keywords(self.rawphot)
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

        # pickle of results (file may not actually exist)
        self.pickle = self.pmn_base + '.pkl'
        
        # see if multinest needs to be re-run
        # TODO: pickling saves no time as is, either remove or change
        run = update
        if not os.path.exists(self.pickle):
            run = True
        elif os.path.getmtime(self.pickle) < os.path.getmtime(rawphot):
            run = True

        # if it does, delete the output and run
        if run:
            self.delete_multinest()
            
            # go there, multinest only takes 100 char paths
            with utils.pushd(self.pmn_dir):
                fitting.multinest( self.obs,self.models,'.' )

            a = pmn.Analyzer(outputfiles_basename=self.pmn_base,
                             n_params=self.model_info['ndim'])
            pfile = open(self.pickle,'wb')
            pickle.dump(a,pfile)
            pfile.close()

        # load the results
        try:
            self.analyzer = a
        except NameError:
            pfile = open(self.pickle,'rb')
            self.analyzer = pickle.load(pfile)
            pfile.close()
        
        # when the results were finished
        self.mtime = os.path.getmtime(self.pickle)

        # parameter corner plot if needed
        self.chain_plot = self.pmn_base+'corner.png'
        run = False
        if not os.path.exists(self.chain_plot):
            run = True
        else:
            if os.path.getmtime(self.chain_plot) < self.mtime:
                run = True
            
        if run:
            d = self.analyzer.get_data()
            fig = corner.corner(d[:,2:],show_titles=True,
                                labels=self.model_info['parameters'],
                                quantiles=[0.16, 0.5, 0.84])
            fig.savefig(self.chain_plot)
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
        
        self.n_parameters = len(self.parameters)
        self.comp_best_params = ()
        self.comp_best_params_1sig = ()
        self.comp_parameters = ()
        i0 = 0
        for comp in self.models:
            nparam = len(comp[0].parameters)+1
            self.comp_parameters += (comp[0].parameters,)
            self.comp_best_params += (self.best_params[i0:i0+nparam],)
            self.comp_best_params_1sig += (self.best_params_1sig[i0:i0+nparam],)
            i0 += nparam
        
        # fluxes etc., this is largely copied from fitting.residual
        
        # observed fluxes
        tmp = fitting.concat_obs(self.obs)
        self.obs_fnujy,self.obs_e_fnujy,self.obs_upperlim,self.filters_ignore,\
            obs_ispec,obs_nel,self.wavelengths,self.filters,self.obs_bibcode = tmp
        spec_norm = np.take(self.best_params+[1.0],obs_ispec)
        self.obs_fnujy = self.obs_fnujy * spec_norm
        self.obs_e_fnujy = self.obs_fnujy * spec_norm

        # model photometry and residuals, including filling of
        # colours/indices
        self.model_fnujy,self.model_comp_fnujy                      \
            = model.model_fluxes(self.models,self.best_params,obs_nel)
        self.residuals,_,_ = fitting.residual(self.best_params,
                                              self.obs,self.models)
        self.chisq = np.sum( np.square( self.residuals ) )
        self.dof = len(self.wavelengths)-len(self.parameters)-1

        # star/disk photometry for all filters, obs_nel is just no of
        # filters (p_all.nphot)
        star_comps = ()
        star_params = []
        disk_comps = ()
        disk_params = []
        for i,comp in enumerate(self.model_comps):
            if self.star_or_disk[i] == 'star':
                star_comps += (comp,)
                star_params += self.comp_best_params[i]
            elif self.star_or_disk[i] == 'disk':
                disk_comps += (comp,)
                disk_params += self.comp_best_params[i]
        
        p_all = photometry.Photometry(filters=filter.Filter.all)
        if len(star_comps) > 0:
            star_mod,_ = model.get_models((p_all,),star_comps)
            self.star_phot,_ = model.model_fluxes(star_mod,star_params,[p_all.nphot])
        else:
            self.star_phot = None
        
        if len(disk_comps) > 0:
            disk_mod,_ = model.get_models((p_all,),disk_comps)
            self.disk_phot,_ = model.model_fluxes(disk_mod,disk_params,[p_all.nphot])
        else:
            self.disk_phot = None
        
        self.all_filters = p_all.filters
        
        # total photometry
        if self.star_phot is not None:
            self.all_phot = self.star_phot
            if self.disk_phot is not None:
                self.all_phot += self.disk_phot
        
        else:
            if self.disk_phot is not None:
                self.all_phot = self.disk_phot
            else:
                self.all_phot = None

        # ObsSpectrum for each component
        wave = cfg.models['default_wave']
        star_spec = np.zeros(len(wave))
        disk_spec = np.zeros(len(wave))
        self.comp_spectra = ()
        for i,comp in enumerate(self.pl_models):
            for m in comp:
                if not isinstance(m,model.SpecModel):
                    continue
                s = spectrum.ObsSpectrum(wavelength=m.wavelength,
                                         fnujy=m.fnujy(self.comp_best_params[i]))
                s.fill_irradiance()
                self.comp_spectra += (s,)
    
                mtmp = m.copy()
                mtmp.interp_to_wavelengths(wave)
                if self.star_or_disk[i] == 'star':
                    star_spec += mtmp.fnujy(self.comp_best_params[i])
                elif self.star_or_disk[i] == 'disk':
                    disk_spec += mtmp.fnujy(self.comp_best_params[i])

        # and star/disk spectra
        # TODO: add realistic uncertainties
        if np.max(star_spec) > 0:
            self.star_spec = spectrum.ObsSpectrum(wavelength=wave,fnujy=star_spec)
        else:
            self.star_spec = None

        if np.max(disk_spec) > 0:
            self.disk_spec = spectrum.ObsSpectrum(wavelength=wave,fnujy=disk_spec)
        else:
            self.disk_spec = None

        # model-specifics, also combined into a single tuple
        self.star = self.star_results()
        self.disk_r = self.disk_r_results()
        self.main_results = self.star + self.disk_r


    def star_results(self):
        """Return tuple of dicts of star-specifics, if result has star."""

        star = ()
        for i,comp in enumerate(self.model_comps):
            if comp in cfg.models['star']:
                star = star + (self.star_results_one(i),)

        return star


    def star_results_one(self,i):
        """Return dict of star-specifics for ith model component."""

        star = {}
        for j,par in enumerate(self.comp_parameters[i]):
            star[par] = self.best_params[j]
            star['e_'+par] = self.best_params_1sig[j]
        
        # stellar luminosity at 1pc, uncertainty is normalisation
        frac_norm = np.log(10) * self.comp_best_params_1sig[i][-1]
        star['lstar_1pc'] = self.comp_spectra[i].irradiance \
                    * 4 * np.pi * (u.pc.to(u.m))**2 / u.L_sun.to(u.W)
        star['e_lstar_1pc'] = star['lstar_1pc'] * frac_norm

        # distance-dependent params
        if self.obs_keywords['plx_value'] is not None:
            if self.obs_keywords['plx_value'] > 0:
                plx_arcsec = self.obs_keywords['plx_value'] / 1e3
                
                if self.obs_keywords['plx_err'] is not None:
                    e_plx_arcsec = self.obs_keywords['plx_err'] / 1e3
                else:
                    e_plx_arcsec = plx_arcsec / 3.
                
                star['lstar'] = star['lstar_1pc'] / plx_arcsec**2
                star['e_lstar'] = star['lstar'] * \
                                  np.sqrt( frac_norm**2
                                          + (2*e_plx_arcsec/plx_arcsec)**2 )
                                  
                star['rstar'] = np.sqrt(cfg.ssr *
                                        10**self.comp_best_params[i][-1]/np.pi) \
                        * u.pc.to(u.m) / plx_arcsec / u.R_sun.to(u.m)
                star['e_rstar'] = star['rstar'] * \
                                  np.sqrt( frac_norm**2
                                          + (2*e_plx_arcsec/plx_arcsec)**2)

        return star


    def disk_r_results(self):
        """Return tuple of dicts of disk-specifics, if result has disk_r."""

        disk_r = ()
        for i,comp in enumerate(self.model_comps):
            if comp in cfg.models['disk_r']:
                disk_r= disk_r + (self.disk_r_results_one(i),)

        return disk_r
    
    
    def disk_r_results_one(self,i):
        """Return dict of disk_r-specifics for ith model component."""

        disk_r = {}
        for j,par in enumerate(self.comp_parameters[i]):
            if 'log_' in par:
                par_in = par.replace('log_','')
                disk_r[par_in] = 10**self.comp_best_params[i][j]
                disk_r['e_'+par_in] = 10**self.comp_best_params_1sig[i][j]
            else:
                disk_r[par] = self.comp_best_params[i][j]
                disk_r['e_'+par] = self.comp_best_params_1sig[i][j]

        # disk and fractional luminosity
        frac_norm = np.log(10) * self.comp_best_params_1sig[i][-1]
        disk_r['ldisk_1pc'] = self.comp_spectra[i].irradiance \
                    * 4 * np.pi * (u.pc.to(u.m))**2 / u.L_sun.to(u.W)
        disk_r['e_ldisk_1pc'] = disk_r['ldisk_1pc'] * frac_norm

        # sum stellar luminosity
        lstar_1pc = 0.0
        e_lstar_1pc = 0.0
        if isinstance(self.star,tuple):
            for star in self.star:
                lstar_1pc += star['lstar_1pc']
                e_lstar_1pc = np.sqrt(e_lstar_1pc**2 + star['e_lstar_1pc']**2)

        if lstar_1pc > 0.0:
            frac_lstar_1pc = e_lstar_1pc / lstar_1pc
            disk_r['ldisk_lstar'] = disk_r['ldisk_1pc']/lstar_1pc
            disk_r['e_ldisk_lstar'] = disk_r['ldisk_lstar'] * \
                            np.sqrt(frac_norm**2 + frac_lstar_1pc**2)

            # distance (and stellar L)-dependent params
            if self.obs_keywords['plx_value'] is not None:
                plx_arcsec = self.obs_keywords['plx_value'] / 1e3
                disk_r['rdisk_bb'] = (lstar_1pc/plx_arcsec**2)**0.5 * \
                                        (278.3/disk_r['Temp'])**2
                disk_r['e_rdisk_bb'] = disk_r['rdisk_bb'] * \
                            np.sqrt( (0.5*frac_lstar_1pc)**2
                                    + (2*(disk_r['e_Temp']/disk_r['Temp']))**2 )

        return disk_r


    def main_results_text(self):
        """Return nicely formatted tuple of text of results."""
    
        # the array sets the order, and the dict the conversion
        text_ord = ['Teff','lstar','rstar',
                    'Temp','ldisk_lstar','rdisk_bb','lam0','beta']
        text_sub = {'Teff':['T<sub>star</sub>','K'],
                    'MH':['[M/H]',''],
                    'logg':['logg',''],
                    'lstar':['L<sub>star</sub>','L<sub>Sun</sub>'],
                    'rstar':['R<sub>star</sub>','R<sub>Sun</Sun>'],
                    'Temp':['T<sub>dust</sub>','K'],
                    'lam0':['&lambda;<sub>0</sub>','&mu;m'],
                    'beta':['&beta;',''],
                    'ldisk_lstar':['L<sub>disk</sub>/L<sub>star</sub>',''],
                    'rdisk_bb':['R<sub>BB</sub>','au']}
    
        text = ()
        for res in self.main_results:

            string = ''
            for par in text_ord:
                if par in res.keys():
                    unc,meas = utils.rnd1sf([res['e_'+par],res[par]])
                    string += '{} = {:g} &plusmn; {:g} {} , '.format(text_sub[par][0],meas,unc,text_sub[par][1])

            text = text + (string,)
                
        return text
    
    
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
