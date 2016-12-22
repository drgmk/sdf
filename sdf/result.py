from functools import lru_cache
import os.path
import pickle
import glob

import numpy as np
import pymultinest as pmn
import matplotlib.pyplot as plt
import corner

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

        # parameter corner plot
        self.chain_plot = self.pmn_base+'corner.png'
        if not os.path.exists(self.chain_plot):
            d = self.analyzer.get_data()
            fig = corner.corner(d[:,2:])
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
        
        # path to SED plot (which probably doesn't exist yet)
        self.sed_plot = self.path + '/index.html'


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
