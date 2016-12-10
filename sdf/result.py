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
from . import fitting
from . import utils
from . import config as cfg


class Result(object):
    """Class to compute and handle multinest results."""

    @lru_cache(maxsize=128)
    def __init__(self,rawphot,model_comps,update=False,nospec=False):
        """Take photometry file and model_name, and fill the rest."""
        
        self.rawphot = rawphot
        self.model_comps = model_comps
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
                                                   module_split=True)
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
        self.best_params = self.analyzer.get_best_fit()['parameters']
        self.n_parameters = len(self.parameters)
        self.comp_best_params = ()
        self.comp_parameters = ()
        i0 = 0
        for comp in self.models:
            nparam = len(comp[0].parameters)+1
            self.comp_parameters += (comp[0].parameters,)
            self.comp_best_params += (self.best_params[i0:nparam],)
            i0 += nparam
        
        # photometry etc., this is largely copied from fitting.residual
        
        # concatenated fluxes
        tmp = fitting.concat_obs(self.obs)
        self.obs_fnujy,self.obs_e_fnujy,self.obs_upperlim,self.filters_ignore,\
            obs_ispec,obs_nel,self.wavelengths,self.filters,self.obs_bibcode = tmp

        spec_norm = np.take(self.best_params+[1.0],obs_ispec)
        self.obs_fnujy = self.obs_fnujy * spec_norm
        self.obs_e_fnujy = self.obs_fnujy * spec_norm

        # get model fluxes, including filling of colours/indices
        self.model_fnujy,self.model_comp_fnujy                      \
            = model.model_fluxes(self.models,self.best_params,obs_nel)
        
        # residuals, a bit neater to use fitting.residual
        self.residuals,_,_ = fitting.residual(self.best_params,
                                              self.obs,self.models)
        self.chisq = np.sum( np.square( self.residuals ) )
        self.dof = len(self.wavelengths)-len(self.parameters)-1

        # ObsSpectrum for each component
        self.comp_spectra = ()
        for i,comp in enumerate(self.pl_models):
            for m in comp:
                if not isinstance(m,model.SpecModel):
                    continue
                s = spectrum.ObsSpectrum(wavelength=m.wavelength,
                                         fnujy=m.fnujy(self.comp_best_params[i]))
                s.fill_irradiance()
                self.comp_spectra += (s,)

        # path to SED plot (which probably doesn't exist yet)
        self.sed_plot = self.path + '/' + self.id + cfg.pl['sed_suffix']


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
