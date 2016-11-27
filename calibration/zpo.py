import os
import sys
sys.path.insert(0, os.path.abspath('..'))

from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.table import Table,operations

from sdf import photometry
from sdf import model
from sdf import spectrum
from sdf import plotting
from sdf import fitting
from sdf import filter
from sdf import config as cfg

def mvb_zpo(fname,flam=None):
    """Get the Mann & von Braun ZPO for the CALSPEC Vega spectrum.
        
    Values are scaled so that the ZPO at Johnson V is 0.027.
    """

    f_vega = filter.Filter.get(fname)

    # this assumes the SVO has the correct zero points
    svn = filter.filter_info.filters[fname]['svo_name']
    f_svo = filter.Filter.svo_get(svn)
    zpo_flux = f_svo.zero_point / f_vega.zero_point
    zpo = 2.5 * np.log10( zpo_flux ) + 0.027 - 0.064605028
    print(fname)
    print('SVO: ',f_svo.zero_point)
    print('Vega:',f_vega.zero_point)
    print('ZPO: ',zpo)
    print('')
    
    # this uses the pivot wavelength to convert F_lam
    fnu = flam * u.Unit('erg/(cm**2*s*AA)').\
            to('Jy',equivalencies=u.spectral_density(f_vega.pivot_wavelength()*u.micron))
    
    zpo_flux = fnu / f_vega.zero_point
    zpo = 2.5 * np.log10( zpo_flux ) #+ 0.027 - 0.04777335


    print('MvB: ',fnu)
    print('Vega:',f_vega.zero_point)
    print('ZPO: ',zpo)


"""This was supposed to be a way to set zero point offsets using some
    set of reference stars, in particular Gaia benchmarks, but it didn't
    look like it was going to work very well so has been shelved
    """
    
def chisq_one(obs,mod,par,mult,mfilt):
    fac = np.array([])
    for f in obs.filters:
        fac = np.append(fac,mult[f==mfilt])

    res = minimize( fitting.chisq,par,args=( (obs,), (mod,) ) )
    return res['x'],res['fun']

# get the info for the targets
t1 = Table.read('calibration/gaia_benchmark/gb.csv')
t2 = Table.read('calibration/gaia_benchmark/ids.csv')
t = operations.join(t1,t2)
t = t[0:3]

# grab the photometry and the models
obs = []
mod = []
par = []
allfilt = []
for i,id in enumerate(t['sdbid']):

    obs.append( photometry.Photometry.read_sdb_file('calibration/gaia_benchmark/'+id+'-rawphot.txt') )
    _,tmp = model.get_models((obs[-1],),('phoenix_m',))
    mod.append(tmp)

    fnujy = tmp[0][0].fnujy( par0+[-0.5] )
    par.append( np.array( [t[i]['Teff'],t[i]['logg'],t[i]['[Fe/H]'],
                           np.log10( np.median(obs[-1].fnujy/fnujy) )] ) )
    allfilt += obs[-1].filters.tolist()
    
# set up filter multipliers
zp_filt = np.unique(allfilt)
zp_fac = np.ones(len(zp_filt))
               
par_tmp,chisq = chisq_one(obs[-1],mod[-1],par[-1],zp_fac,zp_filt)
               
plt.semilogx(obs[-1].fnujy/mod[-1],fnujy(par) + 2*i,'.')

plt.show()
