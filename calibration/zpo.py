import os
import sys
sys.path.insert(0, os.path.abspath('..'))

#from scipy.optimize import minimize
import numpy as np
#import matplotlib.pyplot as plt
#from astropy.table import Table,operations

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
    zpo = 2.5 * np.log10( zpo_flux ) + 0.027 - 0.04777335


    print('MvB: ',fnu)
    print('Vega:',f_vega.zero_point)
    print('ZPO: ',zpo)


"""This was supposed to be a way to set zero point offsets using some
    set of reference start, in particular Gaia benchmarks, but it didn't
    look like it was going to work very well so has been shelved
    
def chisq_norm(norm,*args):
    par0,obs,mod = args
    return fitting.chisq( par0+[norm],(obs,), (mod,) )

# get the info for the targets
t1 = Table.read('gaia_benchmark/gb.csv')
t2 = Table.read('gaia_benchmark/ids.csv')
t = operations.join(t1,t2)

# grab the photometry and the models
obs = []
mod_fnujy = []
mod_norm = np.array([])
mod_filt = []
all_filt = np.array([])
allmod = model.PhotModel.read_model('phoenix_m')
for i,id in enumerate(t['sdbid']):

    obs.append( photometry.Photometry.read_sdb_file('gaia_benchmark/'+id+'-rawphot.txt') )
    
    tmp = allmod.copy()
    tmp.keep_filters( obs[-1].filters )
    par0 = [t[i]['Teff'],t[i]['logg'],t[i]['[Fe/H]']]
    fnujy = tmp.fnujy( par0+[-0.5] )
    mod_filt.append( filter.mean_wavelength(tmp.filters) )
    all_filt = np.append( all_filt, tmp.filters )

    res = minimize(chisq_norm, np.log10( np.median(obs[-1].fnujy/fnujy) ),
                   args=(par0,obs[-1],tmp))
    mod_norm = np.append(mod_norm,res['x'][0])
    mod_fnujy.append( tmp.fnujy( par0+[res['x'][0]] ) )

#    plt.loglog(mod_filt[-1],obs[-1].fnujy)
#    plt.loglog(mod_filt[-1],mod_fnujy[-1])
#    plt.show()

    plt.semilogx(mod_filt[-1],(obs[-1].fnujy/mod_fnujy[-1]) + 2*i,'.')

all_filt = np.unique(all_filt)
plt.show()
"""
