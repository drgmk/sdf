import os
import sys
sys.path.insert(0, os.path.abspath('..'))

from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table,operations

from sdf import photometry
from sdf import model
from sdf import spectrum
from sdf import plotting
from sdf import fitting
from sdf import filter
from sdf import config as cfg

def chisq_norm(norm,*args):
    """Quick chisq for only norm."""
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
