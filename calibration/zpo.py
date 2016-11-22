import os
import sys
sys.path.insert(0, os.path.abspath('..'))

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

# get the info for the targets
t1 = Table.read('gaia_benchmark/gb.csv')
t2 = Table.read('gaia_benchmark/ids.csv')
t = operations.join(t1,t2)

# grab the photometry and the models
obs = []
mod_fnujy = []
mod_filt = []
allmod = model.PhotModel.read_model('kurucz_m')
for i,id in enumerate(t['sdbid']):

    obs.append( photometry.Photometry.read_sdb_file('gaia_benchmark/'+id+'-rawphot.txt') )
    tmp = allmod.copy()
    tmp.keep_filters( obs[-1].filters )
    fnujy = tmp.fnujy( [t[i]['Teff'],t[i]['logg'],t[i]['[Fe/H]'],0.0] )
    norm = np.median(obs[-1].fnujy/fnujy)
    mod_fnujy.append( tmp.fnujy( [t[i]['Teff'],t[i]['logg'],t[i]['[Fe/H]'],np.log10(norm)] ) )
    mod_filt.append( filter.mean_wavelength(tmp.filters) )

    plt.semilogx( mod_filt[-1], (obs[-1].fnujy-mod_fnujy[-1])/obs[-1].fnujy + i,'.' )

plt.show()
