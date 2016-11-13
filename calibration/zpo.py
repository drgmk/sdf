import os
import sys
sys.path.insert(0, os.path.abspath('..'))

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
