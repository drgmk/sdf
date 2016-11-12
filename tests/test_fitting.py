import pytest
import numpy as np
import astropy.units as u

from .context import sdf

def test_fitting_fill_colours():
    m = sdf.model.PhotModel.read_model('kurucz')
    fs = ['BS_YS','MIPS24']
    m.keep_filters(fs,colour_bases=True)
    mod_fnu = sdf.model.model_fluxes( (m,), [6300,4.5,-3], [2] )
