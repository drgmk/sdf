import pytest
import numpy as np
import astropy.units as u

from .context import sdf

def test_model_keep_filters():
    m = sdf.model.PhotModel.read_model('kurucz')
    fs = ['BS_YS','MIPS24']
    m1 = m.copy()
    m1.keep_filters(fs,colour_bases=True)
    assert( np.all( np.equal(m1.filters,['BS_YS','MIPS24','BS','YS']) ) )
    m.keep_filters(fs,colour_bases=False)
    assert( np.all( np.equal(m1.filters,['BS_YS','MIPS24']) ) )

def test_read_kurucz():
    ms,ts,ls,mhs = sdf.spectrum.ModelSpectrum.read_kurucz(sdf.config.file['kurucz_models']+'fp00k2odfnew.pck')
    for i,m in enumerate(ms):
        assert(m.param_values['Teff']==ts[i])
        assert(m.param_values['[M/H]']==mhs[i])
        assert(m.param_values['logg']==ls[i])
