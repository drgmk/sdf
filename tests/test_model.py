import numpy as np

from .context import sdf

def test_read_model_keep_filters():
    m = sdf.model.SpecModel.read_model('kurucz-0.0')
    m = sdf.model.PhotModel.read_model('kurucz-0.0')
    fs = ['BS_YS','MIPS24']
    m1 = m.copy()
    m1.keep_filters(fs,colour_bases=True)
    assert( np.all( np.equal(m1.filters,
                             np.array(['BS_YS','MIPS24','BS','YS'])) ) )
    m.keep_filters(fs,colour_bases=False)
    assert( np.all( np.equal(m1.filters,
                             np.array(['BS_YS','MIPS24'])) ) )
