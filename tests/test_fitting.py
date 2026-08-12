import numpy as np

from .context import sdf

def test_fitting_fill_colours():
    m = sdf.model.PhotModel.read_model('kurucz-0.0')
    fs = ['BS_YS','MIPS24']
    m.keep_filters(fs,colour_bases=True)
    mod_fnu = sdf.model.model_fluxes( (m,), [6300,4.5,-3], [2] )

def test_fitting_pmn_pc():
    n = 100
    prob = np.ones(n)
    samples = np.arange(n)+1
    pcs = np.random.randint(99,size=10)+1
    pc = sdf.fitting.pmn_pc(prob,samples,pcs)
    assert(np.all(pc==pcs))
    samples=np.tile(np.arange(n)+1,(3,1))
    pc = sdf.fitting.pmn_pc(prob,samples,pcs,axis=1)
    assert(np.all(pc.T==np.tile(pcs,(3,1))))
    samples=np.tile(np.arange(n)+1,(3,3,1))
    pc = sdf.fitting.pmn_pc(prob,samples,pcs,axis=2)
    assert(np.all(pc.T==np.tile(pcs,(3,3,1))))


def test_concat_obs_keeps_filters_from_multiple_photometry_objects():
    first = sdf.photometry.Photometry(
        filters=np.array(['2MJ']),
        fnujy=np.array([1.0]),
        e_fnujy=np.array([0.1]),
        upperlim=np.array([False]),
        ignore=np.array([False]),
        bibcode=np.array(['first']),
    )
    second = sdf.photometry.Photometry(
        filters=np.array(['2MH']),
        fnujy=np.array([2.0]),
        e_fnujy=np.array([0.2]),
        upperlim=np.array([False]),
        ignore=np.array([False]),
        bibcode=np.array(['second']),
    )

    values = sdf.fitting.concat_obs((first, second))

    assert values[7].tolist() == ['2MJ', '2MH']
    assert values[0].tolist() == [1.0, 2.0]
