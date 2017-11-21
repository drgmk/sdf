import pytest
import numpy as np

from .context import sdf

def test_filter_common_x():
    x1 = np.array([1,2,3])
    y1 = np.array([1,2,3])
    x2 = np.array([1.5,2.5,3.5])
    y2 = np.array([4,5,6])
    x,newy1,newy2 = sdf.filter.common_x(x1,y1,x2,y2)
    assert(np.all(np.equal(x,np.array([1.5,2,2.5,3]))))
    assert(np.all(np.equal(newy1,np.array([1.5,2.,2.5,3]))))
    assert(np.all(np.equal(newy2,np.array([4.,4.5,5,5.5]))))

def test_filter_arrays_are_floats():
    f = sdf.filter.Filter(nu_hz=[1,2,3],response=[1,2,3])
    assert(f.nu_hz.dtype==float)
    assert(f.response.dtype==float)

def test_filter_get_filter():
    fs = sdf.filter.Filter.all
    for f in fs:
        if not sdf.filter.iscolour(f):
            filt = sdf.filter.Filter.get(f)
            assert( isinstance(filt,sdf.filter.Filter) )

def test_filter_zero_point():
    f = sdf.filter.Filter(zero_point=1234.0,zero_point_offset=0.0)
    assert(f.mag2flux(0.)==1234.0)
    assert(f.mag2flux(2.5)==1234.0/10.0)
    assert(f.mag2flux(-2.5)==1234.0*10.0)

def test_filter_flux_conversion():
    f = sdf.filter.Filter.get('VJ')
    assert(np.abs(f.flux2mag(f.mag2flux(0.0)))<1e-16)

def test_filter_sort():
    f = sdf.filter.Filter(nu_hz=[3,2,1],response=[1,2,3])
    f.sort()
    assert(np.all(np.equal(f.nu_hz,np.array([1.,2.,3]))))
    assert(np.all(np.equal(f.response,np.array([3.,2.,1.]))))

def test_filter_normalisation():
    f = sdf.filter.Filter(nu_hz=[1,3,5],response=[0,1,0])
    f.normalise_response()
    assert(np.all(np.equal(f.response,np.array([0.,0.5,0.]))))

def test_filter_actual_flux():
    nu = [1,100]
    fl = [1,10000]
    s = sdf.spectrum.ObsSpectrum(nu_hz=nu,fnujy=fl,e_fnujy=fl)
    f = sdf.filter.Filter(ref_nu_hz=10.0)
    assert(f.actual_flux(s)==100)

def test_filter_fill_ref_nu_hz():
    f = sdf.filter.Filter(ref_wavelength=1.0)
    f.fill_ref_nu_hz()
    assert(f.ref_nu_hz==sdf.filter.c_micron/1.0)

def test_filter_synthphot_no_ref():
    wv = np.arange(1,3,0.1)
    fl = np.ones(len(wv))
    s = sdf.spectrum.ModelSpectrum(wavelength=wv,fnujy_sr=fl)
    s.fill_wave2hz()
    f = sdf.filter.Filter.get('2MH')
    flux,cc = f.synthphot(s)
    assert(np.abs(flux-1.0)<1e15)
    assert(cc is None)

def test_filter_synthphot_ref():
    wv = np.arange(1,3,0.1)
    fl = np.ones(len(wv))
    ref = lambda x: 1
    s = sdf.spectrum.ModelSpectrum(wavelength=wv,fnujy_sr=fl)
    s.fill_wave2hz()
    f = sdf.filter.Filter.get('2MH')
    f.ref_spectrum = ref
    f.ref_nu_hz = sdf.filter.c_micron/f.mean_wavelength
    f.fill_cc_denom()
    flux,cc = f.synthphot(s)
    assert(np.abs(flux-1.0)<1e15)
    assert(np.abs(cc-1.0)<1e15)

def test_filter_cc_irac():
    filt = sdf.filter.Filter.get('IRAC3P6')
    temps = [5e3,2e3,1500,1e3,800,600,400,200]
    # ccs from Table 4.4 at
    # https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/18/
    ccs = [1.0063,0.9990,0.9959,0.9933,0.9953,1.0068,1.0614,1.5138]
    for i,t in enumerate(temps):
        b = sdf.spectrum.ModelSpectrum.bnu_nu_hz(filt.nu_hz,t)
        fnu,cc = filt.synthphot(b)
        assert(np.abs(cc/ccs[i]-1)<2.7e-2)

def test_filter_cc_mips24():
    filt = sdf.filter.Filter.get('MIPS24')
    temps = [1e4,5e3,1e3,500,300,200,150,100,70,50,30,20]
    # ccs from Table 4.16 at
    # http://irsa.ipac.caltech.edu/data/SPITZER/docs/mips/mipsinstrumenthandbook/51/
    ccs = [1.000,0.999,0.992,0.983,0.970,0.957,
           0.948,0.947,0.986,1.119,2.031,7.005]
    for i,t in enumerate(temps):
        b = sdf.spectrum.ModelSpectrum.bnu_nu_hz(filt.nu_hz,t)
        fnu,cc = filt.synthphot(b)
        assert(np.abs(cc/ccs[i]-1)<3e-3)

def test_filter_cc_mips70():
    filt = sdf.filter.Filter.get('MIPS70')
    temps = [1e4,5e3,1e3,500,300,200,150,100,70,50,30,20]
    # ccs from Table 4.16 at
    # http://irsa.ipac.caltech.edu/data/SPITZER/docs/mips/mipsinstrumenthandbook/51/
    ccs = [1.000,0.999,0.995,0.989,0.980,0.970,
           0.959,0.938,0.914,0.893,0.901,1.052]
    for i,t in enumerate(temps):
        b = sdf.spectrum.ModelSpectrum.bnu_nu_hz(filt.nu_hz,t)
        fnu,cc = filt.synthphot(b)
        assert(np.abs(cc/ccs[i]-1)<7e-4)

def test_filter_mean_wavelength():
    wv = np.array([1,2,3,4,5])
    rs = np.array([1,2,3,4,5])
    f = sdf.filter.Filter(nu_hz=sdf.filter.c_micron/wv,response=rs)
    f.fill_mean_wavelength()
    assert(f.mean_wavelength==12.0/4.0)

def test_iscolour():
    fs = ('UJ','BS_YS','MIPS24','STROMM1')
    result = sdf.filter.iscolour(fs)
    assert(np.all(np.equal(result,[False,True,False,True])))

def test_colour_get():
    fs = sdf.filter.Filter.all
    ntry = 10
    for f1 in [fs[i] for i in np.random.randint(len(fs),size=ntry)]:
        if sdf.filter.iscolour(f1):
            continue
        for f2 in [fs[i] for i in np.random.randint(len(fs),size=ntry)]:
            if sdf.filter.iscolour(f2):
                continue
            name = f1+'_'+f2
            if f1 == f2:
                with pytest.raises(sdf.utils.SdfError):
                    c = sdf.filter.Colour.get(name)
            else:
                c = sdf.filter.Colour.get(name)
                assert(len(c.filters)==2)
                assert(len(c.weights)==2)

    for f1 in fs:
        if sdf.filter.iscolour(f1):
            sdf.filter.Colour.get(f1)
