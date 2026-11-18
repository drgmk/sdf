import numpy as np

from .context import sdf

def test_read_model_keep_filters():
    m = sdf.model.SpecModel.read_model('kurucz-0.0')
    m = sdf.model.PhotModel.read_model('kurucz-0.0')
    fs = ['BS_YS','MIPS24']
    m1 = m.copy()
    m1.keep_filters(fs,colour_bases=True)
    assert( np.array_equal(m1.filters,['BS_YS','MIPS24','BS','YS']) )
    m1.keep_filters(fs,colour_bases=False)
    assert( np.array_equal(m1.filters,['BS_YS','MIPS24']) )


def test_prepare_spectral_model_convolves_higher_resolution_grid():
    wavelength = np.geomspace(8.0, 12.0, 2000)
    flux = np.ones(len(wavelength))
    flux[np.argmin(np.abs(wavelength - 10.0))] = 100.0
    model = sdf.model.SpecModel(
        name='narrow-feature', parameters=[], param_values={},
        wavelength=wavelength, fnujy_sr=flux,
    )
    observed_wavelength = np.geomspace(8.5, 11.5, 80)

    model.prepare_for_spectrum(observed_wavelength, resolving_power=80.0)

    assert np.array_equal(model.wavelength, observed_wavelength)
    assert model.resolution_match['method'] == 'gaussian_convolution'
    assert not model.resolution_match['under_resolved']
    assert np.max(model.fnujy_sr) < 10.0
    assert np.min(model.fnujy_sr) > 0.99


def test_prepare_spectral_model_records_under_resolved_fallback():
    wavelength = np.geomspace(8.0, 12.0, 30)
    model = sdf.model.SpecModel(
        name='coarse', parameters=[], param_values={},
        wavelength=wavelength, fnujy_sr=np.ones(len(wavelength)),
    )
    observed_wavelength = np.geomspace(8.5, 11.5, 20)

    model.prepare_for_spectrum(observed_wavelength, resolving_power=600.0)

    assert model.resolution_match['method'] == 'interpolate'
    assert model.resolution_match['under_resolved']
    assert np.allclose(model.fnujy_sr, 1.0)
