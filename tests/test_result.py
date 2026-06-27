import numpy as np
import pytest

from .context import sdf


def test_sort_identical_components_by_temperature():
    samples = np.array([
        [8000., -2., 2.0, 1.5, 0.5, 2.6, 2.6, 2.0, 1.5, 0.1],
        [8100., -2., 2.7, 2.1, 0.2, 2.1, 1.8, 0.8, 2.5, 1.0],
    ])
    model_comps = ('star', 'disk', 'disk')
    comp_parameters = (('Teff',), ('log_Temp', 'log_lam0', 'beta'),
                       ('log_Temp', 'log_lam0', 'beta'))

    sorted_samples = sdf.result.Result._sort_identical_components(
        samples, model_comps, comp_parameters
    )

    assert np.all(sorted_samples[:, 2] >= sorted_samples[:, 6])
    assert np.array_equal(sorted_samples[0, 2:6], samples[0, 6:10])
    assert np.array_equal(sorted_samples[0, 6:10], samples[0, 2:6])
    assert np.array_equal(sorted_samples[1], samples[1])
    assert np.array_equal(samples[0, 2:6], [2.0, 1.5, 0.5, 2.6])


def test_sort_does_not_mix_different_models():
    samples = np.array([[2.0, 1.0, 3.0, 4.0]])

    sorted_samples = sdf.result.Result._sort_identical_components(
        samples, ('disk_a', 'disk_b'), (('log_Temp',), ('log_Temp',))
    )

    assert np.array_equal(sorted_samples, samples)


def test_sort_rejects_different_parameters_for_same_model():
    with pytest.raises(sdf.utils.SdfError):
        sdf.result.Result._sort_identical_components(
            np.zeros((1, 6)), ('disk', 'disk'),
            (('log_Temp', 'beta'), ('log_Temp',))
        )
