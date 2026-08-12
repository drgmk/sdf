import numpy as np
import pytest
from astropy.table import Table

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


def test_result_id_removes_exact_rawphot_suffix(tmp_path):
    rawphot = tmp_path / 'target-rawphot.txt'
    rawphot.write_text('')

    result = sdf.result.BaseResult(str(rawphot), ('phoenix_sol',))

    assert result.id == 'target'


def test_result_loads_single_file_through_fit_dataset(tmp_path):
    available = sdf.model.PhotModel.read_model('kurucz-0.0').filters
    filter_name = next(name for name in available
                       if not sdf.filter.iscolour(name))
    rawphot = tmp_path / 'target-rawphot.txt'
    table = Table(
        rows=[(filter_name, 1.0, 0.1, 0.0, 0, 'Jy', 'test', '', 0)],
        names=('Band', 'Phot', 'Err', 'Sys', 'Lim', 'Unit', 'bibcode',
               'Note1', 'exclude'),
    )
    table.meta['keywords'] = {
        'id': {'value': 'sdb-test'},
        'main_id': {'value': 'HD 123'},
    }
    table.write(rawphot, format='ascii.ipac')
    result = sdf.result.BaseResult(str(rawphot), ('kurucz-0.0',))

    result.fill_data_models(nospec=True)

    assert result.dataset.targets[0].key == 'simbad:HD 123'
    assert result.model_candidate.model_names == ('kurucz-0.0',)
    assert result.fit_problem.observations == result.obs
    assert result.fit_problem.models == result.models
