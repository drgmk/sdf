import json

import numpy as np
import pytest

from .context import sdf
from .test_dataset import _write_rawphot
import sdf.analysis
import sdf.candidate_library
import sdf.diagnostic_plots
import sdf.likelihood
import sdf.model_selection
import sdf.samplers


def _completed_candidate(tmp_path, offset=0.0):
    rawphot = tmp_path / 'target-rawphot.txt'
    if not rawphot.exists():
        _write_rawphot(rawphot)
    dataset = sdf.dataset.FitDataset.from_rawphot(rawphot)
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('kurucz-0.0',)
    )
    problem = sdf.likelihood.FitProblem.compile(dataset, candidate)
    samples = np.arange(3 * problem.ndim, dtype=float).reshape(
        3, problem.ndim
    ) + offset
    posterior = sdf.samplers.PosteriorRun(
        samples=samples,
        parameter_names=problem.parameter_names,
        sampler='test',
        log_likelihoods=np.array([-3.0, -2.0, -1.0]),
        evidence=np.float64(12.5 + offset),
        evidence_error=np.float64(0.2),
        options={'live_points': 10},
    )
    return dataset, candidate, posterior


def test_candidate_library_round_trips_compact_result(tmp_path):
    dataset, candidate, posterior = _completed_candidate(tmp_path)
    library = sdf.candidate_library.CandidateLibrary(tmp_path)

    stored = library.commit(
        dataset, candidate, posterior,
        analysis={'chi_squared': 1.2},
        analysis_arrays={'prediction': np.array([1.0, 2.0])},
        provenance={'sdf_revision': 'abc123'},
        completed_at='2026-08-12T00:00:00+00:00',
    )

    assert stored.candidate_id == candidate.canonical_id
    assert stored.metadata['display_name'] == 'kurucz-0.0'
    assert stored.metadata['analysis']['chi_squared'] == 1.2
    assert stored.metadata['provenance']['sdf_revision'] == 'abc123'
    assert stored.metadata['posterior']['evidence'] == 12.5
    assert stored.analysis_arrays['prediction'].tolist() == [1.0, 2.0]
    assert np.array_equal(stored.posterior.samples, posterior.samples)
    assert np.array_equal(
        stored.posterior.log_likelihoods, posterior.log_likelihoods
    )
    assert stored.is_current(dataset)
    assert library.candidates()[0].candidate_id == candidate.canonical_id


def test_posterior_artifact_describes_fitted_derived_and_component_axes(
        tmp_path):
    dataset, candidate, posterior = _completed_candidate(tmp_path)
    component = candidate.components[0]
    derived = [{
        'name': component.component_id + '.radius_rsun',
        'label': 'stellar radius',
        'unit': 'R_sun',
        'owner_key': component.owner_key,
        'owner_keys': [component.owner_key],
        'component_id': component.component_id,
        'role': component.role,
    }]
    stored = sdf.candidate_library.CandidateLibrary(tmp_path).commit(
        dataset, candidate, posterior,
        analysis={
            'derived_parameters': derived,
            'targets': {
                component.owner_key: {'main_id': 'HD 123'},
            },
        },
        analysis_arrays={
            'derived_samples': np.arange(3, dtype=float)[:, None],
            'component_spectrum_median': np.ones((1, 2)),
            'spectrum_wavelength_micron': np.array([1.0, 2.0]),
        },
    )
    path = stored.path / stored.metadata['posterior']['file']

    with np.load(path) as artifact:
        assert artifact['artifact_type'].item() == 'sdf.posterior'
        assert artifact['schema_version'].item() == 1
        assert artifact['parameter_names'].tolist() == list(
            posterior.parameter_names
        )
        assert artifact['derived_parameter_names'].tolist() == [
            derived[0]['name']
        ]
        assert artifact['derived_parameter_labels'].tolist() == [
            'stellar radius'
        ]
        assert artifact['derived_parameter_units'].tolist() == ['R_sun']
        assert artifact['component_ids'].tolist() == [
            component.component_id
        ]
        assert artifact['component_roles'].tolist() == ['stellar']
        assert artifact['derived_samples'].shape == (3, 1)
        assert all(artifact[name].dtype.kind != 'O' for name in artifact.files)
        metadata = json.loads(artifact['metadata_json'].item())
        assert metadata['targets'] == [{
            'key': component.owner_key, 'main_id': 'HD 123',
        }]
        assert metadata['array_axes']['samples'] == [
            'sample', 'parameter_names'
        ]
        assert 'component_spectrum_median' not in artifact.files

    sed_path = stored.path / stored.metadata['sed_data']['file']
    with np.load(sed_path) as artifact:
        assert artifact['artifact_type'].item() == 'sdf.sed-data'
        assert artifact['schema_version'].item() == 1
        assert artifact['component_spectrum_median'].shape == (1, 2)
        assert artifact['spectrum_wavelength_micron'].tolist() == [1.0, 2.0]
        assert 'derived_samples' not in artifact.files
        metadata = json.loads(artifact['metadata_json'].item())
        assert metadata['array_axes']['component_spectrum_median'] == [
            'component_ids', 'spectrum_wavelength_micron'
        ]


def test_legacy_posterior_can_load_and_refresh_without_refitting(tmp_path):
    dataset, candidate, posterior = _completed_candidate(tmp_path)
    library = sdf.candidate_library.CandidateLibrary(tmp_path)
    stored = library.commit(
        dataset, candidate, posterior,
        analysis_arrays={'prediction': np.array([1.0, 2.0])},
    )
    path = stored.path / stored.metadata['posterior']['file']
    sed_path = stored.path / stored.metadata['sed_data']['file']
    with np.load(path) as current, np.load(sed_path) as sed_data:
        legacy = {
            'samples': current['samples'],
            'log_likelihoods': current['log_likelihoods'],
            'analysis__prediction': sed_data['prediction'],
        }
    np.savez_compressed(path, **legacy)
    sed_path.unlink()
    result_path = stored.path / 'result.json'
    result = json.loads(result_path.read_text(encoding='utf-8'))
    result['posterior'].pop('artifact_schema')
    result.pop('sed_data')
    result_path.write_text(json.dumps(result), encoding='utf-8')

    loaded = library.load(stored.path)
    assert loaded.analysis_arrays['prediction'].tolist() == [1.0, 2.0]
    assert library.refresh_artifacts() == 1
    assert library.refresh_artifacts() == 0
    with np.load(path) as refreshed:
        assert 'prediction' not in refreshed.files
        assert 'analysis__prediction' not in refreshed.files
        assert refreshed['parameter_names'].tolist() == list(
            posterior.parameter_names
        )
    with np.load(stored.path / 'sed-data.npz') as refreshed:
        assert refreshed['prediction'].tolist() == [1.0, 2.0]


def test_refit_replaces_exact_candidate_without_result_history(tmp_path):
    dataset, candidate, first = _completed_candidate(tmp_path)
    library = sdf.candidate_library.CandidateLibrary(tmp_path)
    library.commit(dataset, candidate, first)
    _, _, second = _completed_candidate(tmp_path, offset=100.0)

    stored = library.commit(dataset, candidate, second)

    assert len(library.candidates()) == 1
    assert np.array_equal(stored.posterior.samples, second.samples)
    assert len(list(stored.path.glob('posterior-*.npz'))) == 1


def test_selection_is_separate_from_retained_candidates(tmp_path):
    dataset, candidate, posterior = _completed_candidate(tmp_path)
    library = sdf.candidate_library.CandidateLibrary(tmp_path)
    stored = library.commit(dataset, candidate, posterior)
    selection = sdf.candidate_library.FitSelection(
        default=stored.candidate_id,
        displayed=(stored.candidate_id,),
        labels={stored.candidate_id: 'Preferred'},
        recommended=stored.candidate_id,
        default_source='operator',
    )

    library.write_selection(selection)

    assert library.read_selection() == selection
    with pytest.raises(ValueError, match='unknown candidates'):
        library.write_selection(sdf.candidate_library.FitSelection(
            default='missing', displayed=('missing',)
        ))


def test_input_change_marks_stored_candidate_non_current(tmp_path):
    dataset, candidate, posterior = _completed_candidate(tmp_path)
    library = sdf.candidate_library.CandidateLibrary(tmp_path)
    stored = library.commit(dataset, candidate, posterior)

    with (tmp_path / 'target-rawphot.txt').open('a', encoding='utf-8') as handle:
        handle.write('\n\\ changed input\n')

    assert not stored.is_current(dataset)


def test_posterior_validates_sample_shape():
    with pytest.raises(ValueError, match='width'):
        sdf.samplers.PosteriorRun(
            samples=np.ones((2, 3)),
            parameter_names=('a', 'b'),
            sampler='test',
        )


def test_candidate_diagnostic_plots_are_regenerable_files(tmp_path):
    dataset, candidate, posterior = _completed_candidate(tmp_path)
    rng = np.random.default_rng(4)
    posterior = sdf.samplers.PosteriorRun(
        samples=rng.normal(size=(48, len(posterior.parameter_names))),
        parameter_names=posterior.parameter_names,
        sampler='test',
    )
    metadata = [{
        'name': '{}.radius_rsun'.format(candidate.components[0].component_id),
        'owner_key': dataset.targets[0].key,
        'component_id': candidate.components[0].component_id,
        'role': 'stellar',
        'label': 'stellar radius',
        'unit': 'R_sun',
    }]
    stored = sdf.candidate_library.CandidateLibrary(tmp_path).commit(
        dataset, candidate, posterior,
        analysis={
            'schema': sdf.analysis.ANALYSIS_SCHEMA,
            'targets': {
                dataset.targets[0].key: {'main_id': 'HD 123'}
            },
            'derived_parameters': metadata,
        },
        analysis_arrays={
            'derived_samples': rng.lognormal(size=(48, 1)),
        },
    )

    products = sdf.diagnostic_plots.ensure_candidate_diagnostic_plots(stored)

    assert set(products) == {'fitted', 'derived'}
    assert all(path.stat().st_size > 1000 for path in products.values())
