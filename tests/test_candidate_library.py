import numpy as np
import pytest

from .context import sdf
from .test_dataset import _write_rawphot
import sdf.candidate_library
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
