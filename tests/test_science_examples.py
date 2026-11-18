import json
from pathlib import Path

import numpy as np
import pytest

from .context import sdf
import sdf.dataset
import sdf.likelihood
import sdf.model_selection


EXAMPLES = Path(__file__).resolve().parents[1] / 'examples' / 'fitting'


def _problem(directory):
    dataset = sdf.dataset.FitDataset.from_path(directory)
    candidate = sdf.model_selection.ModelCandidate.from_owned_models(
        (target.key, 'phoenix_sol') for target in dataset.targets
    )
    return sdf.likelihood.FitProblem.compile(dataset, candidate)


def _truth(directory, problem):
    payload = json.loads((directory / 'truth.json').read_text())
    by_name = payload['components']
    by_key = {
        target.key: by_name[target.main_id] for target in problem.dataset.targets
    }
    values = []
    for component in problem.candidate.components:
        values.extend(by_key[component.owner_key])
    values.extend(payload.get('spectrum_norms', ()))
    return np.asarray(values, dtype=float)


@pytest.mark.parametrize(
    'name', ('single-star', 'faint-binary', 'equal-binary', 'spectrum-star')
)
def test_synthetic_example_is_exact_at_recorded_truth(name):
    directory = EXAMPLES / name
    problem = _problem(directory)

    residuals, _, _ = problem.residual(_truth(directory, problem))

    assert np.max(np.abs(residuals)) < 1e-10
    assert problem.log_likelihood(_truth(directory, problem)) == pytest.approx(
        0.0, abs=1e-20
    )


def test_faint_binary_bright_component_reuses_single_star_photometry():
    single = sdf.dataset.FitDataset.from_path(EXAMPLES / 'single-star')
    binary = sdf.dataset.FitDataset.from_path(EXAMPLES / 'faint-binary')
    observation_a = next(
        observation for observation in binary.observations
        if observation.observation_id == 'A'
    )

    np.testing.assert_array_equal(
        single.observations[0].photometry.filters,
        observation_a.photometry.filters,
    )
    np.testing.assert_allclose(
        single.observations[0].photometry.fnujy,
        observation_a.photometry.fnujy,
        rtol=0.0, atol=0.0,
    )


def test_equal_binary_has_identical_components_and_summed_composite():
    dataset = sdf.dataset.FitDataset.from_path(EXAMPLES / 'equal-binary')
    observations = {
        observation.observation_id: observation.photometry
        for observation in dataset.observations
    }

    np.testing.assert_allclose(
        observations['A'].fnujy, observations['B'].fnujy,
        rtol=0.0, atol=0.0,
    )
    np.testing.assert_allclose(
        observations['AB'].fnujy,
        observations['A'].fnujy + observations['B'].fnujy,
        rtol=5e-15, atol=0.0,
    )


def test_spectrum_example_has_one_portable_segment_and_normalization():
    problem = _problem(EXAMPLES / 'spectrum-star')
    spectral = tuple(
        observation for observation in problem.dataset.observations
        if observation.kind == 'spectrum'
    )

    assert len(spectral) == 1
    assert spectral[0].product_id == 'synthetic:spectrum'
    assert spectral[0].segment_id == 'spectrum'
    assert len(spectral[0].spectra[0].wavelength) == 60
    assert problem.parameter_names[-1] == (
        'observation:synthetic:spectrum:spectrum.norm'
    )
