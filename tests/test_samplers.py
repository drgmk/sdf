import sys
from types import SimpleNamespace

import numpy as np
import pytest

from .context import sdf
import sdf.samplers


class FakeProblem:
    ndim = 2
    parameter_names = ('one', 'two')

    def __init__(self, offset):
        self.offset = offset
        self.seen = []

    def transform_prior(self, values):
        return np.asarray(values) + self.offset

    def log_likelihood(self, parameters):
        self.seen.append(parameters.copy())
        return float(np.sum(parameters))


def test_multinest_callbacks_capture_their_problem(monkeypatch):
    calls = []

    def run(log_likelihood, prior, ndim, **kwargs):
        cube = np.array([0.25, 0.75])
        prior(cube, ndim, ndim)
        calls.append((cube.copy(), log_likelihood(cube, ndim, ndim), kwargs))

    class Analyzer:
        def __init__(self, outputfiles_basename, n_params):
            self.outputfiles_basename = outputfiles_basename
            self.n_params = n_params

        def get_equal_weighted_posterior(self):
            return np.array([[1.0, 2.0, -3.0], [4.0, 5.0, -6.0]])

        def get_stats(self):
            return {
                'global evidence': 7.5,
                'global evidence error': 0.25,
            }

    monkeypatch.setitem(
        sys.modules, 'pymultinest',
        SimpleNamespace(run=run, Analyzer=Analyzer)
    )
    first = FakeProblem(1.0)
    second = FakeProblem(10.0)

    sampler = sdf.samplers.MultiNestSampler()
    first_run = sampler.run(first, 'first_')
    second_run = sampler.run(second, 'second_')

    assert calls[0][0].tolist() == [1.25, 1.75]
    assert calls[0][1] == 3.0
    assert calls[0][2]['outputfiles_basename'] == 'first_'
    assert calls[1][0].tolist() == [10.25, 10.75]
    assert calls[1][1] == 21.0
    assert first.seen[0].tolist() == [1.25, 1.75]
    assert second.seen[0].tolist() == [10.25, 10.75]
    assert first_run.parameter_names == ('one', 'two')
    assert first_run.samples.tolist() == [[1.0, 2.0], [4.0, 5.0]]
    assert first_run.log_likelihoods.tolist() == [-3.0, -6.0]
    assert first_run.evidence == 7.5
    assert first_run.evidence_error == 0.25
    assert second_run.native_output == 'second_'


def test_multinest_capabilities_support_evidence_and_resume():
    capabilities = sdf.samplers.MultiNestSampler.capabilities

    assert capabilities.evidence
    assert capabilities.weighted_samples
    assert capabilities.resumable


def test_multinest_rejects_fortran_output_names_over_limit():
    with pytest.raises(sdf.utils.SdfError, match='too long'):
        sdf.samplers.MultiNestSampler().run(
            FakeProblem(0.0), 'x' * 90
        )
