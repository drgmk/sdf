import sys
from types import SimpleNamespace

import numpy as np

from .context import sdf
import sdf.samplers


class FakeProblem:
    ndim = 2

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

    monkeypatch.setitem(sys.modules, 'pymultinest', SimpleNamespace(run=run))
    first = FakeProblem(1.0)
    second = FakeProblem(10.0)

    sampler = sdf.samplers.MultiNestSampler()
    sampler.run(first, 'first_')
    sampler.run(second, 'second_')

    assert calls[0][0].tolist() == [1.25, 1.75]
    assert calls[0][1] == 3.0
    assert calls[0][2]['outputfiles_basename'] == 'first_'
    assert calls[1][0].tolist() == [10.25, 10.75]
    assert calls[1][1] == 21.0
    assert first.seen[0].tolist() == [1.25, 1.75]
    assert second.seen[0].tolist() == [10.25, 10.75]


def test_multinest_capabilities_support_evidence_and_resume():
    capabilities = sdf.samplers.MultiNestSampler.capabilities

    assert capabilities.evidence
    assert capabilities.weighted_samples
    assert capabilities.resumable
