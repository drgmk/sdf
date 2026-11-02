import numpy as np
import pytest

from .context import sdf
import sdf.posterior_validation


class Problem:
    ndim = 2
    parameter_names = ('x', 'y')
    parameter_ranges = ((0.0, 1.0), (-2.0, 2.0))

    @staticmethod
    def log_likelihood(parameters):
        return -0.5 * float(np.sum(np.square(parameters)))


def _posterior(samples, log_likelihoods=None):
    samples = np.asarray(samples, dtype=float)
    if log_likelihoods is None:
        log_likelihoods = np.array([
            Problem.log_likelihood(sample) for sample in samples
        ])
    return sdf.samplers.PosteriorRun(
        samples=samples,
        parameter_names=Problem.parameter_names,
        sampler='test',
        log_likelihoods=log_likelihoods,
        evidence=-1.0,
        evidence_error=0.1,
    )


def test_valid_posterior_is_accepted():
    sdf.posterior_validation.validate_posterior(
        Problem(), _posterior([[0.25, -1.0], [0.75, 1.0]])
    )


def test_parameter_outside_model_grid_is_rejected():
    with pytest.raises(sdf.utils.SdfError, match='outside'):
        sdf.posterior_validation.validate_posterior(
            Problem(), _posterior([[0.25, -3.0]])
        )


def test_inconsistent_log_likelihood_is_rejected():
    with pytest.raises(sdf.utils.SdfError, match='inconsistent'):
        sdf.posterior_validation.validate_posterior(
            Problem(), _posterior([[0.25, 1.0]], np.array([5.0]))
        )
