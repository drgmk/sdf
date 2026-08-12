"""MultiNest adapter without process-global likelihood state."""

from __future__ import annotations

import numpy as np

from .base import SamplerCapabilities
from .. import config as cfg


class MultiNestSampler:
    """Run PyMultiNest against a sampler-independent fitting problem."""

    capabilities = SamplerCapabilities(
        evidence=True, weighted_samples=True, resumable=True
    )

    def run(self, problem, outputfiles_basename):
        """Run MultiNest, capturing problem state in local callbacks."""

        import pymultinest

        def prior(cube, ndim, nparams):
            values = problem.transform_prior(
                np.array([cube[i] for i in range(ndim)])
            )
            for i, value in enumerate(values):
                cube[i] = value

        def log_likelihood(cube, ndim, nparams):
            parameters = np.array([cube[i] for i in range(ndim)])
            return problem.log_likelihood(parameters)

        pymultinest.run(
            log_likelihood,
            prior,
            problem.ndim,
            n_live_points=cfg.fitting['n_live'],
            n_iter_before_update=cfg.fitting['n_update'],
            multimodal=True,
            sampling_efficiency=0.3,
            verbose=cfg.fitting['verb'],
            outputfiles_basename=outputfiles_basename,
        )
