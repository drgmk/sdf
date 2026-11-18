"""MultiNest adapter without process-global likelihood state."""

from __future__ import annotations

import os
from pathlib import Path
import numpy as np

from .base import PosteriorRun, SamplerCapabilities
from .. import config as cfg
from .. import utils


class MultiNestSampler:
    """Run PyMultiNest against a sampler-independent fitting problem."""

    capabilities = SamplerCapabilities(
        evidence=True, weighted_samples=True, resumable=True
    )

    def run(self, problem, outputfiles_basename):
        """Run MultiNest and return compact sampler-neutral output."""

        longest_output_name = outputfiles_basename + 'post_equal_weights.dat'
        if len(os.fsencode(longest_output_name)) >= 100:
            raise utils.SdfError(
                'MultiNest output basename is too long for its Fortran '
                'filenames; change into the output directory and use a short '
                'relative basename'
            )

        import pymultinest

        resume_file = Path(outputfiles_basename + 'resume.dat')
        had_resume_state = (
            resume_file.is_file() and resume_file.stat().st_size > 3
        )
        resume = _resume_output_is_usable(problem, outputfiles_basename)

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
            resume=resume,
        )
        posterior = self.read(problem, outputfiles_basename)
        if posterior.options is not None:
            posterior.options['resumed'] = bool(resume and had_resume_state)
        return posterior

    def read(self, problem, outputfiles_basename):
        """Read a completed native MultiNest run into a compact result."""

        import pymultinest

        analyzer = pymultinest.Analyzer(
            outputfiles_basename=outputfiles_basename,
            n_params=problem.ndim,
        )
        equal_weight = np.asarray(
            analyzer.get_equal_weighted_posterior(), dtype=float
        )
        stats = analyzer.get_stats()
        return PosteriorRun(
            samples=equal_weight[:, :-1],
            parameter_names=problem.parameter_names,
            sampler='multinest',
            log_likelihoods=equal_weight[:, -1],
            evidence=stats.get('global evidence'),
            evidence_error=stats.get('global evidence error'),
            options={
                'n_live_points': cfg.fitting['n_live'],
                'n_iter_before_update': cfg.fitting['n_update'],
                'multimodal': True,
                'sampling_efficiency': 0.3,
            },
            diagnostics={'stats': stats},
            native_output=outputfiles_basename,
        )


def _resume_output_is_usable(problem, outputfiles_basename):
    """Start clean when an extant equal-weight posterior is implausible."""

    path = Path(outputfiles_basename + 'post_equal_weights.dat')
    if not path.is_file():
        return True
    # More than two million typical SDF posterior rows is already far beyond
    # useful retained sampling and risks exhausting memory in np.loadtxt.
    if path.stat().st_size > 1024**3:
        return False
    try:
        with path.open('r', encoding='ascii') as handle:
            values = np.fromstring(handle.readline(), sep=' ')
    except (OSError, UnicodeError):
        return False
    if values.shape != (problem.ndim + 1,):
        return False
    sample = values[:-1]
    ranges = np.asarray(problem.parameter_ranges, dtype=float)
    if not np.all(np.isfinite(values)):
        return False
    if np.any(sample < ranges[:, 0]) or np.any(sample > ranges[:, 1]):
        return False
    recomputed = problem.log_likelihood(sample)
    return bool(
        np.isfinite(recomputed) and
        np.isclose(values[-1], recomputed, rtol=1e-7, atol=1e-5)
    )
