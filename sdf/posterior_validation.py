"""Scientific and structural validation for sampler output."""

from __future__ import annotations

import numpy as np

from . import utils


def validate_posterior(problem, posterior):
    """Reject posterior samples inconsistent with their fitting problem."""

    samples = np.asarray(posterior.samples, dtype=float)
    if not len(samples):
        raise utils.SdfError('posterior contains no samples')
    if samples.shape != (len(samples), problem.ndim):
        raise utils.SdfError(
            'posterior shape {} does not match {} fitting parameters'.format(
                samples.shape, problem.ndim
            )
        )
    if tuple(posterior.parameter_names) != problem.parameter_names:
        raise utils.SdfError(
            'posterior parameter names do not match the fitting problem'
        )
    if not np.all(np.isfinite(samples)):
        raise utils.SdfError('posterior contains non-finite parameter values')

    ranges = np.asarray(problem.parameter_ranges, dtype=float)
    if ranges.shape != (problem.ndim, 2):
        raise utils.SdfError('fitting parameter ranges are malformed')
    scale = np.maximum(1.0, np.max(np.abs(ranges), axis=1))
    tolerance = 1e-10 * scale
    too_low = samples < ranges[:, 0] - tolerance
    too_high = samples > ranges[:, 1] + tolerance
    outside = np.argwhere(too_low | too_high)
    if len(outside):
        _, parameter_index = outside[0]
        values = samples[:, parameter_index]
        raise utils.SdfError(
            'posterior parameter {} is outside [{:.8g}, {:.8g}] '
            '(observed {:.8g} to {:.8g})'.format(
                problem.parameter_names[parameter_index],
                ranges[parameter_index, 0], ranges[parameter_index, 1],
                np.min(values), np.max(values),
            )
        )

    log_likelihoods = posterior.log_likelihoods
    if log_likelihoods is not None:
        log_likelihoods = np.asarray(log_likelihoods, dtype=float)
        if not np.all(np.isfinite(log_likelihoods)):
            raise utils.SdfError(
                'posterior contains non-finite log-likelihoods'
            )
        check_indices = {0, len(samples) // 2, len(samples) - 1}
        check_indices.add(int(np.argmax(log_likelihoods)))
        for index in sorted(check_indices):
            recomputed = float(problem.log_likelihood(samples[index]))
            if not np.isfinite(recomputed) or not np.isclose(
                    log_likelihoods[index], recomputed,
                    rtol=1e-7, atol=1e-5):
                raise utils.SdfError(
                    'posterior log-likelihood is inconsistent with the '
                    'fitting problem at sample {} (stored {:.8g}, '
                    'recomputed {:.8g})'.format(
                        index, log_likelihoods[index], recomputed
                    )
                )

    for name, value in (
            ('evidence', posterior.evidence),
            ('evidence error', posterior.evidence_error)):
        if value is not None and not np.isfinite(value):
            raise utils.SdfError('{} is not finite'.format(name))
    if (posterior.evidence_error is not None and
            posterior.evidence_error < 0):
        raise utils.SdfError('evidence error is negative')
