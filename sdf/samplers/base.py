"""Shared descriptions of sampler capabilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SamplerCapabilities:
    """Features supplied by a sampler implementation."""

    evidence: bool
    weighted_samples: bool
    resumable: bool


@dataclass(frozen=True)
class PosteriorRun:
    """Compact sampler-neutral output for one completed candidate."""

    samples: np.ndarray
    parameter_names: tuple[str, ...]
    sampler: str
    log_likelihoods: np.ndarray | None = None
    weights: np.ndarray | None = None
    evidence: float | None = None
    evidence_error: float | None = None
    options: dict | None = None
    diagnostics: dict | None = None
    native_output: str | None = None

    def __post_init__(self):
        samples = np.asarray(self.samples, dtype=float)
        if samples.ndim != 2:
            raise ValueError('posterior samples must be a two-dimensional array')
        if samples.shape[1] != len(self.parameter_names):
            raise ValueError(
                'posterior sample width does not match parameter names'
            )
        object.__setattr__(self, 'samples', samples)
        object.__setattr__(
            self, 'parameter_names', tuple(self.parameter_names)
        )
        for name in ('log_likelihoods', 'weights'):
            values = getattr(self, name)
            if values is None:
                continue
            values = np.asarray(values, dtype=float)
            if values.shape != (len(samples),):
                raise ValueError(
                    '{} must contain one value per sample'.format(name)
                )
            object.__setattr__(self, name, values)
