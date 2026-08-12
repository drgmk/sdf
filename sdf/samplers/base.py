"""Shared descriptions of sampler capabilities."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SamplerCapabilities:
    """Features supplied by a sampler implementation."""

    evidence: bool
    weighted_samples: bool
    resumable: bool
