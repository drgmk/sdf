"""Sampler-independent observations and target identities used by fitting."""

from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np

from . import photometry
from . import spectrum
from . import utils


@dataclass(frozen=True)
class TargetIdentity:
    """Identity and metadata for one physical model contributor."""

    key: str
    main_id: str | None
    sdbid: str | None
    keywords: dict

    @classmethod
    def from_keywords(cls, keywords, source=None):
        """Construct a stable identity, preferring the canonical SIMBAD ID."""

        main_id = _clean_keyword(keywords.get('main_id'))
        sdbid = _clean_keyword(keywords.get('id'))
        if main_id is not None:
            key = 'simbad:' + ' '.join(main_id.split())
        elif sdbid is not None:
            key = 'sdb:' + sdbid
        elif source is not None:
            key = 'local:' + os.path.abspath(source)
        else:
            raise utils.SdfError(
                'target has no SIMBAD main_id, SDB id, or source path'
            )

        return cls(key=key, main_id=main_id, sdbid=sdbid,
                   keywords=dict(keywords))


@dataclass(frozen=True)
class Observation:
    """Photometry/spectra observed for one physical or composite scope."""

    observation_id: str
    source: str
    target_key: str
    contributor_keys: tuple[str, ...]
    photometry: photometry.Photometry
    spectra: tuple[spectrum.ObsSpectrum, ...] = ()

    @property
    def data(self):
        """Return the legacy tuple consumed by SDF's model machinery."""

        return (self.photometry,) + self.spectra


@dataclass(frozen=True)
class FitDataset:
    """Validated observations and topology supplied to a fit."""

    targets: tuple[TargetIdentity, ...]
    observations: tuple[Observation, ...]

    @classmethod
    def from_rawphot(cls, rawphot, nospec=False):
        """Load the ordinary one-file/self-mapped form of a fit dataset."""

        rawphot = os.path.abspath(os.fspath(rawphot))
        observed_photometry = photometry.Photometry.read_sdb_file(rawphot)
        if observed_photometry is None:
            return None
        if np.all(observed_photometry.ignore):
            return None

        keywords = utils.get_sdb_keywords(rawphot)
        target = TargetIdentity.from_keywords(keywords, source=rawphot)
        observed_spectra = ()
        if not nospec:
            values = spectrum.ObsSpectrum.read_sdb_file(
                rawphot, module_split=True
            )
            if values is not None:
                observed_spectra = tuple(values)

        observation_id = target.sdbid or target.main_id or target.key
        observation = Observation(
            observation_id=observation_id,
            source=rawphot,
            target_key=target.key,
            contributor_keys=(target.key,),
            photometry=observed_photometry,
            spectra=observed_spectra,
        )
        return cls(targets=(target,), observations=(observation,))

    @property
    def data(self):
        """Flatten observations into the tuple consumed by legacy fitting."""

        return tuple(
            value for observation in self.observations
            for value in observation.data
        )

    @property
    def primary_keywords(self):
        """Return legacy per-result metadata for a single-target dataset."""

        if len(self.targets) != 1:
            raise utils.SdfError(
                'joint datasets do not have one primary keyword mapping'
            )
        return self.targets[0].keywords


def _clean_keyword(value):
    if value is None:
        return None
    value = str(value).strip()
    if value in ('', 'None'):
        return None
    return value
