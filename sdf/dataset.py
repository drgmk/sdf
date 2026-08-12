"""Sampler-independent observations and target identities used by fitting."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

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

        loaded = _load_rawphot(rawphot, nospec=nospec)
        if loaded is None:
            return None
        target, rawphot, observed_photometry, observed_spectra = loaded
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

    @classmethod
    def from_package(cls, directory, nospec=False):
        """Load rawphot files and topology from a joint-fit package."""

        directory = Path(directory).expanduser().resolve()
        topology_path = directory / 'joint-fit.yml'
        if not topology_path.is_file():
            raise utils.SdfError(
                'joint-fit.yml not found in {}'.format(directory)
            )
        topology = _read_joint_fit(topology_path)

        loaded_by_id = {}
        for rawphot in sorted(directory.glob('*-rawphot.txt')):
            loaded = _load_rawphot(rawphot, nospec=nospec)
            if loaded is None:
                raise utils.SdfError(
                    'joint observation {} has no usable photometry'.
                    format(rawphot)
                )
            target, source, observed_photometry, observed_spectra = loaded
            observation_id = target.sdbid or target.main_id or target.key
            if observation_id in loaded_by_id:
                raise utils.SdfError(
                    'duplicate joint observation id {}'.format(observation_id)
                )
            loaded_by_id[observation_id] = (
                target, source, observed_photometry, observed_spectra
            )

        if not loaded_by_id:
            raise utils.SdfError(
                'no *-rawphot.txt files found in {}'.format(directory)
            )

        referenced = set(topology)
        for contributors in topology.values():
            referenced.update(contributors)
        missing = referenced - set(loaded_by_id)
        if missing:
            raise utils.SdfError(
                'joint-fit.yml references missing observations: {}'.
                format(', '.join(sorted(missing)))
            )

        observations = []
        physical_targets = {}
        for observation_id in sorted(loaded_by_id):
            target, source, observed_photometry, observed_spectra = \
                loaded_by_id[observation_id]
            contributor_ids = topology.get(observation_id, (observation_id,))
            contributor_keys = []
            for contributor_id in contributor_ids:
                contributor = loaded_by_id[contributor_id][0]
                physical_targets[contributor.key] = contributor
                contributor_keys.append(contributor.key)
            observations.append(Observation(
                observation_id=observation_id,
                source=source,
                target_key=target.key,
                contributor_keys=tuple(contributor_keys),
                photometry=observed_photometry,
                spectra=observed_spectra,
            ))

        return cls(
            targets=tuple(
                physical_targets[key] for key in sorted(physical_targets)
            ),
            observations=tuple(observations),
        )

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

    def component_projection(self, candidate):
        """Map each flattened observation object to visible model owners."""

        target_keys = {target.key for target in self.targets}
        unknown = {
            component.owner_key for component in candidate.components
        } - target_keys
        if unknown:
            raise utils.SdfError(
                'candidate has unknown target owners: {}'.
                format(', '.join(sorted(unknown)))
            )

        rows = []
        for observation in self.observations:
            row = [
                component.owner_key in observation.contributor_keys
                for component in candidate.components
            ]
            if not any(row):
                raise utils.SdfError(
                    'candidate has no model component visible to {}'.
                    format(observation.observation_id)
                )
            rows.extend([row] * len(observation.data))
        return np.asarray(rows, dtype=bool)


def _clean_keyword(value):
    if value is None:
        return None
    value = str(value).strip()
    if value in ('', 'None'):
        return None
    return value


def _load_rawphot(rawphot, nospec=False):
    rawphot = os.path.abspath(os.fspath(rawphot))
    observed_photometry = photometry.Photometry.read_sdb_file(rawphot)
    if observed_photometry is None or np.all(observed_photometry.ignore):
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
    return target, rawphot, observed_photometry, observed_spectra


def _read_joint_fit(path):
    try:
        payload = yaml.safe_load(path.read_text(encoding='utf-8'))
    except (OSError, yaml.YAMLError) as error:
        raise utils.SdfError(
            'could not read {}: {}'.format(path, error)
        ) from error
    if not isinstance(payload, dict):
        raise utils.SdfError('joint-fit.yml must contain a mapping')
    if set(payload) - {'version', 'observations'}:
        raise utils.SdfError('joint-fit.yml contains unknown keys')
    if payload.get('version') != 1:
        raise utils.SdfError(
            'unsupported joint-fit version {}'.format(payload.get('version'))
        )
    mappings = payload.get('observations')
    if not isinstance(mappings, dict):
        raise utils.SdfError('joint-fit observations must be a mapping')

    normalized = {}
    for observation_id, contributors in mappings.items():
        if not isinstance(observation_id, str) or not observation_id.strip():
            raise utils.SdfError('joint observation ids must be strings')
        if (not isinstance(contributors, list) or not contributors or
                any(not isinstance(value, str) or not value.strip()
                    for value in contributors)):
            raise utils.SdfError(
                'contributors for {} must be a non-empty list of strings'.
                format(observation_id)
            )
        if len(set(contributors)) != len(contributors):
            raise utils.SdfError(
                'observation {} repeats a contributor'.format(observation_id)
            )
        normalized[observation_id] = tuple(contributors)
    return normalized
