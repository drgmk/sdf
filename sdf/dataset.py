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
from .spectral_manifest import read_spectra_manifest


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
    target_key: str | None
    contributor_keys: tuple[str, ...]
    photometry: photometry.Photometry | None = None
    spectra: tuple[spectrum.ObsSpectrum, ...] = ()
    kind: str = 'photometry'
    product_id: str | None = None
    segment_id: str | None = None

    @property
    def data(self):
        """Return the legacy tuple consumed by SDF's model machinery."""

        values = () if self.photometry is None else (self.photometry,)
        return values + self.spectra


@dataclass(frozen=True)
class FitDataset:
    """Validated observations and topology supplied to a fit."""

    targets: tuple[TargetIdentity, ...]
    observations: tuple[Observation, ...]
    manifests: tuple[str, ...] = ()

    @classmethod
    def from_rawphot(cls, rawphot, nospec=False):
        """Load the ordinary one-file/self-mapped form of a fit dataset."""

        loaded = _load_rawphot(rawphot, nospec=nospec)
        if loaded is None:
            return None
        target, rawphot, observed_photometry, observed_spectra = loaded
        observation_id = target.sdbid or target.main_id or target.key
        observations = _observations_for_loaded_rawphot(
            observation_id, target, rawphot, observed_photometry,
            observed_spectra, (target.key,),
        )
        return cls(targets=(target,), observations=observations)

    @classmethod
    def from_package(cls, directory, nospec=False):
        """Load rawphot files and topology from a joint-fit package."""

        directory = Path(directory).expanduser().resolve()
        topology_path = directory / 'joint-fit.yml'
        topology = (
            _read_joint_fit(topology_path) if topology_path.is_file() else {}
        )
        spectra_path = directory / 'spectra.yml'
        spectra_manifest = (
            None if nospec or not spectra_path.is_file() else
            read_spectra_manifest(spectra_path)
        )

        loaded_by_id = {}
        for rawphot in sorted(directory.glob('*-rawphot.txt')):
            loaded = _load_rawphot(
                rawphot,
                nospec=nospec or spectra_manifest is not None,
            )
            if loaded is None:
                keywords = utils.get_sdb_keywords(rawphot)
                target = TargetIdentity.from_keywords(
                    keywords, source=os.path.abspath(os.fspath(rawphot))
                )
                source = os.path.abspath(os.fspath(rawphot))
                observed_photometry = None
                observed_spectra = ()
            else:
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
            observations.extend(_observations_for_loaded_rawphot(
                observation_id, target, source, observed_photometry,
                observed_spectra, tuple(contributor_keys),
            ))

        if spectra_manifest is not None:
            observations.extend(_manifest_spectral_observations(
                spectra_manifest, loaded_by_id, physical_targets,
            ))

        if not observations:
            raise utils.SdfError(
                'fit package contains no usable photometry or spectra'
            )

        manifest_paths = [topology_path]
        if not nospec:
            manifest_paths.append(spectra_path)
        manifests = tuple(
            str(path) for path in manifest_paths if path.is_file()
        )
        return cls(
            targets=tuple(
                physical_targets[key] for key in sorted(physical_targets)
            ),
            observations=tuple(observations),
            manifests=manifests,
        )

    @classmethod
    def from_path(cls, source, nospec=False):
        """Load a rawphot file or its containing joint package."""

        source = canonical_fit_source(source)
        if source.is_dir():
            return cls.from_package(source, nospec=nospec)
        return cls.from_rawphot(source, nospec=nospec)

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


def canonical_fit_source(source):
    """Resolve any member of a joint package to the package directory."""

    source = Path(source).expanduser().resolve()
    if not source.exists():
        raise utils.SdfError('fit input does not exist: {}'.format(source))
    if source.is_dir():
        rawphot = sorted(source.glob('*-rawphot.txt'))
        marked = next((
            path for path in rawphot
            if 'sdf_package' in _rawphot_package_keywords(path)
        ), None)
        if marked is not None:
            package = _explicit_package_directory(marked)
            if package != source:
                raise utils.SdfError(
                    '{} points at a different fit package'.format(marked)
                )
            return source
        if (
            (source / 'joint-fit.yml').is_file() or
            (source / 'spectra.yml').is_file()
        ):
            return source
        if len(rawphot) == 1:
            return rawphot[0]
        raise utils.SdfError(
            '{} is not a joint package and contains {} rawphot files'.format(
                source, len(rawphot)
            )
        )
    if not source.name.endswith('-rawphot.txt'):
        raise utils.SdfError(
            'fit input is not a rawphot file: {}'.format(source)
        )
    package = _explicit_package_directory(source)
    if package is not None:
        return package
    if (
        (source.parent / 'joint-fit.yml').is_file() or
        (source.parent / 'spectra.yml').is_file()
    ):
        return source.parent
    return source


def _explicit_package_directory(rawphot):
    """Resolve and validate the package pointers exported in rawphot files."""

    keywords = _rawphot_package_keywords(rawphot)
    marker = keywords.get('sdf_package')
    if marker is None:
        if 'spectra_manifest' in keywords:
            raise utils.SdfError(
                '{} has spectra_manifest without sdf_package'.format(rawphot)
            )
        return None
    if marker != '.':
        raise utils.SdfError(
            '{} has unsupported sdf_package {!r}; expected \'.\''.format(
                rawphot, marker
            )
        )

    package = rawphot.parent.resolve()
    members = sorted(package.glob('*-rawphot.txt'))
    member_keywords = {
        path: _rawphot_package_keywords(path) for path in members
    }
    incomplete = [
        path.name for path, values in member_keywords.items()
        if values.get('sdf_package') != '.'
    ]
    if incomplete:
        raise utils.SdfError(
            '{} has rawphot files without the same sdf_package marker: {}'.
            format(package, ', '.join(incomplete))
        )
    spectra_references = {
        values.get('spectra_manifest')
        for values in member_keywords.values()
    }
    if len(spectra_references) != 1:
        raise utils.SdfError(
            '{} has inconsistent spectra_manifest pointers'.format(package)
        )
    spectra_reference = next(iter(spectra_references))
    if spectra_reference is not None:
        if spectra_reference != 'spectra.yml':
            raise utils.SdfError(
                '{} has unsupported spectra_manifest {!r}; expected '
                "'spectra.yml'".format(rawphot, spectra_reference)
            )
        manifest = package / spectra_reference
        if not manifest.is_file():
            raise utils.SdfError(
                '{} references missing spectra manifest {}'.format(
                    rawphot, manifest
                )
            )
    return package


def _rawphot_package_keywords(rawphot):
    """Read the two package pointers without parsing the full IPAC table."""

    values = {}
    with Path(rawphot).open('r', encoding='utf-8') as handle:
        for line in handle:
            if not line.startswith('\\'):
                break
            name, separator, value = line[1:].partition('=')
            name = name.strip()
            if not separator or name not in (
                'sdf_package', 'spectra_manifest'
            ):
                continue
            value = value.strip()
            if (
                len(value) >= 2 and value[0] == value[-1]
                and value[0] in ("'", '"')
            ):
                value = value[1:-1]
            values[name] = value
    return values


def discover_fit_sources(paths):
    """Expand paths and deduplicate joint packages across member files."""

    sources = {}
    for value in paths:
        path = Path(value).expanduser().resolve()
        if (
            path.is_dir() and
            not (path / 'joint-fit.yml').is_file() and
            not (path / 'spectra.yml').is_file()
        ):
            candidates = path.rglob('*-rawphot.txt')
        else:
            candidates = (path,)
        for candidate in candidates:
            source = canonical_fit_source(candidate)
            sources[str(source)] = source
    return tuple(sources[key] for key in sorted(sources))


def _load_rawphot(rawphot, nospec=False):
    rawphot = os.path.abspath(os.fspath(rawphot))
    observed_photometry = photometry.Photometry.read_sdb_file(rawphot)
    keywords = utils.get_sdb_keywords(rawphot)
    target = TargetIdentity.from_keywords(keywords, source=rawphot)
    observed_spectra = ()
    if not nospec:
        values = spectrum.ObsSpectrum.read_sdb_file(
            rawphot, module_split=True
        )
        if values is not None:
            observed_spectra = tuple(values)
    no_photometry = (
        observed_photometry is None or np.all(observed_photometry.ignore)
    )
    if no_photometry and not observed_spectra:
        return None
    return target, rawphot, observed_photometry, observed_spectra


def _observations_for_loaded_rawphot(
        observation_id, target, source, observed_photometry,
        observed_spectra, contributor_keys):
    """Split legacy rawphot contents into independently owned scopes."""

    observations = []
    if observed_photometry is not None:
        observations.append(Observation(
            observation_id=observation_id,
            source=source,
            target_key=target.key,
            contributor_keys=contributor_keys,
            photometry=observed_photometry,
        ))
    used_segments = set()
    for index, value in enumerate(observed_spectra, start=1):
        segment_id = _spectrum_segment_id(value, index, used_segments)
        product_id = '{}:legacy-spectrum'.format(observation_id)
        value.product_id = product_id
        value.segment_id = segment_id
        observations.append(Observation(
            observation_id='{}:{}'.format(product_id, segment_id),
            source=source,
            target_key=target.key,
            contributor_keys=contributor_keys,
            spectra=(value,),
            kind='spectrum',
            product_id=product_id,
            segment_id=segment_id,
        ))
    return tuple(observations)


def _manifest_spectral_observations(
        manifest, loaded_by_id, physical_targets):
    observations = []
    for definition in manifest.spectra:
        if not definition.include:
            continue
        missing = set(definition.contributors) - set(loaded_by_id)
        for override in definition.segment_overrides.values():
            if override.contributors is not None:
                missing.update(
                    set(override.contributors) - set(loaded_by_id)
                )
        if missing:
            raise utils.SdfError(
                'spectrum {} references missing contributors: {}'.format(
                    definition.spectrum_id, ', '.join(sorted(missing))
                )
            )
        values = spectrum.ObsSpectrum.read_file_of_type(
            str(definition.path), type_=definition.format,
            module_split=True,
        )
        if not values:
            raise utils.SdfError(
                'spectrum {} contains no usable segments'.format(
                    definition.spectrum_id
                )
            )
        used_segments = set()
        seen_overrides = set()
        for index, value in enumerate(values, start=1):
            segment_id = _spectrum_segment_id(
                value, index, used_segments
            )
            override = definition.segment_overrides.get(segment_id)
            if override is not None:
                seen_overrides.add(segment_id)
            include = (
                True if override is None or override.include is None else
                override.include
            )
            if not include:
                continue
            contributor_ids = (
                definition.contributors
                if override is None or override.contributors is None else
                override.contributors
            )
            contributor_keys = []
            for contributor_id in contributor_ids:
                target = loaded_by_id[contributor_id][0]
                physical_targets[target.key] = target
                contributor_keys.append(target.key)
            value.product_id = definition.spectrum_id
            value.segment_id = segment_id
            observations.append(Observation(
                observation_id='{}:{}'.format(
                    definition.spectrum_id, segment_id
                ),
                source=str(definition.path),
                target_key=None,
                contributor_keys=tuple(contributor_keys),
                spectra=(value,),
                kind='spectrum',
                product_id=definition.spectrum_id,
                segment_id=segment_id,
            ))
        unknown_overrides = (
            set(definition.segment_overrides) - seen_overrides
        )
        if unknown_overrides:
            raise utils.SdfError(
                'spectrum {} has overrides for missing segments: {}'.format(
                    definition.spectrum_id,
                    ', '.join(sorted(unknown_overrides)),
                )
            )
    return observations


def _spectrum_segment_id(value, index, used):
    candidate = getattr(value, 'segment_id', None)
    if candidate is None:
        candidate = getattr(value, 'segment', None)
    if candidate is None:
        instrument = getattr(value, 'instrument', None)
        if instrument:
            candidate = str(instrument).replace('Spitzer IRS ', '').strip()
    if not candidate:
        candidate = 'segment-{}'.format(index)
    candidate = str(candidate).strip().replace(' ', '-')
    base = candidate
    ordinal = 2
    while candidate in used:
        candidate = '{}-{}'.format(base, ordinal)
        ordinal += 1
    used.add(candidate)
    return candidate


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
