"""Compact reusable candidate artifacts and display selection."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile

import numpy as np

from .model_selection import ModelCandidate
from .samplers import PosteriorRun


POSTERIOR_ARTIFACT_SCHEMA = 1


@dataclass(frozen=True)
class StoredCandidate:
    """One candidate artifact loaded from a package library."""

    path: Path
    metadata: dict
    posterior: PosteriorRun
    analysis_arrays: dict[str, np.ndarray]

    @property
    def candidate_id(self):
        return self.metadata['candidate_id']

    def is_current(self, dataset, fit_state=None):
        if self.metadata['input_fingerprint'] != input_fingerprint(dataset):
            return False
        if fit_state is None:
            return True
        return self.metadata.get('fit_fingerprint') == fit_fingerprint(
            fit_state
        )


@dataclass(frozen=True)
class FitSelection:
    """Default and displayed candidates for one package."""

    default: str
    displayed: tuple[str, ...]
    labels: dict[str, str] | None = None
    recommended: str | None = None
    default_source: str = 'automatic'

    def __post_init__(self):
        displayed = tuple(self.displayed)
        if not displayed:
            raise ValueError('fit selection must display at least one candidate')
        if self.default not in displayed:
            raise ValueError('default candidate must also be displayed')
        if len(set(displayed)) != len(displayed):
            raise ValueError('displayed candidates must be unique')
        recommended = (
            self.default if self.recommended is None else self.recommended
        )
        if self.default_source not in {'automatic', 'operator'}:
            raise ValueError('default source must be automatic or operator')
        object.__setattr__(self, 'displayed', displayed)
        object.__setattr__(self, 'recommended', recommended)


class CandidateLibrary:
    """Filesystem store for completed candidate fits in one package."""

    def __init__(self, package_directory):
        self.package_directory = Path(package_directory).expanduser().resolve()
        self.fits_directory = self.package_directory / 'fits'

    def commit(self, dataset, candidate, posterior, *, analysis=None,
               analysis_arrays=None, policy=None, provenance=None,
               fit_state=None,
               completed_at=None):
        """Atomically expose a completed candidate and compact posterior."""

        if not isinstance(candidate, ModelCandidate):
            raise TypeError('candidate must be a ModelCandidate')
        if not isinstance(posterior, PosteriorRun):
            raise TypeError('posterior must be a PosteriorRun')

        candidate_directory = self.fits_directory / candidate.storage_name
        candidate_directory.mkdir(parents=True, exist_ok=True)
        posterior_digest = _posterior_fingerprint(posterior)
        posterior_name = 'posterior-{}.npz'.format(posterior_digest[:12])
        posterior_path = candidate_directory / posterior_name

        metadata = {
            'candidate_id': candidate.canonical_id,
            'display_name': candidate.candidate_id,
            'input_fingerprint': input_fingerprint(dataset),
            'fit_fingerprint': fit_fingerprint(fit_state or {}),
            'completed_at': completed_at or datetime.now(
                timezone.utc
            ).isoformat(),
            'components': [
                {
                    'component_id': component.component_id,
                    'owner_key': component.owner_key,
                    'role': component.role,
                    'model_name': component.model_name,
                    'ordinal': component.ordinal,
                }
                for component in candidate.components
            ],
            'posterior': {
                'file': posterior_name,
                'n_samples': len(posterior.samples),
                'parameter_names': list(posterior.parameter_names),
                'has_log_likelihoods': posterior.log_likelihoods is not None,
                'has_weights': posterior.weights is not None,
                'evidence': (
                    None if posterior.evidence is None
                    else float(posterior.evidence)
                ),
                'evidence_error': (
                    None if posterior.evidence_error is None
                    else float(posterior.evidence_error)
                ),
                'sampler': posterior.sampler,
                'options': posterior.options or {},
                'diagnostics': posterior.diagnostics or {},
                'native_output': posterior.native_output,
                'artifact_schema': POSTERIOR_ARTIFACT_SCHEMA,
            },
            'analysis': analysis or {},
            'analysis_arrays': sorted((analysis_arrays or {}).keys()),
            'policy': policy or {},
            'provenance': provenance or {},
        }
        arrays = _posterior_artifact_arrays(
            posterior, analysis_arrays or {}, metadata
        )
        _atomic_npz(posterior_path, arrays)
        _atomic_json(candidate_directory / 'result.json', metadata)

        for old_posterior in candidate_directory.glob('posterior-*.npz'):
            if old_posterior != posterior_path:
                old_posterior.unlink()
        return self.load(candidate_directory)

    def load(self, candidate_directory):
        """Load one completed candidate artifact."""

        candidate_directory = Path(candidate_directory)
        metadata = json.loads(
            (candidate_directory / 'result.json').read_text(encoding='utf-8')
        )
        posterior_meta = metadata['posterior']
        with np.load(candidate_directory / posterior_meta['file']) as arrays:
            posterior = PosteriorRun(
                samples=arrays['samples'],
                parameter_names=tuple(posterior_meta['parameter_names']),
                sampler=posterior_meta['sampler'],
                log_likelihoods=(
                    arrays['log_likelihoods']
                    if posterior_meta['has_log_likelihoods'] else None
                ),
                weights=(
                    arrays['weights']
                    if posterior_meta['has_weights'] else None
                ),
                evidence=posterior_meta['evidence'],
                evidence_error=posterior_meta['evidence_error'],
                options=posterior_meta['options'],
                diagnostics=posterior_meta['diagnostics'],
                native_output=posterior_meta.get('native_output'),
            )
            analysis_arrays = {}
            for name in metadata.get('analysis_arrays', ()):
                key = name if name in arrays else 'analysis__' + name
                analysis_arrays[name] = arrays[key]
        return StoredCandidate(
            candidate_directory, metadata, posterior, analysis_arrays
        )

    def candidates(self):
        """Return every complete artifact, ignoring interrupted directories."""

        if not self.fits_directory.is_dir():
            return ()
        return tuple(
            self.load(path.parent)
            for path in sorted(self.fits_directory.glob('*/result.json'))
        )

    def find(self, candidate_id):
        """Return a retained candidate by canonical ID, if present."""

        for candidate in self.candidates():
            if candidate.candidate_id == candidate_id:
                return candidate
        return None

    def refresh_posterior_artifacts(self):
        """Make retained posterior files self-describing without refitting."""

        refreshed = 0
        for candidate in self.candidates():
            posterior_path = (
                candidate.path / candidate.metadata['posterior']['file']
            )
            if _posterior_artifact_is_current(
                    posterior_path, candidate.metadata):
                continue
            metadata = dict(candidate.metadata)
            metadata['posterior'] = dict(metadata['posterior'])
            metadata['posterior'][
                'artifact_schema'
            ] = POSTERIOR_ARTIFACT_SCHEMA
            arrays = _posterior_artifact_arrays(
                candidate.posterior, candidate.analysis_arrays, metadata
            )
            _atomic_npz(posterior_path, arrays)
            _atomic_json(candidate.path / 'result.json', metadata)
            refreshed += 1
        return refreshed

    def write_selection(self, selection):
        """Atomically update the cheap mutable display/default choice."""

        available = {candidate.candidate_id for candidate in self.candidates()}
        referenced = set(selection.displayed) | {
            selection.default, selection.recommended,
        }
        missing = referenced - available
        if missing:
            raise ValueError(
                'selection references unknown candidates: {}'.format(
                    ', '.join(sorted(missing))
                )
            )
        self.fits_directory.mkdir(parents=True, exist_ok=True)
        _atomic_json(self.fits_directory / 'selection.json', {
            'default': selection.default,
            'recommended': selection.recommended,
            'default_source': selection.default_source,
            'displayed': list(selection.displayed),
            'labels': selection.labels or {},
        })

    def read_selection(self):
        payload = json.loads(
            (self.fits_directory / 'selection.json').read_text(
                encoding='utf-8'
            )
        )
        return FitSelection(
            default=payload['default'],
            displayed=tuple(payload['displayed']),
            labels=payload.get('labels') or {},
            recommended=payload.get('recommended', payload['default']),
            default_source=payload.get('default_source', 'automatic'),
        )


def input_fingerprint(dataset):
    """Hash portable package inputs rather than their absolute locations."""

    paths = {Path(observation.source) for observation in dataset.observations}
    paths.update(Path(value) for value in dataset.manifests)
    if dataset.manifests:
        root = Path(dataset.manifests[0]).parent
    elif len(paths) == 1:
        root = next(iter(paths)).parent
    else:
        root = Path(os.path.commonpath([str(path) for path in paths]))
        if root.is_file():
            root = root.parent

    def portable_name(path):
        try:
            return path.relative_to(root).as_posix()
        except ValueError:
            return path.name

    digest = sha256()
    for path in sorted(paths, key=portable_name):
        digest.update(portable_name(path).encode('utf-8'))
        digest.update(b'\0')
        digest.update(path.read_bytes())
        digest.update(b'\0')
    return digest.hexdigest()


def fit_fingerprint(fit_state):
    """Hash material numerical settings independently of display policy."""

    encoded = json.dumps(
        fit_state, sort_keys=True, separators=(',', ':'),
        default=_json_default
    ).encode('utf-8')
    return sha256(encoded).hexdigest()


def _posterior_fingerprint(posterior):
    digest = sha256()
    digest.update(np.ascontiguousarray(posterior.samples).view(np.uint8))
    for values in (posterior.log_likelihoods, posterior.weights):
        if values is not None:
            digest.update(np.ascontiguousarray(values).view(np.uint8))
    digest.update('\0'.join(posterior.parameter_names).encode('utf-8'))
    return digest.hexdigest()


def _posterior_artifact_arrays(posterior, analysis_arrays, metadata):
    """Return safe, self-describing arrays for a retained posterior."""

    arrays = {
        'artifact_type': np.asarray('sdf.posterior'),
        'schema_version': np.asarray(
            POSTERIOR_ARTIFACT_SCHEMA, dtype=np.int64
        ),
        'samples': posterior.samples,
        'parameter_names': _string_array(posterior.parameter_names),
    }
    if posterior.log_likelihoods is not None:
        arrays['log_likelihoods'] = posterior.log_likelihoods
    if posterior.weights is not None:
        arrays['weights'] = posterior.weights
    reserved = set(arrays) | {
        'metadata_json', 'derived_parameter_names',
        'derived_parameter_labels', 'derived_parameter_units',
        'derived_parameter_owner_keys', 'derived_parameter_component_ids',
        'derived_parameter_roles', 'component_ids', 'component_owner_keys',
        'component_roles', 'component_model_names',
    }
    for name, values in analysis_arrays.items():
        if name in reserved:
            raise ValueError(
                'analysis array name is reserved: {}'.format(name)
            )
        values = np.asarray(values)
        if values.dtype.hasobject:
            try:
                values = values.astype(str)
            except (TypeError, ValueError) as error:
                raise TypeError(
                    'analysis array {} cannot use object dtype'.format(name)
                ) from error
        arrays[name] = values

    derived = metadata.get('analysis', {}).get('derived_parameters', ())
    derived_samples = arrays.get('derived_samples')
    if derived_samples is not None and (
            derived_samples.ndim != 2 or
            derived_samples.shape[1] != len(derived)):
        raise ValueError(
            'derived_samples columns do not match derived parameter metadata'
        )
    arrays.update({
        'derived_parameter_names': _metadata_strings(derived, 'name'),
        'derived_parameter_labels': _metadata_strings(derived, 'label'),
        'derived_parameter_units': _metadata_strings(derived, 'unit'),
        'derived_parameter_owner_keys': _metadata_strings(
            derived, 'owner_key'
        ),
        'derived_parameter_component_ids': _metadata_strings(
            derived, 'component_id'
        ),
        'derived_parameter_roles': _metadata_strings(derived, 'role'),
    })

    components = metadata.get('components', ())
    arrays.update({
        'component_ids': _metadata_strings(components, 'component_id'),
        'component_owner_keys': _metadata_strings(components, 'owner_key'),
        'component_roles': _metadata_strings(components, 'role'),
        'component_model_names': _metadata_strings(components, 'model_name'),
    })
    arrays['metadata_json'] = np.asarray(json.dumps(
        _posterior_artifact_metadata(metadata),
        sort_keys=True, separators=(',', ':'), default=_json_default,
    ))
    return arrays


def _posterior_artifact_metadata(metadata):
    analysis = metadata.get('analysis', {})
    targets = []
    for key, values in analysis.get('targets', {}).items():
        targets.append({
            'key': key,
            'main_id': values.get('main_id'),
        })
    return {
        'artifact_type': 'sdf.posterior',
        'schema_version': POSTERIOR_ARTIFACT_SCHEMA,
        'description': (
            'Retained fitted posterior samples and regenerable SDF analysis '
            'arrays. Rows of samples and derived_samples correspond.'
        ),
        'candidate_id': metadata['candidate_id'],
        'display_name': metadata['display_name'],
        'parameter_names': metadata['posterior']['parameter_names'],
        'derived_parameters': analysis.get('derived_parameters', []),
        'components': metadata.get('components', []),
        'targets': targets,
        'array_axes': {
            'samples': ['sample', 'parameter_names'],
            'derived_samples': ['sample', 'derived_parameter_names'],
            'component_prediction_p16': ['component_ids', 'observation'],
            'component_prediction_median': [
                'component_ids', 'observation'
            ],
            'component_prediction_p84': ['component_ids', 'observation'],
            'component_spectrum_median': [
                'component_ids', 'spectrum_wavelength_micron'
            ],
        },
    }


def _metadata_strings(records, key):
    return _string_array(
        '' if record.get(key) is None else record.get(key, '')
        for record in records
    )


def _string_array(values):
    return np.asarray(tuple(str(value) for value in values), dtype=str)


def _posterior_artifact_is_current(path, metadata):
    if metadata['posterior'].get(
            'artifact_schema') != POSTERIOR_ARTIFACT_SCHEMA:
        return False
    required = {
        'artifact_type', 'schema_version', 'samples', 'parameter_names',
        'metadata_json', 'derived_parameter_names',
        'derived_parameter_labels', 'derived_parameter_units',
        'derived_parameter_owner_keys', 'derived_parameter_component_ids',
        'derived_parameter_roles', 'component_ids', 'component_owner_keys',
        'component_roles', 'component_model_names',
    }
    with np.load(path) as arrays:
        return (
            required.issubset(arrays.files) and
            not any(name.startswith('analysis__') for name in arrays.files) and
            str(arrays['artifact_type']) == 'sdf.posterior' and
            int(arrays['schema_version']) == POSTERIOR_ARTIFACT_SCHEMA
        )


def _atomic_npz(path, arrays):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        prefix='.posterior-', suffix='.npz', dir=path.parent, delete=False
    )
    temporary = Path(handle.name)
    handle.close()
    try:
        np.savez_compressed(temporary, **arrays)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _atomic_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode='w', encoding='utf-8', prefix='.result-', suffix='.json',
        dir=path.parent, delete=False
    )
    temporary = Path(handle.name)
    try:
        with handle:
            json.dump(
                payload, handle, indent=2, sort_keys=True,
                default=_json_default
            )
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(
        'value of type {} is not JSON serializable'.format(
            type(value).__name__
        )
    )
