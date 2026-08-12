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


@dataclass(frozen=True)
class StoredCandidate:
    """One candidate artifact loaded from a package library."""

    path: Path
    metadata: dict
    posterior: PosteriorRun

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

    def __post_init__(self):
        displayed = tuple(self.displayed)
        if not displayed:
            raise ValueError('fit selection must display at least one candidate')
        if self.default not in displayed:
            raise ValueError('default candidate must also be displayed')
        if len(set(displayed)) != len(displayed):
            raise ValueError('displayed candidates must be unique')
        object.__setattr__(self, 'displayed', displayed)


class CandidateLibrary:
    """Filesystem store for completed candidate fits in one package."""

    def __init__(self, package_directory):
        self.package_directory = Path(package_directory).expanduser().resolve()
        self.fits_directory = self.package_directory / 'fits'

    def commit(self, dataset, candidate, posterior, *, analysis=None,
               policy=None, provenance=None, fit_state=None,
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

        arrays = {'samples': posterior.samples}
        if posterior.log_likelihoods is not None:
            arrays['log_likelihoods'] = posterior.log_likelihoods
        if posterior.weights is not None:
            arrays['weights'] = posterior.weights
        _atomic_npz(posterior_path, arrays)

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
            },
            'analysis': analysis or {},
            'policy': policy or {},
            'provenance': provenance or {},
        }
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
            )
        return StoredCandidate(candidate_directory, metadata, posterior)

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

    def write_selection(self, selection):
        """Atomically update the cheap mutable display/default choice."""

        available = {candidate.candidate_id for candidate in self.candidates()}
        referenced = set(selection.displayed) | {selection.default}
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
        )


def input_fingerprint(dataset):
    """Hash portable package inputs rather than their absolute locations."""

    paths = {Path(observation.source) for observation in dataset.observations}
    for parent in {path.parent for path in paths}:
        topology = parent / 'joint-fit.yml'
        if topology.is_file():
            paths.add(topology)
    digest = sha256()
    for path in sorted(paths, key=lambda value: value.name):
        digest.update(path.name.encode('utf-8'))
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
