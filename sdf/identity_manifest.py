"""Search-only target identities carried by an SDB fit package."""

from __future__ import annotations

from pathlib import Path

import yaml


IDENTITY_MANIFEST_VERSION = 1


def read_identity_manifest(path):
    """Return exported aliases keyed by SDB ID."""

    path = Path(path)
    try:
        payload = yaml.safe_load(path.read_text(encoding='utf-8'))
    except (OSError, yaml.YAMLError) as error:
        raise ValueError('could not read {}: {}'.format(path, error)) from error
    if not isinstance(payload, dict):
        raise ValueError('identities.yml must contain a mapping')
    unknown = set(payload) - {'version', 'targets'}
    if unknown:
        raise ValueError(
            'identities.yml contains unknown keys: {}'.format(
                ', '.join(sorted(unknown))
            )
        )
    if payload.get('version') != IDENTITY_MANIFEST_VERSION:
        raise ValueError(
            'unsupported identities.yml version {}'.format(
                payload.get('version')
            )
        )
    rows = payload.get('targets')
    if not isinstance(rows, dict):
        raise ValueError('identities.yml targets must be a mapping')

    result = {}
    for raw_sdbid, row in rows.items():
        sdbid = _nonempty_text(raw_sdbid, 'target SDB ID')
        if not isinstance(row, dict):
            raise ValueError(
                'identity for {} must be a mapping'.format(sdbid)
            )
        unknown = set(row) - {'aliases'}
        if unknown:
            raise ValueError(
                'identity for {} contains unknown keys: {}'.format(
                    sdbid, ', '.join(sorted(unknown))
                )
            )
        aliases = row.get('aliases')
        if not isinstance(aliases, list):
            raise ValueError(
                'aliases for {} must be a list'.format(sdbid)
            )
        result[sdbid] = tuple(dict.fromkeys(
            _nonempty_text(value, 'alias for {}'.format(sdbid))
            for value in aliases
        ))
    return result


def _nonempty_text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('{} must be non-empty text'.format(label))
    return value.strip()
