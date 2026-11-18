"""Provider-neutral links to products accompanying exported observations."""

from __future__ import annotations

from pathlib import Path

import yaml


RESOURCE_MANIFEST_VERSION = 1
_LIST_FIELDS = ('observation_ids', 'observation_sdbids', 'target_sdbids')


def read_resource_manifest(path):
    """Read and validate an optional SDB ``resources.yml`` sidecar."""

    path = Path(path)
    document = yaml.safe_load(path.read_text(encoding='utf-8'))
    if not isinstance(document, dict):
        raise ValueError('{} is not a resource manifest'.format(path))
    if document.get('version') != RESOURCE_MANIFEST_VERSION:
        raise ValueError(
            '{} has unsupported resource manifest version {}'.format(
                path, document.get('version')
            )
        )
    raw_resources = document.get('resources')
    if not isinstance(raw_resources, list):
        raise ValueError('{} resources must be an array'.format(path))

    resources = []
    for index, raw in enumerate(raw_resources):
        if not isinstance(raw, dict):
            raise ValueError(
                '{} resource {} must be an object'.format(path, index)
            )
        provider = _required_text(path, index, raw, 'provider')
        label = _required_text(path, index, raw, 'label')
        url = _required_text(path, index, raw, 'url')
        values = {
            field: _string_list(path, index, raw.get(field, ()), field)
            for field in _LIST_FIELDS
        }
        resources.append({
            'provider': provider,
            'label': label,
            'url': url,
            **values,
        })
    return tuple(resources)


def _required_text(path, index, row, field):
    value = row.get(field)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(
            '{} resource {} {} must be non-empty text'.format(
                path, index, field
            )
        )
    return value.strip()


def _string_list(path, index, value, field):
    numeric_ids = field == 'observation_ids'
    if not isinstance(value, (list, tuple)):
        raise ValueError(
            '{} resource {} {} must be an array of non-empty strings'.format(
                path, index, field
            )
        )
    cleaned = []
    for item in value:
        if numeric_ids and isinstance(item, int) and not isinstance(item, bool):
            item = str(item)
        if not isinstance(item, str) or not item.strip():
            raise ValueError(
                '{} resource {} {} must be an array of non-empty strings'.
                format(path, index, field)
            )
        cleaned.append(item.strip())
    return tuple(dict.fromkeys(cleaned))
