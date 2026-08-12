"""Model-family metadata and per-target fitting policy."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

import yaml

from . import config as cfg
from . import utils


@dataclass(frozen=True)
class ModelFamily:
    name: str
    role: str


class ModelFamilyRegistry:
    """Exact model-family classification during the config transition."""

    def __init__(self, families):
        families = tuple(families)
        self._families = {family.name: family for family in families}
        if len(self._families) != len(families):
            raise utils.SdfError('model family names must be unique')

    @classmethod
    def from_config(cls):
        names = set(cfg.models['star']) | set(cfg.models['disk'])
        overlap = set(cfg.models['star']) & set(cfg.models['disk'])
        if overlap:
            raise utils.SdfError(
                'models classified as both stellar and disk: {}'.format(
                    ', '.join(sorted(overlap))
                )
            )
        return cls(
            ModelFamily(
                name=name,
                role='stellar' if name in cfg.models['star'] else 'disk',
            )
            for name in sorted(names)
        )

    def require(self, name, role=None):
        try:
            family = self._families[name]
        except KeyError as error:
            raise utils.SdfError(
                'unknown model family {}'.format(name)
            ) from error
        if role is not None and family.role != role:
            raise utils.SdfError(
                '{} is a {} model, not {}'.format(
                    name, family.role, role
                )
            )
        return family

    def names(self, role=None):
        return tuple(
            name for name, family in sorted(self._families.items())
            if role is None or family.role == role
        )


@dataclass(frozen=True)
class TargetModelPolicy:
    """Resolved model search limits for one physical contributor."""

    target_key: str
    stellar_models: tuple[str, ...]
    disk_models: tuple[str, ...] = ('modbb_disk_r',)
    max_disks: int = 2
    search_mode: str = 'evidence-guided'
    reddening: bool = False
    origin: str = 'automatic'

    def __post_init__(self):
        object.__setattr__(self, 'stellar_models', tuple(self.stellar_models))
        object.__setattr__(self, 'disk_models', tuple(self.disk_models))
        if not self.stellar_models:
            raise utils.SdfError('at least one stellar model is required')
        if self.max_disks < 0:
            raise utils.SdfError('max_disks cannot be negative')
        if self.max_disks and not self.disk_models:
            raise utils.SdfError(
                'disk models are required when max_disks is non-zero'
            )
        if self.search_mode not in {
                'evidence-guided', 'exhaustive', 'exact', 'posterior-only'}:
            raise utils.SdfError(
                'unsupported search mode {}'.format(self.search_mode)
            )


def automatic_policy(target, registry=None):
    """Infer conservative defaults from the target spectral type."""

    registry = registry or ModelFamilyRegistry.from_config()
    spectral_type = target.keywords.get('sp_type')
    spectral_type = '' if spectral_type is None else str(spectral_type).strip()
    upper = spectral_type.upper()

    if upper.startswith(('DA', 'DB', 'DC')):
        stellar_models = ('koester_wd',)
        disk_models = ('bb_disk_r',)
    elif upper.startswith(('L', 'T', 'Y')):
        stellar_models = ('phoenix_cool',)
        disk_models = ('modbb_disk_r',)
    elif _is_late_m(upper):
        stellar_models = ('phoenix_cool', 'phoenix_sol')
        disk_models = ('modbb_disk_r',)
    else:
        stellar_models = ('phoenix_sol',)
        disk_models = ('modbb_disk_r',)

    for name in stellar_models:
        registry.require(name, 'stellar')
    for name in disk_models:
        registry.require(name, 'disk')
    return TargetModelPolicy(
        target_key=target.key,
        stellar_models=stellar_models,
        disk_models=disk_models,
    )


def resolve_package_policies(dataset, package_directory, registry=None):
    """Apply optional portable `sdf-fit.yml` over automatic defaults."""

    registry = registry or ModelFamilyRegistry.from_config()
    policies = {
        target.key: automatic_policy(target, registry=registry)
        for target in dataset.targets
    }
    path = Path(package_directory) / 'sdf-fit.yml'
    if not path.is_file():
        return tuple(policies[key] for key in sorted(policies))
    payload = _read_policy(path)
    by_name = {
        _normalize_name(target.main_id): target
        for target in dataset.targets if target.main_id is not None
    }
    unknown = set(payload) - set(by_name)
    if unknown:
        raise utils.SdfError(
            'sdf-fit.yml contains targets not in this dataset: {}'.format(
                ', '.join(sorted(unknown))
            )
        )
    for name, values in payload.items():
        target = by_name[name]
        current = policies[target.key]
        stellar_models = tuple(
            values.get('stellar_models', current.stellar_models)
        )
        disk_models = tuple(values.get('disk_models', current.disk_models))
        for model_name in stellar_models:
            registry.require(model_name, 'stellar')
        for model_name in disk_models:
            registry.require(model_name, 'disk')
        policies[target.key] = replace(
            current,
            stellar_models=stellar_models,
            disk_models=disk_models,
            max_disks=values.get('max_disks', current.max_disks),
            search_mode=values.get('search_mode', current.search_mode),
            reddening=values.get('reddening', current.reddening),
            origin='package',
        )
    return tuple(policies[key] for key in sorted(policies))


def _is_late_m(spectral_type):
    value = spectral_type[1:] if spectral_type.startswith('D') else spectral_type
    return len(value) > 1 and value.startswith('M') and value[1] in '56789'


def _read_policy(path):
    try:
        payload = yaml.safe_load(path.read_text(encoding='utf-8'))
    except (OSError, yaml.YAMLError) as error:
        raise utils.SdfError(
            'could not read {}: {}'.format(path, error)
        ) from error
    if not isinstance(payload, dict) or set(payload) != {'targets'}:
        raise utils.SdfError('sdf-fit.yml must contain only a targets mapping')
    if not isinstance(payload['targets'], dict):
        raise utils.SdfError('sdf-fit targets must be a mapping')
    allowed = {
        'stellar_models', 'disk_models', 'max_disks', 'search_mode',
        'reddening'
    }
    normalized = {}
    for target_name, values in payload['targets'].items():
        if not isinstance(target_name, str) or not isinstance(values, dict):
            raise utils.SdfError('target policies must be named mappings')
        extra = set(values) - allowed
        if extra:
            raise utils.SdfError(
                'unknown policy fields for {}: {}'.format(
                    target_name, ', '.join(sorted(extra))
                )
            )
        for field in ('stellar_models', 'disk_models'):
            if field in values and (
                    not isinstance(values[field], list) or
                    any(not isinstance(value, str) for value in values[field])):
                raise utils.SdfError('{} must be a list of model names'.format(field))
        if 'max_disks' in values and (
                not isinstance(values['max_disks'], int) or
                isinstance(values['max_disks'], bool)):
            raise utils.SdfError('max_disks must be an integer')
        normalized_name = _normalize_name(target_name)
        if normalized_name in normalized:
            raise utils.SdfError(
                'duplicate normalized target name {}'.format(normalized_name)
            )
        normalized[normalized_name] = values
    return normalized


def _normalize_name(value):
    return ' '.join(value.split())
