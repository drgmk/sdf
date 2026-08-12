"""Sampler-independent descriptions of model components and candidates."""

from __future__ import annotations

from dataclasses import dataclass

from . import config as cfg
from . import utils


@dataclass(frozen=True)
class ModelComponentSpec:
    """One numerical model component owned by one physical target."""

    component_id: str
    owner_key: str
    role: str
    model_name: str
    ordinal: int


@dataclass(frozen=True)
class ModelCandidate:
    """A complete model hypothesis to evaluate for one dataset."""

    components: tuple[ModelComponentSpec, ...]

    @classmethod
    def for_single_target(cls, target_key, model_names):
        """Describe a legacy tuple of model names with explicit ownership."""

        return cls.from_owned_models(
            (target_key, model_name) for model_name in tuple(model_names)
        )

    @classmethod
    def from_owned_models(cls, owned_models):
        """Describe an ordered set of `(target_key, model_name)` pairs."""

        counts = {}
        components = []
        for target_key, model_name in tuple(owned_models):
            if model_name in cfg.models['star']:
                role = 'stellar'
            elif model_name in cfg.models['disk']:
                role = 'disk'
            else:
                raise utils.SdfError(
                    "couldn't assign model {} to a stellar or disk role".
                    format(model_name)
                )
            key = (target_key, role, model_name)
            counts[key] = counts.get(key, 0) + 1
            ordinal = counts[key]
            components.append(ModelComponentSpec(
                component_id='{}:{}:{}:{}'.format(
                    target_key, role, model_name, ordinal
                ),
                owner_key=target_key,
                role=role,
                model_name=model_name,
                ordinal=ordinal,
            ))
        if not components:
            raise utils.SdfError('a model candidate needs at least one component')
        return cls(components=tuple(components))

    @property
    def model_names(self):
        return tuple(component.model_name for component in self.components)

    @property
    def candidate_id(self):
        """Return the familiar readable name for a single-owner candidate."""

        owners = {component.owner_key for component in self.components}
        if len(owners) == 1:
            return cfg.fitting['model_join'].join(self.model_names)
        return '+'.join(
            '{}-{}'.format(component.owner_key, component.model_name)
            for component in self.components
        )
