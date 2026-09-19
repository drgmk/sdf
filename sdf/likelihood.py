"""Compiled fitting problems independent of any particular sampler."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import model


@dataclass
class FitProblem:
    """Numerical observations, models, priors, and likelihood for a candidate."""

    dataset: object
    candidate: object
    observations: tuple
    models: tuple
    plot_models: tuple
    model_info: dict
    projection: np.ndarray

    @classmethod
    def compile(cls, dataset, candidate):
        observations = dataset.data
        models, plot_models = model.get_models(
            observations, candidate.model_names
        )
        return cls(
            dataset=dataset,
            candidate=candidate,
            observations=observations,
            models=models,
            plot_models=plot_models,
            model_info=model.models_info(models),
            projection=dataset.component_projection(candidate),
        )

    @classmethod
    def from_models(cls, observations, models):
        """Wrap already loaded legacy observations/models."""

        return cls(
            dataset=None,
            candidate=None,
            observations=tuple(observations),
            models=tuple(models),
            plot_models=(),
            model_info=model.models_info(models),
            projection=np.ones(
                (len(observations), len(models)), dtype=bool
            ),
        )

    @property
    def ndim(self):
        return self.model_info['ndim']

    @property
    def parameter_ranges(self):
        return tuple(self.model_info['p_rng'])

    @property
    def parameter_names(self):
        """Return deterministic owner-qualified parameter labels."""

        if self.candidate is None:
            return tuple(self.model_info['parameters'])
        names = []
        for component, component_models in zip(
                self.candidate.components, self.models):
            model_parameters = tuple(component_models[0].parameters) + ('norm',)
            names.extend(
                '{}.{}'.format(component.component_id, parameter)
                for parameter in model_parameters
            )
        for observation in self.dataset.observations:
            if observation.kind == 'spectrum' and observation.segment_id:
                names.extend(
                    'observation:{}.norm'.format(observation.observation_id)
                    for _ in observation.spectra
                )
            else:
                names.extend(
                    'observation:{}.spectrum:{}.norm'.format(
                        observation.observation_id, index + 1
                    )
                    for index in range(len(observation.spectra))
                )
        if len(names) != self.ndim:
            raise ValueError(
                'qualified parameter count does not match fitting dimension'
            )
        return tuple(names)

    def transform_prior(self, unit_values):
        """Transform unit-cube values to the model parameter ranges."""

        values = np.asarray(unit_values, dtype=float).copy()
        for i, (lower, upper) in enumerate(self.parameter_ranges):
            values[i] = lower + values[i] * (upper - lower)
        return values

    def residual(self, parameters):
        from . import fitting

        return fitting.residual(
            parameters, self.observations, self.models,
            projection=self.projection
        )

    def chi_squared(self, parameters):
        residuals, _, _ = self.residual(parameters)
        return np.sum(np.square(residuals))

    def log_likelihood(self, parameters):
        chi_squared = self.chi_squared(parameters)
        if np.isfinite(chi_squared):
            return -0.5 * chi_squared
        return -np.inf
