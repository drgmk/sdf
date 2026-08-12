"""Sampler-neutral scientific analysis of a completed fitting problem."""

from __future__ import annotations

from dataclasses import dataclass

import astropy.units as u
import numpy as np
from scipy.stats import truncnorm

from . import config as cfg
from . import fitting
from . import model


@dataclass(frozen=True)
class FitAnalysis:
    """Compact summaries plus arrays needed by plots and inspection."""

    summary: dict
    arrays: dict[str, np.ndarray]

    @classmethod
    def run(cls, problem, posterior):
        samples = posterior.samples
        if samples.shape[1] != problem.ndim:
            raise ValueError(
                'posterior sample width does not match fitting problem'
            )
        if tuple(posterior.parameter_names) != problem.parameter_names:
            raise ValueError(
                'posterior parameter names do not match fitting problem'
            )
        parameter_summary = {
            name: _summary(samples[:, index])
            for index, name in enumerate(problem.parameter_names)
        }
        median_parameters = np.median(samples, axis=0)
        residuals, wavelengths, filters = problem.residual(median_parameters)
        observed = fitting.concat_obs(problem.observations)
        obs_nel = observed[5].astype(int)

        prediction_samples = np.empty((len(samples), len(residuals)))
        component_samples = np.empty((
            len(samples), len(problem.candidate.components), len(residuals)
        ))
        for index, parameters in enumerate(samples):
            prediction, components = model.model_fluxes(
                problem.models, parameters, obs_nel,
                projection=problem.projection,
            )
            prediction_samples[index] = prediction
            component_samples[index] = components

        prediction = np.percentile(
            prediction_samples, (16.0, 50.0, 84.0), axis=0
        )
        component_prediction = np.percentile(
            component_samples, (16.0, 50.0, 84.0), axis=0
        )
        summary = {
            'parameters': parameter_summary,
            'chi_squared': float(np.sum(np.square(residuals))),
            'n_residuals': int(len(residuals)),
            'data_scopes': _data_scopes(problem.dataset, obs_nel),
            'components': _component_summaries(problem, samples),
            'targets': _target_summaries(problem, samples),
        }
        arrays = {
            'observed_flux': observed[0],
            'observed_uncertainty': observed[1],
            'observed_upper_limit': observed[2],
            'observed_ignore': observed[3],
            'wavelength_micron': wavelengths,
            'filter': filters,
            'prediction_p16': prediction[0],
            'prediction_median': prediction[1],
            'prediction_p84': prediction[2],
            'component_prediction_p16': component_prediction[0],
            'component_prediction_median': component_prediction[1],
            'component_prediction_p84': component_prediction[2],
            'residual': residuals,
        }
        return cls(summary=summary, arrays=arrays)


def _component_summaries(problem, samples):
    summaries = []
    offset = 0
    for component, component_models, plot_models in zip(
            problem.candidate.components, problem.models,
            problem.plot_models):
        parameter_names = tuple(component_models[0].parameters) + ('norm',)
        width = len(parameter_names)
        component_samples = samples[:, offset:offset + width]
        summaries.append({
            'component_id': component.component_id,
            'owner_key': component.owner_key,
            'role': component.role,
            'model_name': component.model_name,
            'parameters': {
                name: _summary(component_samples[:, index])
                for index, name in enumerate(parameter_names)
            },
        })
        offset += width
    return summaries


def _target_summaries(problem, samples):
    targets = {target.key: target for target in problem.dataset.targets}
    summaries = {
        key: {
            'main_id': target.main_id,
            'sdbid': target.sdbid,
            'stellar_components': [],
            'disk_components': [],
        }
        for key, target in targets.items()
    }
    parallax_samples = {
        key: _parallax_samples(target, len(samples))
        for key, target in targets.items()
    }
    for key, values in parallax_samples.items():
        if values is not None:
            summaries[key]['parallax_arcsec'] = _summary(values)
    layouts = []
    offset = 0
    for component, component_models, plot_models in zip(
            problem.candidate.components, problem.models,
            problem.plot_models):
        parameters = tuple(component_models[0].parameters) + ('norm',)
        width = len(parameters)
        component_samples = samples[:, offset:offset + width]
        layouts.append((component, plot_models, parameters,
                        component_samples))
        offset += width

    stellar_luminosity = {
        key: np.zeros(len(samples)) for key in targets
    }
    for component, _, parameters, component_samples in layouts:
        if component.role == 'stellar':
            derived, luminosity_at_1pc = _stellar_derived(
                component.model_name, parameters, component_samples,
                parallax_samples[component.owner_key]
            )
            derived['component_id'] = component.component_id
            summaries[component.owner_key]['stellar_components'].append(derived)
            if luminosity_at_1pc is not None:
                stellar_luminosity[component.owner_key] += luminosity_at_1pc

    for component, component_models, parameters, component_samples in layouts:
        if component.role == 'disk':
            derived = _disk_derived(
                component_models, component.model_name, parameters,
                component_samples,
                stellar_luminosity[component.owner_key],
                parallax_samples[component.owner_key],
            )
            derived['component_id'] = component.component_id
            summaries[component.owner_key]['disk_components'].append({
                **derived,
            })
    return summaries


def _stellar_derived(model_name, parameters, samples, parallax_samples):
    values = {'model_name': model_name}
    if 'norm' not in parameters:
        return values, None
    norm = samples[:, parameters.index('norm')]
    radius_at_1pc_m = np.sqrt(cfg.ssr * 10**norm / np.pi) * u.pc.to(u.m)
    values['radius_at_1pc_m'] = _summary(radius_at_1pc_m)

    if 'Teff' in parameters:
        temperature = samples[:, parameters.index('Teff')]
        luminosity_at_1pc = (
            4 * np.pi * radius_at_1pc_m**2 * temperature**4 *
            5.670374419e-8 / u.L_sun.to(u.W)
        )
        values['luminosity_at_1pc_lsun'] = _summary(luminosity_at_1pc)
    else:
        luminosity_at_1pc = None

    if parallax_samples is not None:
        values['radius_rsun'] = _summary(
            radius_at_1pc_m / parallax_samples / u.R_sun.to(u.m)
        )
        if luminosity_at_1pc is not None:
            values['luminosity_lsun'] = _summary(
                luminosity_at_1pc / parallax_samples**2
            )
    return values, luminosity_at_1pc


def _disk_derived(component_models, model_name, parameters, samples,
                  owner_luminosity_at_1pc, parallax_samples):
    values = {'model_name': model_name}
    spectral_model = next(
        (value for value in component_models
         if isinstance(value, model.SpecModel)),
        None,
    )
    if spectral_model is None:
        return values

    frequency_hz = u.micron.to(
        u.Hz, equivalencies=u.spectral()
    ) / spectral_model.wavelength
    order = np.argsort(frequency_hz)
    luminosity_at_1pc = np.empty(len(samples))
    for index, parameters_one in enumerate(samples):
        flux_jy = spectral_model.fnujy(parameters_one)
        irradiance = 1e-26 * np.trapezoid(
            flux_jy[order], frequency_hz[order]
        )
        luminosity_at_1pc[index] = (
            irradiance * 4 * np.pi * u.pc.to(u.m)**2 /
            u.L_sun.to(u.W)
        )
    values['luminosity_at_1pc_lsun'] = _summary(luminosity_at_1pc)

    has_stellar_luminosity = np.all(owner_luminosity_at_1pc > 0)
    if has_stellar_luminosity:
        values['fractional_luminosity'] = _summary(
            luminosity_at_1pc / owner_luminosity_at_1pc
        )

    temperature = None
    if 'Temp' in parameters:
        temperature = samples[:, parameters.index('Temp')]
    elif 'log_Temp' in parameters:
        temperature = 10**samples[:, parameters.index('log_Temp')]
    if temperature is not None:
        values['temperature_k'] = _summary(temperature)
        if has_stellar_luminosity and parallax_samples is not None:
            luminosity = owner_luminosity_at_1pc / parallax_samples**2
            values['blackbody_radius_au'] = _summary(
                np.sqrt(luminosity) * (278.3 / temperature)**2
            )
    return values


def _parallax_samples(target, count):
    parallax_mas = target.keywords.get('plx_value')
    if parallax_mas is None or float(parallax_mas) <= 0:
        return None
    mean = float(parallax_mas) / 1000.0
    error_mas = target.keywords.get('plx_err')
    if error_mas is None or float(error_mas) <= 0:
        return np.repeat(mean, count)
    scale = float(error_mas) / 1000.0
    quantiles = (np.arange(count, dtype=float) + 0.5) / count
    return truncnorm.ppf(
        quantiles, -mean / scale, np.inf, loc=mean, scale=scale
    )


def _data_scopes(dataset, obs_nel):
    scopes = []
    data_index = 0
    element_start = 0
    for observation in dataset.observations:
        for value in observation.data:
            element_stop = element_start + int(obs_nel[data_index])
            scopes.append({
                'observation_id': observation.observation_id,
                'contributors': list(observation.contributor_keys),
                'kind': type(value).__name__,
                'element_start': element_start,
                'element_stop': element_stop,
            })
            data_index += 1
            element_start = element_stop
    return scopes


def _summary(values):
    lower, median, upper = np.percentile(values, (16.0, 50.0, 84.0))
    return {
        'median': float(median),
        'minus': float(median - lower),
        'plus': float(upper - median),
    }
