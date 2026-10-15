"""Sampler-neutral scientific analysis of a completed fitting problem."""

from __future__ import annotations

from dataclasses import dataclass

import astropy.units as u
import numpy as np
from scipy.stats import truncnorm

from . import config as cfg
from . import fitting
from . import model

ANALYSIS_SCHEMA = 4
PREDICTION_SAMPLES_MAX = 256


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

        # Spectrum calibration parameters multiply the observed segment in
        # the likelihood. Retain that fitted representation for plots and
        # residual inspection, while preserving the native arrays separately.
        spectrum_norm = np.take(
            np.append(median_parameters, 1.0), observed[4]
        )
        observed_flux = observed[0] * spectrum_norm
        observed_uncertainty = observed[1] * spectrum_norm

        prediction_indices = _prediction_indices(len(samples))
        prediction_samples = np.empty((
            len(prediction_indices), len(residuals)
        ))
        component_samples = np.empty((
            len(prediction_indices), len(problem.candidate.components),
            len(residuals)
        ))
        for index, sample_index in enumerate(prediction_indices):
            parameters = samples[sample_index]
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
        spectrum_wavelength, component_spectrum = _median_component_spectra(
            problem, median_parameters
        )
        target_summaries, derived_parameters, derived_samples = \
            _target_summaries(problem, samples)
        summary = {
            'schema': ANALYSIS_SCHEMA,
            'parameters': parameter_summary,
            'chi_squared': float(np.sum(np.square(residuals))),
            'n_residuals': int(len(residuals)),
            'prediction_samples': int(len(prediction_indices)),
            'data_scopes': _data_scopes(problem.dataset, obs_nel),
            'components': _component_summaries(problem, samples),
            'targets': target_summaries,
            'derived_parameters': derived_parameters,
        }
        arrays = {
            'observed_flux': observed_flux,
            'observed_uncertainty': observed_uncertainty,
            'observed_flux_native': observed[0],
            'observed_uncertainty_native': observed[1],
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
            'spectrum_wavelength_micron': spectrum_wavelength,
            'component_spectrum_median': component_spectrum,
            'residual': residuals,
            'derived_samples': derived_samples,
        }
        return cls(summary=summary, arrays=arrays)


def _median_component_spectra(problem, parameters):
    wavelength = np.asarray(cfg.models['default_wave'], dtype=float)
    spectra = np.zeros((len(problem.candidate.components), len(wavelength)))
    offset = 0
    for index, component_models in enumerate(problem.plot_models):
        spectral_model = next(
            (value for value in component_models
             if isinstance(value, model.SpecModel)),
            None,
        )
        if spectral_model is None:
            raise ValueError(
                'component {} has no spectral model'.format(index)
            )
        width = len(spectral_model.parameters) + 1
        if (len(spectral_model.wavelength) != len(wavelength) or
                not np.array_equal(spectral_model.wavelength, wavelength)):
            spectral_model = spectral_model.copy()
            spectral_model.interp_to_wavelengths(wavelength)
        spectra[index] = spectral_model.fnujy(
            parameters[offset:offset + width]
        )
        offset += width
    return wavelength, spectra


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
    derived_parameters = []
    derived_columns = []
    parallax_groups = []
    for key, values in parallax_samples.items():
        if values is not None:
            summaries[key]['parallax_arcsec'] = _summary(values)
            group = next((
                item for item in parallax_groups
                if np.array_equal(item['values'], values)
            ), None)
            if group is None:
                parallax_groups.append({
                    'values': values, 'owner_keys': [key],
                })
            else:
                group['owner_keys'].append(key)
    for group in parallax_groups:
        owner_keys = group['owner_keys']
        _add_derived_parameter(
            derived_parameters, derived_columns, group['values'],
            name='{}.parallax_arcsec'.format('+'.join(owner_keys)),
            owner_key=owner_keys[0], owner_keys=owner_keys,
            component_id=None, role='astrometry',
            label='parallax', unit='arcsec',
        )
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
            derived, distributions, luminosity_at_1pc = _stellar_derived(
                component.model_name, parameters, component_samples,
                parallax_samples[component.owner_key]
            )
            derived['component_id'] = component.component_id
            summaries[component.owner_key]['stellar_components'].append(derived)
            for quantity, values in distributions.items():
                label, unit = _DERIVED_LABELS[quantity]
                _add_derived_parameter(
                    derived_parameters, derived_columns, values,
                    name='{}.{}'.format(component.component_id, quantity),
                    owner_key=component.owner_key,
                    component_id=component.component_id, role='stellar',
                    label=label, unit=unit,
                )
            if luminosity_at_1pc is not None:
                stellar_luminosity[component.owner_key] += luminosity_at_1pc

    for component, component_models, parameters, component_samples in layouts:
        if component.role == 'disk':
            derived, distributions = _disk_derived(
                component_models, component.model_name, parameters,
                component_samples,
                stellar_luminosity[component.owner_key],
                parallax_samples[component.owner_key],
            )
            derived['component_id'] = component.component_id
            summaries[component.owner_key]['disk_components'].append({
                **derived,
            })
            for quantity, values in distributions.items():
                label, unit = _DERIVED_LABELS[quantity]
                _add_derived_parameter(
                    derived_parameters, derived_columns, values,
                    name='{}.{}'.format(component.component_id, quantity),
                    owner_key=component.owner_key,
                    component_id=component.component_id, role='disk',
                    label=label, unit=unit,
                )
    derived_samples = (
        np.column_stack(derived_columns)
        if derived_columns else np.empty((len(samples), 0))
    )
    return summaries, derived_parameters, derived_samples


def _stellar_derived(model_name, parameters, samples, parallax_samples):
    values = {'model_name': model_name}
    distributions = {}
    if 'norm' not in parameters:
        return values, distributions, None
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
        radius_rsun = (
            radius_at_1pc_m / parallax_samples / u.R_sun.to(u.m)
        )
        values['radius_rsun'] = _summary(radius_rsun)
        distributions['radius_rsun'] = radius_rsun
        if luminosity_at_1pc is not None:
            luminosity_lsun = luminosity_at_1pc / parallax_samples**2
            values['luminosity_lsun'] = _summary(luminosity_lsun)
            distributions['luminosity_lsun'] = luminosity_lsun
    return values, distributions, luminosity_at_1pc


def _disk_derived(component_models, model_name, parameters, samples,
                  owner_luminosity_at_1pc, parallax_samples):
    values = {'model_name': model_name}
    distributions = {}
    spectral_model = next(
        (value for value in component_models
         if isinstance(value, model.SpecModel)),
        None,
    )
    if spectral_model is None:
        return values, distributions

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
    if parallax_samples is not None:
        luminosity_lsun = luminosity_at_1pc / parallax_samples**2
        values['luminosity_lsun'] = _summary(luminosity_lsun)
        distributions['disk_luminosity_lsun'] = luminosity_lsun

    has_stellar_luminosity = np.all(owner_luminosity_at_1pc > 0)
    if has_stellar_luminosity:
        fractional_luminosity = (
            luminosity_at_1pc / owner_luminosity_at_1pc
        )
        values['fractional_luminosity'] = _summary(fractional_luminosity)
        distributions['fractional_luminosity'] = fractional_luminosity

    temperature = None
    if 'Temp' in parameters:
        temperature = samples[:, parameters.index('Temp')]
    elif 'log_Temp' in parameters:
        temperature = 10**samples[:, parameters.index('log_Temp')]
    if temperature is not None:
        values['temperature_k'] = _summary(temperature)
        distributions['temperature_k'] = temperature
        if has_stellar_luminosity and parallax_samples is not None:
            luminosity = owner_luminosity_at_1pc / parallax_samples**2
            blackbody_radius_au = (
                np.sqrt(luminosity) * (278.3 / temperature)**2
            )
            values['blackbody_radius_au'] = _summary(blackbody_radius_au)
            distributions['blackbody_radius_au'] = blackbody_radius_au
    return values, distributions


_DERIVED_LABELS = {
    'radius_rsun': ('stellar radius', 'R_sun'),
    'luminosity_lsun': ('stellar luminosity', 'L_sun'),
    'disk_luminosity_lsun': ('disk luminosity', 'L_sun'),
    'fractional_luminosity': ('disk / stellar luminosity', ''),
    'temperature_k': ('disk temperature', 'K'),
    'blackbody_radius_au': ('blackbody radius', 'au'),
}


def _add_derived_parameter(parameters, columns, values, *, name, owner_key,
                           component_id, role, label, unit, owner_keys=None):
    """Keep a plotted derived series aligned with its descriptive metadata."""

    parameters.append({
        'name': name,
        'owner_key': owner_key,
        'owner_keys': list(owner_keys or (owner_key,)),
        'component_id': component_id,
        'role': role,
        'label': label,
        'unit': unit,
    })
    columns.append(np.asarray(values, dtype=float))


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
            scope = {
                'observation_id': observation.observation_id,
                'contributors': list(observation.contributor_keys),
                'kind': type(value).__name__,
                'element_start': element_start,
                'element_stop': element_stop,
            }
            if observation.kind == 'spectrum':
                scope.update({
                    'product_id': observation.product_id,
                    'segment_id': observation.segment_id,
                    'instrument': getattr(value, 'instrument', None),
                    'resolving_power': getattr(
                        value, 'resolving_power', None
                    ),
                })
            scopes.append(scope)
            data_index += 1
            element_start = element_stop
    return scopes


def _prediction_indices(count):
    """Use a deterministic bounded posterior subset for prediction bands."""

    if count <= PREDICTION_SAMPLES_MAX:
        return np.arange(count, dtype=int)
    return np.linspace(
        0, count - 1, PREDICTION_SAMPLES_MAX, dtype=int
    )


def _summary(values):
    lower, median, upper = np.percentile(values, (16.0, 50.0, 84.0))
    return {
        'median': float(median),
        'minus': float(median - lower),
        'plus': float(upper - median),
    }
