"""Regenerable posterior-distribution products for retained candidates."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile

import numpy as np


MINIMUM_SAMPLES = 16


def ensure_candidate_diagnostic_plots(candidate, *, force=False):
    """Write fitted and derived corner plots when meaningful samples exist."""

    products = {}
    fitted = candidate.path / 'corner.png'
    if force or not fitted.is_file():
        if _enough_samples(candidate.posterior.samples):
            _write_corner(
                candidate.posterior.samples,
                _parameter_labels(candidate),
                fitted,
                weights=candidate.posterior.weights,
            )
    if fitted.is_file():
        products['fitted'] = fitted

    derived = candidate.path / 'distributions.png'
    samples = candidate.analysis_arrays.get('derived_samples')
    metadata = candidate.metadata.get('analysis', {}).get(
        'derived_parameters', ()
    )
    if (samples is not None and np.asarray(samples).ndim == 2 and
            np.asarray(samples).shape[1] == len(metadata) and metadata):
        if force or not derived.is_file():
            if _enough_samples(samples):
                _write_corner(
                    samples, _derived_labels(candidate, metadata), derived
                )
    elif force and derived.is_file():
        derived.unlink()
    if derived.is_file():
        products['derived'] = derived
    return products


def _enough_samples(samples):
    samples = np.asarray(samples)
    return (
        samples.ndim == 2 and samples.shape[1] > 0 and
        samples.shape[0] >= MINIMUM_SAMPLES
    )


def _write_corner(samples, labels, path, *, weights=None):
    # Keep the plotting stack out of catalog/server startup.  Fit workers and
    # first-time page renders pay this cost only when a product is needed.
    import matplotlib
    matplotlib.use('agg')
    import corner
    import matplotlib.pyplot as plt

    samples = np.asarray(samples, dtype=float)
    finite = np.all(np.isfinite(samples), axis=1)
    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        finite &= np.isfinite(weights) & (weights >= 0)
        weights = weights[finite]
        if not np.any(weights > 0):
            weights = None
    samples = samples[finite]
    if not _enough_samples(samples):
        return None

    ranges = [_plot_range(samples[:, index]) for index in range(samples.shape[1])]
    figure = corner.corner(
        samples,
        labels=labels,
        weights=weights,
        range=ranges,
        show_titles=True,
        title_fmt='.3g',
        max_n_ticks=4,
        label_kwargs={'fontsize': 9},
        title_kwargs={'fontsize': 9},
    )
    figure.suptitle(
        'Fitted parameters' if path.name == 'corner.png'
        else 'Derived quantities',
        fontsize=13,
    )
    _atomic_savefig(figure, path)
    plt.close(figure)
    return path


def _plot_range(values):
    lower = float(np.min(values))
    upper = float(np.max(values))
    if upper > lower:
        padding = 0.025 * (upper - lower)
        return lower - padding, upper + padding
    padding = max(abs(lower) * 1e-6, 1e-9)
    return lower - padding, upper + padding


def _parameter_labels(candidate):
    components = candidate.metadata.get('components', ())
    component_by_id = {
        component['component_id']: component for component in components
    }
    labels = []
    for name in candidate.posterior.parameter_names:
        match = next((
            (component_id, component)
            for component_id, component in component_by_id.items()
            if name.startswith(component_id + '.')
        ), None)
        if match is None:
            labels.append(_observation_parameter_label(name))
            continue
        component_id, component = match
        parameter = name[len(component_id) + 1:]
        labels.append('{} {} · {}'.format(
            _owner_label(candidate, component['owner_key']),
            _component_label(component, components),
            _PARAMETER_LABELS.get(parameter, parameter),
        ))
    return labels


def _derived_labels(candidate, parameters):
    components = candidate.metadata.get('components', ())
    component_by_id = {
        component['component_id']: component for component in components
    }
    labels = []
    for parameter in parameters:
        owner_keys = parameter.get('owner_keys', (parameter['owner_key'],))
        owner = (
            'system' if len(owner_keys) > 1 else
            _owner_label(candidate, owner_keys[0])
        )
        component = component_by_id.get(parameter.get('component_id'))
        part = '' if component is None else ' ' + _component_label(
            component, components
        )
        unit = parameter.get('unit')
        label = '{}{} · {}'.format(owner, part, parameter['label'])
        if unit:
            label += ' ({})'.format(_UNIT_LABELS.get(unit, unit))
        labels.append(label)
    return labels


def _owner_label(candidate, owner_key):
    target = candidate.metadata.get('analysis', {}).get(
        'targets', {}
    ).get(owner_key, {})
    return target.get('main_id') or target.get('sdbid') or owner_key


def _component_label(component, components):
    role = component['role']
    siblings = [
        value for value in components
        if value['owner_key'] == component['owner_key'] and
        value['role'] == role
    ]
    symbol = 'star' if role == 'stellar' else 'disk'
    if len(siblings) > 1:
        symbol += ' {}'.format(component.get('ordinal', 1))
    return symbol


def _observation_parameter_label(name):
    if name.startswith('observation:') and name.endswith('.norm'):
        return '{} · scale'.format(name[len('observation:'):-len('.norm')])
    return name


def _atomic_savefig(figure, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        prefix='.' + path.stem + '-', suffix=path.suffix,
        dir=path.parent, delete=False,
    )
    temporary = Path(handle.name)
    handle.close()
    try:
        figure.savefig(
            temporary, dpi=135, bbox_inches='tight', facecolor='white'
        )
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


_PARAMETER_LABELS = {
    'Teff': 'T_eff (K)',
    'logg': 'log(g)',
    'norm': 'log normalisation',
    'log_Temp': 'log T_disk',
    'Temp': 'T_disk (K)',
    'log_lam0': 'log lambda_0',
    'beta': 'beta',
}

_UNIT_LABELS = {
    'R_sun': 'R_sun',
    'L_sun': 'L_sun',
}
