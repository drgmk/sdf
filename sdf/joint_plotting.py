"""Bokeh presentation for compact joint-fit candidate artifacts."""

from __future__ import annotations

from html import escape
from pathlib import Path

import numpy as np
from bokeh.embed import components
from bokeh.layouts import column, gridplot
from bokeh.models import (
    ColumnDataSource, Div, HoverTool, InlineStyleSheet, Range1d, Span,
    TabPanel, Tabs,
)
from bokeh.palettes import Category10, Category20
from bokeh.plotting import figure

from . import filter as filter_module
from .analytics import BB_Disk
from .candidate_library import CandidateLibrary
from .dataset import FitDataset, canonical_fit_source
from .diagnostic_plots import ensure_candidate_diagnostic_plots
from .model_policy import ModelFamilyRegistry, resolve_package_policies
from .sed_page import (
    add_operator_policy,
    candidate_page_rows,
    observation_page_rows,
    render_sed_page,
    system_page_rows,
    target_page_rows,
)


GOLDEN_RATIO = (1 + np.sqrt(5)) / 2
DATA_MARKER_SIZE = 12
MODEL_MARKER_SIZE = 14
UPPER_LIMIT_MARKER_SIZE = 10
EXCLUDED_COLOUR = '#9ca3af'
LEGEND_COLOUR = '#4c78a8'


def write_joint_results_page(package_directory, output=None, *, dataset=None,
                             inline=True, candidate_ids=None, operator=False,
                             policy_store=None, operator_api=None,
                             operator_token=None, operator_actor=None,
                             operator_actor_locked=False, navigation=None):
    """Write a standalone candidate-review page for one fitting package."""

    source = canonical_fit_source(package_directory)
    package_directory = source if source.is_dir() else source.parent
    if dataset is None:
        dataset = FitDataset.from_path(source)
    library = CandidateLibrary(package_directory)
    candidates = library.candidates()
    if not candidates:
        raise ValueError('joint package has no retained candidates to plot')

    try:
        selection = library.read_selection()
    except FileNotFoundError:
        selection = None
    ranked_ids = [
        candidate.candidate_id for candidate in sorted(
            candidates, key=_sortable_evidence, reverse=True,
        )
    ]
    public_default = (
        selection.default if selection is not None else ranked_ids[0]
    )
    recommended = (
        selection.recommended if selection is not None else ranked_ids[0]
    )
    if candidate_ids is not None:
        order = list(candidate_ids)
        if not order:
            raise ValueError('at least one candidate is required for preview')
        default = public_default
    else:
        if selection is not None:
            order = list(selection.displayed)
            default = selection.default
        else:
            order = ranked_ids
            default = order[0]

    by_id = {candidate.candidate_id: candidate for candidate in candidates}
    missing = set(order) - set(by_id)
    if missing:
        raise ValueError('unknown candidates for plot: {}'.format(
            ', '.join(sorted(missing))
        ))
    ordered = [by_id[candidate_id] for candidate_id in order]
    for candidate in ordered:
        ensure_candidate_diagnostic_plots(candidate)
    owner_labels = _owner_labels(ordered[0])
    owner_colours = _owner_colour_map(owner_labels)
    external_route_root = (
        None if navigation is None else navigation.get('external_route_root')
    )
    target_rows = target_page_rows(
        dataset, package_directory, external_route_root=external_route_root,
    )
    for target in target_rows:
        target['colour'] = owner_colours.get(target['key'], '#64748b')
    observation_rows = observation_page_rows(
        dataset, package_directory, external_route_root=external_route_root,
    )
    systems, component_labels = system_page_rows(
        target_rows, observation_rows
    )
    title = _system_title(systems, owner_labels, package_directory)

    panels = []
    for candidate in ordered:
        panels.append(TabPanel(
            child=_candidate_panel(
                candidate, dataset, owner_labels, component_labels,
                owner_colours,
                is_default=candidate.candidate_id == default,
            ),
            title=_candidate_shorthand_text(
                candidate, owner_labels, component_labels
            ),
        ))
    panels.append(_limits_tab(
        by_id[default], owner_labels, component_labels, owner_colours,
    ))
    active = next((
        index for index, candidate in enumerate(ordered)
        if candidate.candidate_id == default
    ), 0)
    tabs = Tabs(
        tabs=panels, active=active, sizing_mode='stretch_width',
        stylesheets=[InlineStyleSheet(css='''
          :host { --base-font: ui-sans-serif, system-ui, sans-serif; }
        ''')],
    )
    plot_script, plot_div = components(tabs)
    candidate_rows = candidate_page_rows(
        ordered, owner_labels, default, package_directory, dataset,
        displayed_ids=order,
        component_labels=component_labels, owner_colours=owner_colours,
        recommended_candidate=recommended,
    )
    available_rows = ()
    queued_jobs = ()
    if operator:
        registry = ModelFamilyRegistry.from_config()
        policies = resolve_package_policies(
            dataset, package_directory, registry=registry,
            policy_store=policy_store,
        )
        if operator_actor is None and policy_store is not None:
            operator_actor = policy_store.last_actor()
        add_operator_policy(
            target_rows, policies, registry, policy_store=policy_store,
            package_directory=package_directory, actor=operator_actor,
        )
        available_rows = candidate_page_rows(
            sorted(candidates, key=_sortable_evidence, reverse=True),
            owner_labels, default, package_directory, dataset,
            displayed_ids=order,
            component_labels=component_labels, owner_colours=owner_colours,
            recommended_candidate=recommended,
        )
        if policy_store is not None:
            queued_jobs = policy_store.jobs(
                package=package_directory, status='queued',
            )
    html = render_sed_page(
        title=title,
        plot_script=plot_script,
        plot_div=plot_div,
        targets=target_rows,
        candidates=candidate_rows,
        observations=observation_rows,
        systems=systems,
        default_candidate=default,
        inline=inline,
        operator=operator,
        available_fits=available_rows,
        operator_api=operator_api,
        operator_token=operator_token,
        operator_actor=operator_actor,
        operator_actor_locked=operator_actor_locked,
        queued_jobs=queued_jobs,
        navigation=navigation,
    )
    default_name = 'operator.html' if operator else 'index.html'
    output = (
        Path(output) if output is not None else package_directory / default_name
    )
    output.write_text(html, encoding='utf-8')
    return output


def _limits_tab(candidate, owner_labels, component_labels, owner_colours):
    """Build the system sensitivity view from the default retained fit."""

    owner_panels = []
    for owner_key in owner_labels:
        limits = _detection_limits(candidate, owner_key)
        if limits is None:
            continue
        label = (
            component_labels.get(owner_key) or owner_labels.get(owner_key) or
            owner_key
        )
        owner_panels.append(TabPanel(
            child=_limits_panel(
                limits, label, owner_colours.get(owner_key, '#285ea8')
            ),
            title=str(label),
        ))

    if not owner_panels:
        child = Div(
            text='''
<section style="padding:16px;font:13px ui-sans-serif,system-ui,sans-serif">
  <strong>Detection limits are not available</strong>
  <p style="color:#526071">A stellar luminosity and usable photometry beyond
  8 μm are needed to derive this view.</p>
</section>''',
            sizing_mode='stretch_width',
        )
    elif len(owner_panels) == 1:
        child = owner_panels[0].child
    else:
        child = Tabs(
            tabs=owner_panels, sizing_mode='stretch_width',
            stylesheets=[InlineStyleSheet(css='''
              :host { --base-font: ui-sans-serif, system-ui, sans-serif; }
            ''')],
        )
    return TabPanel(child=child, title='Limits')


def _limits_panel(limits, label, colour):
    """Render one physical contributor's fractional-luminosity limits."""

    radius_available = limits['radius_au'] is not None
    x = limits['radius_au'] if radius_available else limits['temperature_k']
    x_label = (
        'Blackbody radius (au)' if radius_available else
        'Disk temperature (K)'
    )
    fig = figure(
        title='Detection limits · {}'.format(label),
        x_axis_label=x_label, x_axis_type='log',
        y_axis_label='Fractional luminosity', y_axis_type='log',
        x_range=(float(np.min(x)), float(np.max(x))),
        y_range=(1e-7, 2e-1),
        tools='pan,wheel_zoom,box_zoom,reset,save',
        active_scroll='wheel_zoom', aspect_ratio=GOLDEN_RATIO,
        sizing_mode='stretch_width',
    )
    envelope_source = ColumnDataSource({
        'x': x,
        'limit': limits['envelope'],
        'top': np.repeat(2e-1, len(x)),
    })
    fig.varea(
        'x', 'limit', 'top', source=envelope_source,
        fill_color='#dbeafe', fill_alpha=0.42,
        legend_label='detectable',
    )

    colours = Category20[20]
    for index in range(limits['curves'].shape[1]):
        band = str(limits['filters'][index])
        basis = str(limits['basis'][index])
        band_colour = colours[index % len(colours)]
        source = ColumnDataSource({
            'x': x,
            'limit': limits['curves'][:, index],
            'band': np.repeat(band, len(x)),
            'basis': np.repeat(basis, len(x)),
        })
        renderer = fig.line(
            'x', 'limit', source=source, color=band_colour,
            line_width=2.2,
            line_dash='dashed' if basis == 'upper limit' else 'solid',
            muted_alpha=0.15, legend_label=band,
        )
        fig.add_tools(HoverTool(
            renderers=[renderer],
            tooltips=[('band', '@band'), ('limit from', '@basis')],
        ))

    for curve, curve_label, dash in limits['reference_curves']:
        fig.line(
            x, curve, color='#94a3b8', line_width=1.8,
            line_dash=dash, line_alpha=0.9, legend_label=curve_label,
        )

    for disk in limits['disk_samples']:
        disk_x = (
            disk['radius_au'] if radius_available else disk['temperature_k']
        )
        ok = (
            np.isfinite(disk_x) & (disk_x > 0) &
            np.isfinite(disk['fractional_luminosity']) &
            (disk['fractional_luminosity'] > 0)
        )
        if not np.any(ok):
            continue
        fig.scatter(
            disk_x[ok], disk['fractional_luminosity'][ok],
            marker='circle', size=5, color=colour, alpha=0.12,
            legend_label='fitted disk',
        )

    fig.grid.grid_line_alpha = 0.18
    fig.outline_line_color = '#d1d5db'
    if fig.legend:
        fig.legend.location = 'top_left'
        fig.legend.orientation = 'vertical'
        fig.legend.click_policy = 'mute'
        fig.legend.background_fill_alpha = 0.82
        fig.legend.border_line_alpha = 0.25
        fig.legend.label_text_font_size = '9pt'
        fig.legend.ncols = 2
        fig.legend.padding = 4
        fig.legend.spacing = 4

    note = '''
<section style="padding:8px 10px 2px;font:13px ui-sans-serif,system-ui,sans-serif">
  <div style="display:flex;align-items:baseline;gap:12px;flex-wrap:wrap">
    <strong style="font-size:15px">Fractional-luminosity limits</strong>
    <span style="color:#526071">based on the default stellar fit</span>
  </div>
  <p style="color:#526071;margin:3px 0 0">The blue region is detectable by at
  least one observed band. Coloured solid curves use 3σ measurement limits;
  coloured dashed curves are catalogue upper limits. Grey curves are reference
  sensitivities. Points show the posterior distributions of disk components
  in the default fit.</p>
</section>'''
    return column(
        Div(text=note, sizing_mode='stretch_width'), fig,
        sizing_mode='stretch_width',
    )


def _detection_limits(candidate, owner_key, min_wavelength=8.0):
    """Return plotting arrays for one target's detectable disk space."""

    arrays = candidate.analysis_arrays
    required = {
        'observed_flux', 'observed_uncertainty', 'observed_upper_limit',
        'observed_ignore', 'wavelength_micron', 'filter',
    }
    if not required.issubset(arrays):
        return None

    analysis = candidate.metadata.get('analysis', {})
    target = analysis.get('targets', {}).get(owner_key)
    if target is None:
        return None
    luminosity = sum(
        _summary_median(component.get('luminosity_lsun')) or 0.0
        for component in target.get('stellar_components', ())
    )
    luminosity_at_1pc = sum(
        _summary_median(component.get('luminosity_at_1pc_lsun')) or 0.0
        for component in target.get('stellar_components', ())
    )
    parallax = _summary_median(target.get('parallax_arcsec'))
    distance = (
        None if parallax is None or parallax <= 0 else 1.0 / parallax
    )
    if luminosity <= 0:
        luminosity = None
    if luminosity_at_1pc <= 0:
        luminosity_at_1pc = None
    if luminosity is None and luminosity_at_1pc is None:
        return None

    components = candidate.metadata.get('components', ())
    stellar_indices = np.asarray([
        index for index, component in enumerate(components)
        if component.get('role') == 'stellar'
    ], dtype=int)
    if not len(stellar_indices):
        return None

    selected = []
    for scope in analysis.get('data_scopes', ()):
        if (
            scope.get('kind') != 'Photometry' or
            owner_key not in scope.get('contributors', ())
        ):
            continue
        selected.extend(range(
            int(scope['element_start']), int(scope['element_stop'])
        ))
    if not selected:
        return None
    selected = np.asarray(selected, dtype=int)

    wavelength = np.asarray(arrays['wavelength_micron'][selected], dtype=float)
    filters = np.asarray(arrays['filter'][selected]).astype(str)
    observed = np.asarray(arrays['observed_flux'][selected], dtype=float)
    uncertainty = np.asarray(
        arrays['observed_uncertainty'][selected], dtype=float
    )
    upper_limit = np.asarray(
        arrays['observed_upper_limit'][selected], dtype=bool
    )
    ignored = np.asarray(arrays['observed_ignore'][selected], dtype=bool)
    stellar_uncertainty = np.zeros_like(observed)
    if {
        'component_prediction_p16', 'component_prediction_p84'
    }.issubset(arrays):
        stellar_low = np.sum(np.asarray(
            arrays['component_prediction_p16'], dtype=float
        )[stellar_indices][:, selected], axis=0)
        stellar_high = np.sum(np.asarray(
            arrays['component_prediction_p84'], dtype=float
        )[stellar_indices][:, selected], axis=0)
        stellar_uncertainty = 0.5 * np.abs(stellar_high - stellar_low)

    usable = (
        ~ignored & np.isfinite(wavelength) &
        (wavelength >= min_wavelength) & np.isfinite(observed) &
        (observed > 0) & np.isfinite(uncertainty)
    )
    usable &= np.asarray([
        not filter_module.iscolour(name) for name in filters
    ], dtype=bool)
    if not np.any(usable):
        return None
    wavelength = wavelength[usable]
    filters = filters[usable]
    observed = observed[usable]
    uncertainty = uncertainty[usable]
    upper_limit = upper_limit[usable]
    stellar_uncertainty = stellar_uncertainty[usable]
    combined_uncertainty = np.sqrt(
        uncertainty**2 + stellar_uncertainty**2
    )
    valid_uncertainty = np.isfinite(combined_uncertainty) & \
        (combined_uncertainty > 0)
    usable_limit = upper_limit | valid_uncertainty
    if not np.any(usable_limit):
        return None
    wavelength = wavelength[usable_limit]
    filters = filters[usable_limit]
    observed = observed[usable_limit]
    upper_limit = upper_limit[usable_limit]
    combined_uncertainty = combined_uncertainty[usable_limit]

    flux_limit = np.where(upper_limit, observed, 3.0 * combined_uncertainty)
    disk = BB_Disk(lstar=luminosity, distance=distance)
    kwargs = {}
    if distance is None:
        kwargs['lstar_1pc'] = luminosity_at_1pc
    curves = disk.f_limits(
        wavelength, flux_limits=flux_limit, **kwargs
    )

    reference_curves = []
    teff = _owner_teff(candidate, owner_key)
    if teff is not None and teff > 0:
        reference_curves.append((
            0.05 * disk.temperatures**3 / teff**3,
            '5% photometry', 'solid',
        ))
    try:
        alma = disk.f_limits(
            np.asarray([880.0]), flux_limits=np.asarray([20e-6]),
            fwhm=np.asarray([1.0]), **kwargs,
        )[:, 0]
    except (TypeError, ValueError):
        alma = None
    if alma is not None:
        reference_curves.append((alma, '1 mm / 1″ / 20 μJy', 'dashed'))

    return {
        'temperature_k': disk.temperatures,
        'radius_au': (
            disk.blackbody_radii() if luminosity is not None else None
        ),
        'curves': curves,
        'envelope': np.min(curves, axis=1),
        'filters': filters,
        'basis': np.where(upper_limit, 'upper limit', '3σ uncertainty'),
        'reference_curves': tuple(reference_curves),
        'disk_samples': _disk_limit_samples(candidate, owner_key),
    }


def _disk_limit_samples(candidate, owner_key, maximum=750):
    """Extract aligned derived disk samples for the limits overlay."""

    samples = candidate.analysis_arrays.get('derived_samples')
    parameters = candidate.metadata.get('analysis', {}).get(
        'derived_parameters', ()
    )
    if samples is None or not len(samples):
        return ()
    samples = np.asarray(samples, dtype=float)
    by_component = {}
    for index, parameter in enumerate(parameters):
        if (
            parameter.get('owner_key') != owner_key or
            parameter.get('role') != 'disk' or
            not parameter.get('component_id')
        ):
            continue
        name = str(parameter.get('name', '')).rsplit('.', 1)[-1]
        by_component.setdefault(parameter['component_id'], {})[name] = \
            samples[:, index]
    extracted = []
    for values in by_component.values():
        if not {'fractional_luminosity', 'temperature_k'}.issubset(values):
            continue
        count = len(values['fractional_luminosity'])
        indices = (
            np.arange(count) if count <= maximum else
            np.linspace(0, count - 1, maximum, dtype=int)
        )
        extracted.append({
            'fractional_luminosity': values['fractional_luminosity'][indices],
            'temperature_k': values['temperature_k'][indices],
            'radius_au': values.get(
                'blackbody_radius_au',
                np.repeat(np.nan, count),
            )[indices],
        })
    return tuple(extracted)


def _owner_teff(candidate, owner_key):
    for component in candidate.metadata.get('analysis', {}).get(
            'components', ()):
        if component.get('owner_key') != owner_key:
            continue
        value = _summary_median(component.get('parameters', {}).get('Teff'))
        if value is not None:
            return value
    return None


def _summary_median(summary):
    if not isinstance(summary, dict):
        return None
    value = summary.get('median')
    if value is None or not np.isfinite(value):
        return None
    return float(value)


def _candidate_panel(candidate, dataset, owner_labels, component_labels,
                     owner_colours, is_default=False):
    arrays = candidate.analysis_arrays
    required = {
        'observed_flux', 'observed_uncertainty', 'observed_upper_limit',
        'observed_ignore', 'wavelength_micron', 'filter',
        'prediction_median', 'residual', 'spectrum_wavelength_micron',
        'component_spectrum_median',
    }
    missing = required - set(arrays)
    if missing:
        raise ValueError(
            '{} needs re-analysis before plotting; missing {}'.format(
                candidate.metadata['display_name'], ', '.join(sorted(missing))
            )
        )

    scopes = candidate.metadata['analysis']['data_scopes']
    scope_labels = _scope_labels(scopes, candidate)
    scope_colours = _scope_colours(scopes, owner_colours)
    x_range, y_range = _plot_ranges(arrays)
    sed = figure(
        x_axis_type='log', y_axis_type='log',
        x_axis_label='Wavelength (μm)', y_axis_label='Flux density (Jy)',
        x_range=x_range, y_range=y_range,
        tools='pan,wheel_zoom,box_zoom,reset,save',
        active_scroll='wheel_zoom', aspect_ratio=GOLDEN_RATIO,
        sizing_mode='stretch_width',
    )
    plot_residuals = _display_residual_values(arrays)
    plot_wavelengths = np.asarray(arrays['wavelength_micron'], dtype=float)
    residual_low, residual_high = _residual_bounds(plot_residuals[
        np.isfinite(plot_wavelengths) & (plot_wavelengths > 0)
    ])
    residual = figure(
        x_axis_type='log', x_axis_location='above', x_range=sed.x_range,
        y_axis_label='Residual (σ)',
        y_range=Range1d(residual_low, residual_high),
        tools='pan,wheel_zoom,box_zoom,reset,save',
        active_scroll='wheel_zoom', height=160,
        sizing_mode='stretch_width',
    )
    for value, dash in ((0, 'solid'), (-3, 'dashed'), (3, 'dashed')):
        residual.add_layout(Span(
            location=value, dimension='width', line_color='#9ca3af',
            line_dash=dash, line_alpha=0.65,
        ))

    _add_spectra(
        sed, candidate, scopes, scope_colours, owner_colours,
    )
    for scope in scopes:
        _add_scope_data(
            sed, residual, candidate, scope,
            scope_labels[scope['observation_id']],
            scope_colours[scope['observation_id']],
        )
    _add_symbol_legend(sed)
    for plot in (sed, residual):
        plot.grid.grid_line_alpha = 0.18
        plot.outline_line_color = '#d1d5db'
    if sed.legend:
        sed.legend.location = 'top_right'
        sed.legend.orientation = 'horizontal'
        sed.legend.background_fill_alpha = 0.78
        sed.legend.border_line_alpha = 0.25
        sed.legend.label_text_font_size = '9pt'
        sed.legend.ncols = 2
        sed.legend.glyph_width = 15
        sed.legend.glyph_height = 12
        sed.legend.label_standoff = 3
        sed.legend.padding = 3
        sed.legend.margin = 4
        sed.legend.spacing = 5

    plots = gridplot(
        [[sed], [residual]], toolbar_location='below', merge_tools=True,
        sizing_mode='stretch_width',
    )

    return column(
        Div(
            text=_candidate_header(
                candidate, dataset, owner_labels, component_labels,
                owner_colours, is_default=is_default,
            ),
            sizing_mode='stretch_width', styles={'width': '100%'},
            stylesheets=[InlineStyleSheet(css='''
              .bk-clearfix { display: block !important; width: 100%; }
            ''')],
        ),
        plots,
        sizing_mode='stretch_width',
    )


def _add_spectra(fig, candidate, scopes, scope_colours, owner_colours):
    arrays = candidate.analysis_arrays
    wavelength = arrays['spectrum_wavelength_micron']
    components = arrays['component_spectrum_median']
    metadata = candidate.metadata['components']

    for index, component in enumerate(metadata):
        owner = component['owner_key']
        fig.line(
            wavelength, components[index],
            color=owner_colours.get(owner, '#6b7280'),
            line_width=1.5, line_alpha=0.6,
            line_dash='dashed' if component['role'] == 'stellar' else 'dotted',
            legend_label='star' if component['role'] == 'stellar' else 'disk',
        )

    seen = set()
    for scope in scopes:
        observation_id = scope['observation_id']
        if observation_id in seen:
            continue
        seen.add(observation_id)
        contributors = set(scope['contributors'])
        visible = np.array([
            component['owner_key'] in contributors for component in metadata
        ])
        combined = np.sum(components[visible], axis=0)
        fig.line(
            wavelength, combined,
            color=scope_colours[observation_id], line_width=3,
            line_alpha=0.82,
            legend_label='total',
        )


def _add_scope_data(sed, residual, candidate, scope, label, colour):
    arrays = candidate.analysis_arrays
    start = scope['element_start']
    stop = scope['element_stop']
    indices = np.arange(start, stop)
    wavelength = np.asarray(arrays['wavelength_micron'][indices], dtype=float)
    observed = np.asarray(arrays['observed_flux'][indices], dtype=float)
    uncertainty = np.asarray(
        arrays['observed_uncertainty'][indices], dtype=float
    )
    predicted = np.asarray(arrays['prediction_median'][indices], dtype=float)
    residual_values = _display_residual_values(arrays)[indices]
    upper_limit = np.asarray(
        arrays['observed_upper_limit'][indices], dtype=bool
    )
    ignored = np.asarray(arrays['observed_ignore'][indices], dtype=bool)
    filters = np.asarray(arrays['filter'][indices]).astype(str)
    usable_x = np.isfinite(wavelength) & (wavelength > 0)

    if scope['kind'] == 'ObsSpectrum':
        _add_spectral_scope_data(
            sed, residual, wavelength, observed, uncertainty, predicted,
            residual_values, ignored, usable_x, label, colour,
        )
        return

    normal = usable_x & np.isfinite(observed) & (observed > 0) & ~upper_limit
    active = normal & ~ignored
    if np.any(active):
        source = _point_source(
            wavelength[active], observed[active], uncertainty[active],
            predicted[active], residual_values[active], filters[active], label,
        )
        renderer = sed.scatter(
            'wavelength', 'observed', source=source, marker='circle',
            size=DATA_MARKER_SIZE,
            color=colour, line_color='white', line_width=0.7,
        )
        _add_hover(sed, renderer)
        error_ok = source.data['lower'] > 0
        if np.any(error_ok):
            error_source = ColumnDataSource({
                key: np.asarray(value)[error_ok]
                for key, value in source.data.items()
            })
            sed.segment(
                'wavelength', 'lower', 'wavelength', 'upper',
                source=error_source, color=colour, line_alpha=0.65,
            )

    excluded = normal & ignored
    if np.any(excluded):
        source = _point_source(
            wavelength[excluded], observed[excluded], uncertainty[excluded],
            predicted[excluded], residual_values[excluded], filters[excluded],
            label,
        )
        renderer = sed.scatter(
            'wavelength', 'observed', source=source, marker='circle',
            size=DATA_MARKER_SIZE, color=EXCLUDED_COLOUR,
            line_color='white', line_width=0.7, alpha=0.7,
        )
        _add_hover(sed, renderer)

    limits = usable_x & np.isfinite(observed) & (observed > 0) & upper_limit
    for limit_mask, limit_colour, limit_alpha in (
            (limits & ~ignored, colour, 0.9),
            (limits & ignored, EXCLUDED_COLOUR, 0.7)):
        if not np.any(limit_mask):
            continue
        source = _point_source(
            wavelength[limit_mask], observed[limit_mask],
            uncertainty[limit_mask], predicted[limit_mask],
            residual_values[limit_mask], filters[limit_mask], label,
        )
        renderer = sed.scatter(
            'wavelength', 'observed', source=source,
            marker='inverted_triangle', size=UPPER_LIMIT_MARKER_SIZE,
            color=limit_colour, alpha=limit_alpha,
        )
        _add_hover(sed, renderer)

    predicted_ok = usable_x & np.isfinite(predicted) & (predicted > 0)
    if np.any(predicted_ok):
        source = _point_source(
            wavelength[predicted_ok], observed[predicted_ok],
            uncertainty[predicted_ok], predicted[predicted_ok],
            residual_values[predicted_ok], filters[predicted_ok], label,
        )
        renderer = sed.scatter(
            'wavelength', 'predicted', source=source,
            marker='circle', size=MODEL_MARKER_SIZE, fill_alpha=0,
            line_color=colour, line_width=1.4,
        )
        _add_hover(sed, renderer)

    residual_groups = (
        (normal & ~ignored, 'circle', DATA_MARKER_SIZE, colour, 1.0),
        (normal & ignored, 'circle', DATA_MARKER_SIZE,
         EXCLUDED_COLOUR, 0.7),
        (limits & ~ignored, 'inverted_triangle', UPPER_LIMIT_MARKER_SIZE,
         colour, 0.9),
        (limits & ignored, 'inverted_triangle', UPPER_LIMIT_MARKER_SIZE,
         EXCLUDED_COLOUR, 0.7),
    )
    for mask, marker, size, marker_colour, alpha in residual_groups:
        mask = mask & np.isfinite(residual_values)
        if not np.any(mask):
            continue
        source = _point_source(
            wavelength[mask], observed[mask], uncertainty[mask],
            predicted[mask], residual_values[mask], filters[mask], label,
        )
        renderer = residual.scatter(
            'wavelength', 'residual', source=source, marker=marker, size=size,
            color=marker_colour, alpha=alpha,
        )
        _add_hover(residual, renderer)


def _add_symbol_legend(fig):
    """Add a stable, compact key for photometry point encodings."""

    # A real (but non-plottable) row is needed for Bokeh to draw the legend
    # glyph. An empty source leaves the label visible but omits its symbol.
    source = ColumnDataSource({'wavelength': [np.nan], 'value': [np.nan]})
    fig.scatter(
        'wavelength', 'value', source=source, marker='circle',
        size=DATA_MARKER_SIZE, color=LEGEND_COLOUR,
        line_color='white', line_width=0.7, legend_label='data',
    )
    fig.scatter(
        'wavelength', 'value', source=source, marker='inverted_triangle',
        size=UPPER_LIMIT_MARKER_SIZE, color=LEGEND_COLOUR,
        legend_label='upper limit',
    )
    fig.scatter(
        'wavelength', 'value', source=source, marker='circle',
        size=DATA_MARKER_SIZE, color=EXCLUDED_COLOUR,
        line_color='white', line_width=0.7, alpha=0.7,
        legend_label='excluded',
    )
    fig.scatter(
        'wavelength', 'value', source=source, marker='circle',
        size=MODEL_MARKER_SIZE, fill_alpha=0,
        line_color=LEGEND_COLOUR, line_width=1.4, legend_label='model',
    )


def _add_spectral_scope_data(
        sed, residual, wavelength, observed, uncertainty, predicted,
        residual_values, ignored, usable_x, label, colour):
    """Draw a spectral segment as traces rather than photometric points."""

    usable = (
        usable_x & np.isfinite(observed) & (observed > 0) &
        np.isfinite(predicted) & (predicted > 0) & ~ignored
    )
    if not np.any(usable):
        return
    source = _point_source(
        wavelength[usable], observed[usable], uncertainty[usable],
        predicted[usable], residual_values[usable],
        np.repeat('spectrum', np.sum(usable)), label,
    )
    observed_renderer = sed.line(
        'wavelength', 'observed', source=source,
        color=colour, line_width=1.5, line_alpha=0.78,
        legend_label='spectrum',
    )
    predicted_renderer = sed.line(
        'wavelength', 'predicted', source=source,
        color=colour, line_width=2.2, line_alpha=0.95,
        line_dash='dashed',
    )
    residual_renderer = residual.line(
        'wavelength', 'residual', source=source,
        color=colour, line_width=1.35, line_alpha=0.85,
    )
    for plot, renderer in (
        (sed, observed_renderer), (sed, predicted_renderer),
        (residual, residual_renderer),
    ):
        _add_hover(plot, renderer)


def _point_source(wavelength, observed, uncertainty, predicted, residual,
                  filters, scope):
    lower = observed - uncertainty
    upper = observed + uncertainty
    return ColumnDataSource({
        'wavelength': wavelength,
        'observed': observed,
        'uncertainty': uncertainty,
        'lower': lower,
        'upper': upper,
        'predicted': predicted,
        'residual': residual,
        'filter': filters,
        'scope': np.repeat(scope, len(wavelength)),
    })


def _add_hover(fig, renderer):
    fig.add_tools(HoverTool(
        renderers=[renderer],
        tooltips=[
            ('scope', '@scope'), ('band', '@filter'),
            ('wavelength', '@wavelength{0.000} μm'),
            ('observed', '@observed{0.0000} Jy'),
            ('predicted', '@predicted{0.0000} Jy'),
            ('residual', '@residual{0.00} σ'),
        ],
    ))


def _candidate_header(candidate, dataset, owner_labels, component_labels,
                      owner_colours, is_default=False):
    evidence = candidate.posterior.evidence
    evidence_error = candidate.posterior.evidence_error
    evidence_text = 'not available' if evidence is None else '{:.2f}'.format(evidence)
    if evidence_error is not None:
        evidence_text += ' ± {:.2f}'.format(evidence_error)
    analysis = candidate.metadata['analysis']
    rows = []
    analysis_components = {
        component['component_id']: component
        for component in analysis.get('components', ())
    }
    target_summaries = analysis.get('targets', {})
    for owner, label in owner_labels.items():
        components = [
            component for component in candidate.metadata['components']
            if component['owner_key'] == owner
        ]
        if not components:
            continue
        summary = target_summaries.get(owner, {})
        derived = {
            component['component_id']: component
            for role in ('stellar_components', 'disk_components')
            for component in summary.get(role, ())
        }
        details = []
        disk_index = 0
        disk_count = sum(
            component['role'] == 'disk' for component in components
        )
        for component in components:
            if component['role'] == 'stellar':
                prefix = '★'
                values = _stellar_summary(
                    analysis_components.get(component['component_id'], {}),
                    derived.get(component['component_id'], {}),
                )
            else:
                disk_index += 1
                prefix = '○{}'.format(disk_index) if disk_count > 1 else '○'
                values = _disk_summary(
                    analysis_components.get(component['component_id'], {}),
                    derived.get(component['component_id'], {})
                )
            details.append(
                '<div style="margin:5px 0 0 14px">'
                '<div><strong>{}</strong> '
                '<span style="color:#475569">{}</span></div>'
                '<div style="margin:2px 0 0 19px;color:#334155">{}</div>'
                '</div>'.format(
                    prefix, escape(component['model_name']),
                    ' · '.join(
                        '<span style="white-space:nowrap">{}</span>'.format(
                            value
                        )
                        for value in values
                    ) if values else '—',
                )
            )
        short = component_labels.get(owner) or label
        symbols = '★' * sum(
            component['role'] == 'stellar' for component in components
        ) + '○' * sum(component['role'] == 'disk' for component in components)
        colour = owner_colours.get(owner, '#64748b')
        rows.append(
            '<div style="border-left:4px solid {};padding:4px 9px;min-width:0">'
            '<div><strong>{} {}</strong>'
            '<span style="color:#64748b"> · {}</span></div>'
            '{}</div>'.format(
                colour, escape(str(short)), symbols, escape(label),
                ''.join(details),
            )
        )
    badge = (
        ' <span style="background:#dff3e4;color:#176b37;border-radius:10px;'
        'padding:2px 7px;font-size:11px">default</span>'
        if is_default else ''
    )
    stale_warning = '' if candidate.is_current(dataset) else (
        ' · <strong style="color:#a13b2b">fit input has changed</strong>'
    )
    product_links = _candidate_product_links(candidate)
    return '''
<section style="padding:8px 10px 6px;margin:0;width:100%;max-width:none;
font:13px ui-sans-serif,system-ui,sans-serif">
  <div style="display:flex;align-items:baseline;gap:12px;flex-wrap:wrap">
    <strong style="font-size:15px">Fit results{}</strong>
    <span style="margin-left:auto;text-align:right">{}</span>
  </div>
  <div style="color:#526071;margin-top:2px">
    ln Z {} · χ² {:.1f} · {} residuals{}
  </div>
  <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));
gap:7px;margin-top:7px">{}</div>
  <div style="color:#64748b;margin-top:7px;font-size:12px">
  ★ stellar component · ○ disk component</div>
</section>
'''.format(
        badge, product_links, evidence_text, analysis['chi_squared'],
        analysis['n_residuals'], stale_warning, ''.join(rows),
    )


def _candidate_product_links(candidate):
    directory = 'fits/{}/'.format(candidate.path.name)
    links = [
        '<a style="color:#285ea8;text-decoration:none" href="{}">{}</a>'.format(
            escape(directory + 'result.json'), 'metadata (JSON)'
        ),
        '<a style="color:#285ea8;text-decoration:none" href="{}">{}</a>'.format(
            escape(directory + candidate.metadata['posterior']['file']),
            'posterior data (NPZ)',
        ),
    ]
    for label, title, filename in (
        ('parameters (PNG)', 'Fitted parameter posteriors', 'corner.png'),
        ('derived (PNG)', 'Derived-quantity distributions', 'distributions.png'),
    ):
        if (candidate.path / filename).is_file():
            links.append(
                '<a style="color:#285ea8;text-decoration:none" href="{}" '
                'data-distribution-plot data-plot-title="{}">{}</a>'.format(
                    escape(directory + filename), escape(title), escape(label)
                )
            )
    return ' · '.join(links)


def _stellar_summary(component, derived):
    values = []
    parameters = component.get('parameters', {})
    if 'Teff' in parameters:
        values.append('T<sub>eff</sub> {}'.format(
            _format_summary(parameters['Teff'], 'K', digits=0)
        ))
    if 'radius_rsun' in derived:
        values.append('R {}'.format(
            _format_summary(derived['radius_rsun'], 'R<sub>☉</sub>')
        ))
    if 'luminosity_lsun' in derived:
        values.append('L<sub>★</sub> {}'.format(
            _format_summary(
                derived['luminosity_lsun'], 'L<sub>☉</sub>', scientific=True,
            )
        ))
    return values


def _disk_summary(component, derived):
    values = []
    parameters = component.get('parameters', {})
    if 'temperature_k' in derived:
        values.append('T {}'.format(
            _format_summary(derived['temperature_k'], 'K', digits=0)
        ))
    if 'beta' in parameters:
        values.append('β {}'.format(
            _format_summary(parameters['beta'], '')
        ))
    if 'log_lam0' in parameters:
        values.append('λ<sub>0</sub> {}'.format(
            _format_summary(
                _pow10_summary(parameters['log_lam0']), 'μm'
            )
        ))
    elif 'lam0' in parameters:
        values.append('λ<sub>0</sub> {}'.format(
            _format_summary(parameters['lam0'], 'μm')
        ))
    if 'blackbody_radius_au' in derived:
        values.append('R<sub>BB</sub> {}'.format(
            _format_summary(derived['blackbody_radius_au'], 'au')
        ))
    if 'fractional_luminosity' in derived:
        values.append('L<sub>d</sub>/L<sub>★</sub> {}'.format(
            _format_summary(
                derived['fractional_luminosity'], '', scientific=True,
            )
        ))
    return values


def _pow10_summary(summary):
    """Transform a base-10 logarithmic quantile summary to linear units."""

    median = summary.get('median')
    minus = summary.get('minus')
    plus = summary.get('plus')
    if median is None or not np.isfinite(median):
        return {'median': None}
    linear_median = 10**median
    transformed = {'median': linear_median}
    if (minus is not None and plus is not None and
            np.isfinite(minus) and np.isfinite(plus)):
        transformed.update({
            'minus': linear_median - 10**(median - minus),
            'plus': 10**(median + plus) - linear_median,
        })
    return transformed


def _format_summary(summary, unit, digits=2, scientific=False):
    median = summary.get('median')
    if median is None or not np.isfinite(median):
        return '—'
    minus = summary.get('minus')
    plus = summary.get('plus')
    uncertainty = None
    if (minus is not None and plus is not None and
            np.isfinite(minus) and np.isfinite(plus)):
        uncertainty = (minus + plus) / 2
    if scientific:
        number = '{:.2e}'.format(median)
        if uncertainty is not None:
            number += ' ± {:.2e}'.format(uncertainty)
    elif abs(median) >= 1e4 or (median != 0 and abs(median) < 1e-2):
        number = '{:.2g}'.format(median)
    else:
        number = '{:.{}f}'.format(median, digits)
        if uncertainty is not None:
            number += ' ± {:.{}f}'.format(uncertainty, digits)
    return '{}{}'.format(number, ' ' + unit if unit else '')


def _owner_labels(candidate):
    targets = candidate.metadata.get('analysis', {}).get('targets', {})
    return {
        key: values.get('main_id') or values.get('sdbid') or key
        for key, values in targets.items()
    }


def _owner_colour_map(owner_labels):
    return {
        owner: Category10[10][index % 10]
        for index, owner in enumerate(owner_labels)
    }


def _candidate_shorthand_text(candidate, owner_labels, component_labels):
    grouped = {}
    order = []
    for component in candidate.metadata['components']:
        owner = component['owner_key']
        if owner not in grouped:
            grouped[owner] = {'stellar': 0, 'disk': 0}
            order.append(owner)
        grouped[owner][component['role']] += 1
    parts = []
    for owner in order:
        values = grouped[owner]
        label = component_labels.get(owner)
        if label is None and len(order) > 1:
            label = owner_labels.get(owner, owner)
        symbols = '★' * values['stellar'] + '○' * values['disk']
        parts.append('{} {}'.format(label, symbols) if label else symbols)
    return ' · '.join(parts)


def _system_title(systems, owner_labels, package_directory):
    if len(systems) == 1 and systems[0].get('label'):
        return systems[0]['label']
    labels = [system.get('label') for system in systems if system.get('label')]
    return ' + '.join(labels or owner_labels.values()) or package_directory.name


def _candidate_label(candidate, owner_labels):
    grouped = {}
    for component in candidate.metadata['components']:
        grouped.setdefault(component['owner_key'], []).append(
            component['model_name']
        )
    return ' · '.join(
        '{}: {}'.format(owner_labels.get(owner, owner), ' + '.join(models))
        for owner, models in grouped.items()
    )


def _scope_labels(scopes, candidate):
    targets = candidate.metadata['analysis'].get('targets', {})
    by_sdbid = {
        values.get('sdbid'): values.get('main_id') or values.get('sdbid')
        for values in targets.values()
    }
    labels = {}
    for scope in scopes:
        observation_id = scope['observation_id']
        if observation_id in labels:
            continue
        if observation_id in by_sdbid:
            labels[observation_id] = by_sdbid[observation_id]
        else:
            labels[observation_id] = ' + '.join(
                targets.get(key, {}).get('main_id') or key
                for key in scope['contributors']
            )
        if scope.get('segment_id'):
            labels[observation_id] += ' · {}'.format(scope['segment_id'])
    return labels


def _scope_colours(scopes, owner_colours):
    """Keep resolved observations tied to their physical contributor colour."""

    composite_colours = ('#334155', '#7c3aed', '#a16207', '#0f766e')
    colours = {}
    composite_index = 0
    for scope in scopes:
        observation_id = scope['observation_id']
        if observation_id in colours:
            continue
        contributors = tuple(scope['contributors'])
        if len(contributors) == 1 and contributors[0] in owner_colours:
            colours[observation_id] = owner_colours[contributors[0]]
        else:
            colours[observation_id] = composite_colours[
                composite_index % len(composite_colours)
            ]
            composite_index += 1
    return colours


def _plot_ranges(arrays):
    wavelength = np.asarray(arrays['wavelength_micron'], dtype=float)
    observed = np.asarray(arrays['observed_flux'], dtype=float)
    predicted = np.asarray(arrays['prediction_median'], dtype=float)
    valid_x = wavelength[np.isfinite(wavelength) & (wavelength > 0)]
    x_low = np.min(valid_x) / 1.6
    x_high = np.max(valid_x) * 1.6
    fluxes = np.concatenate((
        observed[np.isfinite(observed) & (observed > 0)],
        predicted[np.isfinite(predicted) & (predicted > 0)],
    ))
    fluxes = fluxes[np.isfinite(fluxes) & (fluxes > 0)]
    return (
        Range1d(x_low, x_high),
        Range1d(np.min(fluxes) / 2.5, np.max(fluxes) * 2.5),
    )


def _display_residual_values(arrays):
    """Return residuals suitable for display, including ignored data."""

    values = np.asarray(arrays['residual'], dtype=float).copy()
    observed = np.asarray(arrays['observed_flux'], dtype=float)
    uncertainty = np.asarray(arrays['observed_uncertainty'], dtype=float)
    predicted = np.asarray(arrays['prediction_median'], dtype=float)
    upper_limit = np.asarray(arrays['observed_upper_limit'], dtype=bool)
    ignored = np.asarray(arrays['observed_ignore'], dtype=bool)

    ignored_data = ignored & ~upper_limit
    values[ignored_data] = np.nan
    valid_data = (
        ignored_data & np.isfinite(observed) & np.isfinite(predicted) &
        np.isfinite(uncertainty) & (uncertainty > 0)
    )
    values[valid_data] = (
        observed[valid_data] - predicted[valid_data]
    ) / uncertainty[valid_data]

    ignored_limits = (
        ignored & upper_limit & np.isfinite(observed) & (observed > 0) &
        np.isfinite(predicted)
    )
    values[ignored & upper_limit] = np.nan
    threshold = observed[ignored_limits] / 3.0
    values[ignored_limits] = np.where(
        predicted[ignored_limits] > threshold,
        -predicted[ignored_limits] / threshold,
        0.0,
    )
    return values


def _residual_bounds(values):
    """Autoscale residuals without clipping, retaining at least -4 to 4."""

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return -4.0, 4.0
    low_value = float(np.min(values))
    high_value = float(np.max(values))
    padding = max(0.25, (high_value - low_value) * 0.05)
    low = low_value - padding if low_value < -4.0 else -4.0
    high = high_value + padding if high_value > 4.0 else 4.0
    return low, high


def _sortable_evidence(candidate):
    return candidate.posterior.evidence \
        if candidate.posterior.evidence is not None else -np.inf
