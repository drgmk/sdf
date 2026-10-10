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
from bokeh.palettes import Category10
from bokeh.plotting import figure

from .candidate_library import CandidateLibrary
from .dataset import FitDataset, canonical_fit_source
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


def write_joint_results_page(package_directory, output=None, *, dataset=None,
                             inline=True, candidate_ids=None, operator=False,
                             policy_store=None, operator_api=None,
                             operator_token=None, operator_actor=None,
                             operator_actor_locked=False):
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
    owner_labels = _owner_labels(ordered[0])
    owner_colours = _owner_colour_map(owner_labels)
    target_rows = target_page_rows(dataset, package_directory)
    for target in target_rows:
        target['colour'] = owner_colours.get(target['key'], '#64748b')
    observation_rows = observation_page_rows(dataset, package_directory)
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
    )
    default_name = 'operator.html' if operator else 'index.html'
    output = (
        Path(output) if output is not None else package_directory / default_name
    )
    output.write_text(html, encoding='utf-8')
    return output


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
    residual_limit = _residual_limit(arrays['residual'])
    residual = figure(
        x_axis_type='log', x_axis_location='above', x_range=sed.x_range,
        y_axis_label='Residual (σ)',
        y_range=Range1d(-residual_limit, residual_limit),
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
    for plot in (sed, residual):
        plot.grid.grid_line_alpha = 0.18
        plot.outline_line_color = '#d1d5db'
    if sed.legend:
        sed.legend.location = 'top_right'
        sed.legend.orientation = 'horizontal'
        sed.legend.background_fill_alpha = 0.78
        sed.legend.border_line_alpha = 0.25
        sed.legend.label_text_font_size = '9pt'
        sed.legend.glyph_width = 18
        sed.legend.padding = 4
        sed.legend.spacing = 8

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
    residual_values = np.asarray(arrays['residual'][indices], dtype=float)
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
            'wavelength', 'observed', source=source, marker='circle', size=12,
            color=colour, line_color='white', line_width=0.7,
            legend_label='data',
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
            'wavelength', 'observed', source=source, marker='x', size=13.5,
            color=colour, alpha=0.3,
        )
        _add_hover(sed, renderer)

    limits = usable_x & np.isfinite(observed) & (observed > 0) & upper_limit
    if np.any(limits):
        source = _point_source(
            wavelength[limits], observed[limits], uncertainty[limits],
            predicted[limits], residual_values[limits], filters[limits], label,
        )
        source.data['alpha'] = np.where(ignored[limits], 0.3, 0.9)
        renderer = sed.scatter(
            'wavelength', 'observed', source=source,
            marker='inverted_triangle', size=10, color=colour,
            alpha='alpha',
            legend_label='upper lim.',
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
            marker='circle', size=10.5, fill_alpha=0,
            line_color=colour, line_width=1.4,
            legend_label='model',
        )
        _add_hover(sed, renderer)

    residual_ok = usable_x & np.isfinite(residual_values) & ~ignored
    if np.any(residual_ok):
        source = _point_source(
            wavelength[residual_ok], observed[residual_ok],
            uncertainty[residual_ok], predicted[residual_ok],
            residual_values[residual_ok], filters[residual_ok], label,
        )
        renderer = residual.scatter(
            'wavelength', 'residual', source=source, marker='circle', size=10.5,
            color=colour,
        )
        _add_hover(residual, renderer)


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
                    derived.get(component['component_id'], {})
                )
            details.append(
                '<div style="margin:5px 0 0 14px">'
                '<div><strong>{}</strong> '
                '<span style="color:#475569">{}</span></div>'
                '<div style="margin:2px 0 0 19px;color:#334155">{}</div>'
                '</div>'.format(
                    prefix, escape(component['model_name']),
                    ' · '.join(values) if values else '—',
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
    status = 'current input' if candidate.is_current(dataset) else 'input changed'
    product_links = _candidate_product_links(candidate)
    return '''
<section style="padding:8px 10px 6px;margin:0;width:100%;max-width:none;
font:13px ui-sans-serif,system-ui,sans-serif">
  <div style="display:flex;justify-content:space-between;gap:10px;flex-wrap:wrap">
    <strong style="font-size:15px">Fit results{}</strong>
    <span style="color:#526071">ln Z {} · χ² {:.1f} · {} residuals · {} · {}</span>
  </div>
  <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));
gap:7px;margin-top:7px">{}</div>
  <div style="color:#64748b;margin-top:7px;font-size:12px">
  ★ stellar component · ○ disk component</div>
</section>
'''.format(
        badge, evidence_text, analysis['chi_squared'],
        analysis['n_residuals'], status, product_links, ''.join(rows),
    )


def _candidate_product_links(candidate):
    directory = 'fits/{}/'.format(candidate.path.name)
    links = [
        ('result', directory + 'result.json'),
        ('posterior', directory + candidate.metadata['posterior']['file']),
    ]
    for label, filename in (
        ('fit', 'corner.png'), ('derived', 'distributions.png'),
    ):
        if (candidate.path / filename).is_file():
            links.append((label, directory + filename))
    return ' '.join(
        '<a style="color:#285ea8;text-decoration:none" href="{}">{}</a>'.format(
            escape(url), escape(label)
        ) for label, url in links
    )


def _stellar_summary(component, derived):
    values = []
    parameters = component.get('parameters', {})
    if 'Teff' in parameters:
        values.append('T<sub>eff</sub> {}'.format(
            _format_summary(parameters['Teff'], 'K', digits=0)
        ))
    if 'logg' in parameters:
        values.append('log(g) {}'.format(
            _format_summary(parameters['logg'], '')
        ))
    if 'radius_rsun' in derived:
        values.append('R {}'.format(
            _format_summary(derived['radius_rsun'], 'R<sub>☉</sub>')
        ))
    return values


def _disk_summary(derived):
    values = []
    if 'temperature_k' in derived:
        values.append('T {}'.format(
            _format_summary(derived['temperature_k'], 'K', digits=0)
        ))
    if 'blackbody_radius_au' in derived:
        values.append('R<sub>BB</sub> {}'.format(
            _format_summary(derived['blackbody_radius_au'], 'au')
        ))
    if 'fractional_luminosity' in derived:
        values.append('L<sub>d</sub>/L<sub>★</sub> {}'.format(
            _format_summary(derived['fractional_luminosity'], '')
        ))
    return values


def _format_summary(summary, unit, digits=2):
    median = summary.get('median')
    if median is None or not np.isfinite(median):
        return '—'
    minus = summary.get('minus')
    plus = summary.get('plus')
    if abs(median) >= 1e4 or (median != 0 and abs(median) < 1e-2):
        number = '{:.2g}'.format(median)
        uncertainty = None
    else:
        number = '{:.{}f}'.format(median, digits)
        uncertainty = None
        if (minus is not None and plus is not None and
                np.isfinite(minus) and np.isfinite(plus)):
            uncertainty = (minus + plus) / 2
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


def _residual_limit(values):
    values = np.abs(np.asarray(values, dtype=float))
    values = values[np.isfinite(values)]
    if not len(values):
        return 5.0
    return max(5.0, min(50.0, float(np.percentile(values, 95)) * 1.15))


def _sortable_evidence(candidate):
    return candidate.posterior.evidence \
        if candidate.posterior.evidence is not None else -np.inf
