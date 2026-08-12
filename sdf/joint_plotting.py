"""Bokeh presentation for compact joint-fit candidate artifacts."""

from __future__ import annotations

from html import escape
from pathlib import Path

import numpy as np
from bokeh.embed import file_html
from bokeh.layouts import column
from bokeh.models import (
    ColumnDataSource, Div, HoverTool, Range1d, Span, TabPanel, Tabs,
)
from bokeh.palettes import Category10
from bokeh.plotting import figure
from bokeh.resources import CDN, INLINE

from .candidate_library import CandidateLibrary
from .dataset import FitDataset, canonical_fit_source


def write_joint_results_page(package_directory, output=None, *, dataset=None,
                             inline=True):
    """Write a standalone candidate-review page for one joint package."""

    package_directory = canonical_fit_source(package_directory)
    if not package_directory.is_dir():
        package_directory = package_directory.parent
    if dataset is None:
        dataset = FitDataset.from_package(package_directory)
    library = CandidateLibrary(package_directory)
    candidates = library.candidates()
    if not candidates:
        raise ValueError('joint package has no retained candidates to plot')

    try:
        selection = library.read_selection()
        order = list(selection.displayed)
        order.extend(
            candidate.candidate_id for candidate in candidates
            if candidate.candidate_id not in order
        )
        default = selection.default
    except FileNotFoundError:
        order = [
            candidate.candidate_id for candidate in sorted(
                candidates,
                key=lambda value: _sortable_evidence(value),
                reverse=True,
            )
        ]
        default = order[0]

    by_id = {candidate.candidate_id: candidate for candidate in candidates}
    ordered = [by_id[candidate_id] for candidate_id in order]
    owner_labels = _owner_labels(ordered[0])
    title = ' + '.join(owner_labels.values()) or package_directory.name

    panels = []
    for rank, candidate in enumerate(ordered, start=1):
        label = _candidate_label(candidate, owner_labels)
        panels.append(TabPanel(
            child=_candidate_panel(
                candidate, dataset, owner_labels,
                is_default=candidate.candidate_id == default,
            ),
            title='Fit {}'.format(rank),
        ))
    active = next(
        (index for index, candidate in enumerate(ordered)
         if candidate.candidate_id == default),
        0,
    )
    tabs = Tabs(tabs=panels, active=active, sizing_mode='stretch_width')
    page = column(
        Div(text=_page_header(title, ordered, dataset, default),
            sizing_mode='stretch_width'),
        tabs,
        sizing_mode='stretch_width',
    )
    html = file_html(
        page, INLINE if inline else CDN,
        title='{} — SDF joint fits'.format(title),
    )
    output = Path(output) if output is not None else package_directory / 'index.html'
    output.write_text(html, encoding='utf-8')
    return output


def _candidate_panel(candidate, dataset, owner_labels, is_default=False):
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
    scope_colours = {
        observation_id: Category10[10][index % 10]
        for index, observation_id in enumerate(scope_labels)
    }
    owner_colours = {
        owner: Category10[10][index % 10]
        for index, owner in enumerate(owner_labels)
    }
    x_range, y_range = _plot_ranges(arrays)
    sed = figure(
        x_axis_type='log', y_axis_type='log',
        x_axis_label='Wavelength (μm)', y_axis_label='Flux density (Jy)',
        x_range=x_range, y_range=y_range,
        tools='pan,wheel_zoom,box_zoom,reset,save',
        active_scroll='wheel_zoom', height=560,
        sizing_mode='stretch_width', toolbar_location='above',
    )
    residual_limit = _residual_limit(arrays['residual'])
    residual = figure(
        x_axis_type='log', x_range=sed.x_range,
        y_axis_label='Residual (σ)', x_axis_label='Wavelength (μm)',
        y_range=Range1d(-residual_limit, residual_limit),
        tools='pan,wheel_zoom,box_zoom,reset,save',
        active_scroll='wheel_zoom', height=260,
        sizing_mode='stretch_width', toolbar_location='above',
    )
    for value, dash in ((0, 'solid'), (-3, 'dashed'), (3, 'dashed')):
        residual.add_layout(Span(
            location=value, dimension='width', line_color='#9ca3af',
            line_dash=dash, line_alpha=0.65,
        ))

    _add_spectra(
        sed, candidate, scopes, scope_labels, scope_colours,
        owner_labels, owner_colours,
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
        if plot.legend:
            plot.legend.click_policy = 'hide'
            plot.legend.label_text_font_size = '9pt'
            plot.add_layout(plot.legend[0], 'right')

    return column(
        Div(text=_candidate_header(
            candidate, dataset, owner_labels, is_default=is_default
        ), sizing_mode='stretch_width'),
        sed,
        residual,
        sizing_mode='stretch_width',
    )


def _add_spectra(fig, candidate, scopes, scope_labels, scope_colours,
                 owner_labels, owner_colours):
    arrays = candidate.analysis_arrays
    wavelength = arrays['spectrum_wavelength_micron']
    components = arrays['component_spectrum_median']
    metadata = candidate.metadata['components']

    for index, component in enumerate(metadata):
        owner = component['owner_key']
        label = '{} {} {}'.format(
            owner_labels.get(owner, owner),
            component['model_name'],
            component['role'],
        )
        fig.line(
            wavelength, components[index],
            color=owner_colours.get(owner, '#6b7280'),
            line_width=1.5, line_alpha=0.6,
            line_dash='dashed' if component['role'] == 'stellar' else 'dotted',
            legend_label=label,
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
            legend_label='{} model'.format(scope_labels[observation_id]),
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

    normal = usable_x & np.isfinite(observed) & (observed > 0) & ~upper_limit
    active = normal & ~ignored
    if np.any(active):
        source = _point_source(
            wavelength[active], observed[active], uncertainty[active],
            predicted[active], residual_values[active], filters[active], label,
        )
        renderer = sed.scatter(
            'wavelength', 'observed', source=source, marker='circle', size=8,
            color=colour, line_color='white', line_width=0.7,
            legend_label='{} data'.format(label),
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
            'wavelength', 'observed', source=source, marker='x', size=9,
            color=colour, alpha=0.3,
            legend_label='{} ignored'.format(label),
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
            legend_label='{} upper limit'.format(label),
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
            marker='diamond', size=7, fill_alpha=0,
            line_color=colour, line_width=1.4,
            legend_label='{} predicted'.format(label),
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
            'wavelength', 'residual', source=source, marker='circle', size=7,
            color=colour, legend_label=label,
        )
        _add_hover(residual, renderer)


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


def _page_header(title, candidates, dataset, default):
    rows = []
    labels = _owner_labels(candidates[0])
    for index, candidate in enumerate(candidates, start=1):
        evidence = candidate.posterior.evidence
        evidence_text = '—' if evidence is None else '{:.2f}'.format(evidence)
        chi_squared = candidate.metadata['analysis'].get('chi_squared')
        chi_text = '—' if chi_squared is None else '{:.1f}'.format(chi_squared)
        status = 'current' if candidate.is_current(dataset) else 'input changed'
        default_text = ' <span class="default-badge">default</span>' \
            if candidate.candidate_id == default else ''
        rows.append(
            '<tr><td>{}</td><td>{}{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.
            format(
                index, escape(_candidate_label(candidate, labels)),
                default_text, evidence_text, chi_text, status,
            )
        )
    return '''
<style>
  :root {{ color-scheme: light; }}
  body {{ margin: 0; background: #f5f7fa; color: #172033;
          font-family: ui-sans-serif, system-ui, -apple-system, sans-serif; }}
  .joint-head, .candidate-head {{ background: white; border: 1px solid #dfe5ec;
          border-radius: 8px; padding: 14px 18px; margin: 10px 8px; }}
  .joint-head h1 {{ margin: 0 0 4px; font-size: 24px; }}
  .joint-head p, .candidate-head p {{ margin: 4px 0; color: #526071; }}
  .fit-table, .component-table {{ border-collapse: collapse; width: 100%;
          margin-top: 12px; font-size: 13px; }}
  .fit-table th, .fit-table td, .component-table th, .component-table td {{
          border-bottom: 1px solid #e6ebf0; text-align: left; padding: 6px 8px; }}
  .fit-table th, .component-table th {{ color: #526071; font-weight: 600; }}
  .default-badge {{ background: #dff3e4; color: #176b37; border-radius: 10px;
          padding: 2px 7px; font-size: 11px; white-space: nowrap; }}
</style>
<section class="joint-head">
  <h1>{}</h1>
  <p>Retained joint fits. Select a candidate tab to compare its scoped SED and residuals.</p>
  <table class="fit-table">
    <thead><tr><th>Tab</th><th>Candidate</th><th>Evidence</th><th>χ²</th><th>Input</th></tr></thead>
    <tbody>{}</tbody>
  </table>
</section>
'''.format(escape(title), ''.join(rows))


def _candidate_header(candidate, dataset, owner_labels, is_default=False):
    evidence = candidate.posterior.evidence
    evidence_error = candidate.posterior.evidence_error
    evidence_text = 'not available' if evidence is None else '{:.2f}'.format(evidence)
    if evidence_error is not None:
        evidence_text += ' ± {:.2f}'.format(evidence_error)
    analysis = candidate.metadata['analysis']
    rows = []
    for component in candidate.metadata['components']:
        rows.append(
            '<tr><td>{}</td><td>{}</td><td>{}</td></tr>'.format(
                escape(owner_labels.get(
                    component['owner_key'], component['owner_key']
                )),
                escape(component['role']), escape(component['model_name']),
            )
        )
    badge = ' <span class="default-badge">default</span>' if is_default else ''
    status = 'current input' if candidate.is_current(dataset) else 'input changed'
    return '''
<section class="candidate-head">
  <h2 style="margin:0 0 6px;font-size:18px">{}{}</h2>
  <p>Evidence {} · χ² {:.1f} · {} residuals · {}</p>
  <table class="component-table">
    <thead><tr><th>Owner</th><th>Role</th><th>Model</th></tr></thead>
    <tbody>{}</tbody>
  </table>
</section>
'''.format(
        escape(_candidate_label(candidate, owner_labels)), badge,
        evidence_text, analysis['chi_squared'], analysis['n_residuals'],
        status, ''.join(rows),
    )


def _owner_labels(candidate):
    targets = candidate.metadata.get('analysis', {}).get('targets', {})
    return {
        key: values.get('main_id') or values.get('sdbid') or key
        for key, values in targets.items()
    }


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
    return labels


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
