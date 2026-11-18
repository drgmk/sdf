"""Hierarchy-aware photometry exports and disk observability diagnostics."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from bokeh.embed import components
from bokeh.layouts import column
from bokeh.models import (
    ColumnDataSource, HoverTool, OpenURL, Span, TabPanel, Tabs, TapTool,
)
from bokeh.plotting import figure
from bokeh.resources import CDN
import jinja2
import numpy as np

from .candidate_library import CandidateLibrary
from .filter import Filter
from .spectrum import ObsSpectrum


@dataclass(frozen=True)
class FluxSizeBand:
    name: str
    beam_arcsec: float
    point_source_sensitivity_mjy: float


FLUX_SIZE_BANDS = (
    FluxSizeBand('MIRI.F1000W', 0.328, 0.48e-3),
    FluxSizeBand('MIRI.F1800W', 0.591, 0.23e-3),
    FluxSizeBand('MIRI.F2550W', 0.803, 12.2e-3),
    FluxSizeBand('PACS70', 6.0, 2000e-3),
    FluxSizeBand('WAV880', 1.0, 20e-3),
)


PHOTOMETRY_FIELDS = (
    'package', 'candidate_id', 'data_current', 'observation_id',
    'contributors', 'element_index', 'filter', 'wavelength_micron',
    'observed_jy', 'uncertainty_jy', 'upper_limit', 'ignored',
    'model_total_p16_jy', 'model_total_median_jy', 'model_total_p84_jy',
    'component_id', 'owner_key', 'role', 'model_name',
    'component_p16_jy', 'component_median_jy', 'component_p84_jy',
)


def write_photometry_csv(projection, path):
    """Write one long-form row per observed datum and model component."""

    path = Path(path)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=PHOTOMETRY_FIELDS)
        writer.writeheader()
        writer.writerows(photometry_rows(projection))
    return path


def photometry_rows(projection):
    rows = []
    for system in projection.systems:
        candidate = _load_candidate(projection, system)
        arrays = candidate.analysis_arrays
        analysis = candidate.metadata.get('analysis', {})
        fitted_components = analysis.get('components', ())
        for scope in analysis.get('data_scopes', ()):
            if scope.get('kind') != 'Photometry':
                continue
            contributor_keys = set(scope.get('contributors', ()))
            start = int(scope['element_start'])
            stop = int(scope['element_stop'])
            for element in range(start, stop):
                common = {
                    'package': system.package_path,
                    'candidate_id': system.candidate_id,
                    'data_current': system.data_current,
                    'observation_id': scope['observation_id'],
                    'contributors': ';'.join(scope.get('contributors', ())),
                    'element_index': element - start,
                    'filter': str(arrays['filter'][element]),
                    'wavelength_micron': arrays['wavelength_micron'][element],
                    'observed_jy': arrays['observed_flux'][element],
                    'uncertainty_jy': arrays['observed_uncertainty'][element],
                    'upper_limit': bool(arrays['observed_upper_limit'][element]),
                    'ignored': bool(arrays['observed_ignore'][element]),
                    'model_total_p16_jy': arrays['prediction_p16'][element],
                    'model_total_median_jy': arrays[
                        'prediction_median'
                    ][element],
                    'model_total_p84_jy': arrays['prediction_p84'][element],
                }
                for index, component in enumerate(fitted_components):
                    if component['owner_key'] not in contributor_keys:
                        continue
                    rows.append({
                        **common,
                        'component_id': component['component_id'],
                        'owner_key': component['owner_key'],
                        'role': component['role'],
                        'model_name': component['model_name'],
                        'component_p16_jy': arrays[
                            'component_prediction_p16'
                        ][index, element],
                        'component_median_jy': arrays[
                            'component_prediction_median'
                        ][index, element],
                        'component_p84_jy': arrays[
                            'component_prediction_p84'
                        ][index, element],
                    })
    return tuple(rows)


def flux_size_rows(projection, bands=FLUX_SIZE_BANDS):
    """Project retained median component spectra into planning bands."""

    rows = []
    filters = {band.name: Filter.get(band.name) for band in bands}
    for system in projection.systems:
        candidate = _load_candidate(projection, system)
        arrays = candidate.analysis_arrays
        wavelength = arrays.get('spectrum_wavelength_micron')
        spectra = arrays.get('component_spectrum_median')
        if wavelength is None or spectra is None:
            continue
        components = candidate.metadata.get('analysis', {}).get(
            'components', ()
        )
        component_index = {
            component['component_id']: index
            for index, component in enumerate(components)
        }
        for target in system.targets:
            parallax = _median(target.parallax_arcsec)
            if parallax is None or parallax <= 0:
                continue
            stellar_indices = [
                index for index, component in enumerate(components)
                if component['owner_key'] == target.target_key
                and component['role'] == 'stellar'
            ]
            if not stellar_indices:
                continue
            stellar_spectrum = np.sum(spectra[stellar_indices], axis=0)
            for disk in target.disks:
                radius_au = _median(disk.blackbody_radius_au)
                disk_index = component_index.get(disk.component_id)
                if radius_au is None or radius_au <= 0 or disk_index is None:
                    continue
                radius_arcsec = radius_au * parallax
                disk_spectrum = spectra[disk_index]
                for band in bands:
                    disk_flux_jy = _synthetic_flux(
                        filters[band.name], wavelength, disk_spectrum
                    )
                    stellar_flux_jy = _synthetic_flux(
                        filters[band.name], wavelength, stellar_spectrum
                    )
                    if not _positive(
                            radius_arcsec, disk_flux_jy, stellar_flux_jy):
                        continue
                    rows.append({
                        'band': band.name,
                        'beam_arcsec': band.beam_arcsec,
                        'sensitivity_mjy': (
                            band.point_source_sensitivity_mjy
                        ),
                        'system': system.label,
                        'target': target.label,
                        'component_id': disk.component_id,
                        'direct_member': target.direct_member,
                        'url': system.url,
                        'radius_arcsec': radius_arcsec,
                        'disk_flux_mjy': 1e3 * disk_flux_jy,
                        'stellar_flux_mjy': 1e3 * stellar_flux_jy,
                        'disk_contrast_per_beam': (
                            disk_flux_jy / stellar_flux_jy /
                            (radius_arcsec / band.beam_arcsec) ** 2
                        ),
                        'disk_surface_brightness_mjy_arcsec2': (
                            1e3 * disk_flux_jy /
                            (np.pi * radius_arcsec ** 2)
                        ),
                    })
    return tuple(rows)


def write_flux_size_plot(projection, output_directory):
    script, div = flux_size_plot_components(projection)
    env = jinja2.Environment(
        autoescape=True,
        loader=jinja2.PackageLoader('sdf', package_path='www/templates'),
    )
    document = env.get_template('joint_sample_flux_size.html').render(
        projection=projection,
        plot_script=script,
        plot_div=div,
        bokeh_js=CDN.render_js(),
        bokeh_css=CDN.render_css(),
        creation_time=_created(),
    )
    path = Path(output_directory) / 'fnuvsr.html'
    path.write_text(document, encoding='utf-8')
    return path


def flux_size_plot_components(projection):
    """Return embeddable Bokeh components for disk flux/size diagnostics."""

    rows = flux_size_rows(projection)
    panels = {'flux': [], 'contrast': [], 'surface_brightness': []}
    for band in FLUX_SIZE_BANDS:
        selected = [row for row in rows if row['band'] == band.name]
        panels['flux'].append(TabPanel(
            child=_flux_size_figure(
                selected, band, 'disk_flux_mjy', 'Disk flux / mJy',
                sensitivity=band.point_source_sensitivity_mjy,
            ),
            title=band.name,
        ))
        panels['contrast'].append(TabPanel(
            child=_flux_size_figure(
                selected, band, 'disk_contrast_per_beam',
                'Disk/star contrast per beam',
            ),
            title=band.name,
        ))
        panels['surface_brightness'].append(TabPanel(
            child=_flux_size_figure(
                selected, band, 'disk_surface_brightness_mjy_arcsec2',
                'Disk surface brightness / mJy arcsec⁻²',
                sensitivity=(
                    band.point_source_sensitivity_mjy /
                    (1.133 * band.beam_arcsec ** 2)
                ),
            ),
            title=band.name,
        ))
    layout = column(
        Tabs(tabs=panels['flux'], sizing_mode='stretch_width'),
        Tabs(tabs=panels['contrast'], sizing_mode='stretch_width'),
        Tabs(
            tabs=panels['surface_brightness'],
            sizing_mode='stretch_width',
        ),
        sizing_mode='stretch_width',
    )
    return components(layout)


def _flux_size_figure(rows, band, y_name, y_label, sensitivity=None):
    title = '{} — {} disk components'.format(band.name, len(rows))
    plot = figure(
        title=title,
        x_axis_label='Disk blackbody radius / arcsec',
        y_axis_label=y_label,
        x_axis_type='log', y_axis_type='log',
        height=570, sizing_mode='stretch_width',
        tools='pan,wheel_zoom,box_zoom,box_select,tap,reset,save',
        active_scroll='wheel_zoom',
    )
    plot.add_layout(Span(
        location=band.beam_arcsec, dimension='height',
        line_color='#c94b40', line_alpha=0.5, line_width=3,
    ))
    if sensitivity is not None and sensitivity > 0:
        plot.add_layout(Span(
            location=sensitivity, dimension='width',
            line_color='#438a53', line_alpha=0.5, line_width=3,
        ))
    if not rows:
        return plot
    data = {key: [row[key] for row in rows] for key in rows[0]}
    data['colour'] = [
        '#2f6fad' if row['direct_member'] else '#a7b0bb'
        for row in rows
    ]
    data['marker'] = [
        'circle' if row['direct_member'] else 'triangle' for row in rows
    ]
    source = ColumnDataSource(data)
    renderer = plot.scatter(
        'radius_arcsec', y_name, source=source,
        marker='marker', size=10, color='colour', alpha=0.8,
    )
    plot.add_tools(HoverTool(
        renderers=[renderer],
        tooltips=[
            ('target', '@target'), ('system', '@system'),
            ('component', '@component_id'),
        ],
    ))
    tap = plot.select(type=TapTool)
    tap.callback = OpenURL(url='@url')
    return plot


def _load_candidate(projection, system):
    directory = projection.root / Path(system.result_path).parent
    library = CandidateLibrary(projection.root / system.package_path)
    return library.load(directory)


def _synthetic_flux(filter_, wavelength, flux):
    spectrum = ObsSpectrum(wavelength=wavelength, fnujy=flux)
    value, _ = filter_.synthphot(spectrum)
    return float(value)


def _median(value):
    if not isinstance(value, dict) or value.get('median') is None:
        return None
    return float(value['median'])


def _positive(*values):
    return all(value is not None and np.isfinite(value) and value > 0
               for value in values)


def _created():
    return datetime.now(timezone.utc).strftime('%d/%m/%y %X UTC')
