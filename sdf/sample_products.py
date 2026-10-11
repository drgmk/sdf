"""Hierarchy-aware sample summaries built from the retained-fit index."""

from __future__ import annotations

import ast
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path

from astropy.table import Table
from bokeh.embed import components
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, HoverTool, OpenURL, TapTool
from bokeh.plotting import figure
from bokeh.resources import CDN
import jinja2
import numpy as np

from .sample_diagnostics import write_flux_size_plot, write_photometry_csv


@dataclass(frozen=True)
class SampleDisk:
    component_id: str
    model_name: str
    temperature_k: dict | None
    fractional_luminosity: dict | None
    blackbody_radius_au: dict | None


@dataclass(frozen=True)
class SampleStar:
    component_id: str
    model_name: str
    temperature_k: dict | None
    luminosity_lsun: dict | None
    radius_rsun: dict | None


@dataclass(frozen=True)
class SampleTarget:
    target_key: str
    main_id: str | None
    sdbid: str | None
    direct_member: bool
    parallax_arcsec: dict | None
    ra_deg: float | None
    dec_deg: float | None
    spectral_type: str | None
    stars: tuple[SampleStar, ...]
    disks: tuple[SampleDisk, ...]

    @property
    def label(self):
        return self.main_id or self.sdbid or self.target_key

    @property
    def distance_pc(self):
        parallax = _median(self.parallax_arcsec)
        return None if parallax is None or parallax <= 0 else 1.0 / parallax

    @property
    def ra_hours(self):
        return None if self.ra_deg is None else self.ra_deg / 15.0


@dataclass(frozen=True)
class SampleSystem:
    package_path: str
    candidate_id: str
    result_path: str
    display_name: str
    evidence: float | None
    data_current: bool
    matched_observations: tuple[str, ...]
    targets: tuple[SampleTarget, ...]
    url: str

    @property
    def label(self):
        if self.matched_observations:
            return ' / '.join(self.matched_observations)
        return ' + '.join(target.label for target in self.targets)


@dataclass(frozen=True)
class SampleProjection:
    name: str
    identifiers: tuple[str, ...] | None
    systems: tuple[SampleSystem, ...]
    unmatched: tuple[str, ...]
    root: Path

    @property
    def target_count(self):
        return sum(len(system.targets) for system in self.systems)


def project_sample(index, name, identifiers=None, *, output_directory=None,
                   sed_base_url=None):
    """Load hierarchical default results for a sample-sized identifier set."""

    identifiers = (
        None if identifiers is None else
        tuple(dict.fromkeys(str(value) for value in identifiers))
    )
    requested = set(identifiers or ())
    packages = index.default_packages(identifiers)
    observations = index.observations_for_packages(
        row['package_path'] for row in packages
    )
    observations_by_package = {}
    for row in observations:
        package = observations_by_package.setdefault(
            row['package_path'], {}
        )
        package.setdefault(row['observation_id'], []).append(row['target_key'])
    matched = set()
    systems = []
    root = index.root()
    for package in packages:
        payload = json.loads(
            (root / package['result_path']).read_text(encoding='utf-8')
        )
        target_metadata = _package_target_metadata(
            root / package['package_path']
        )
        analysis = payload.get('analysis', {})
        fitted = {
            component['component_id']: component
            for component in analysis.get('components', ())
        }
        targets = []
        for target_key, values in sorted(
                analysis.get('targets', {}).items()):
            aliases = {
                value for value in (
                    target_key, values.get('main_id'), values.get('sdbid')
                ) if value is not None
            }
            direct = (
                bool(aliases & requested) if identifiers is not None else True
            )
            matched.update(aliases & requested)
            stars = tuple(
                _star(component, fitted.get(component['component_id'], {}))
                for component in values.get('stellar_components', ())
            )
            disks = tuple(
                SampleDisk(
                    component_id=component['component_id'],
                    model_name=component.get('model_name', ''),
                    temperature_k=component.get('temperature_k'),
                    fractional_luminosity=component.get(
                        'fractional_luminosity'
                    ),
                    blackbody_radius_au=component.get('blackbody_radius_au'),
                )
                for component in values.get('disk_components', ())
            )
            targets.append(SampleTarget(
                target_key=target_key, main_id=values.get('main_id'),
                sdbid=values.get('sdbid'), direct_member=direct,
                parallax_arcsec=values.get('parallax_arcsec'),
                ra_deg=target_metadata.get(target_key, {}).get('ra_deg'),
                dec_deg=target_metadata.get(target_key, {}).get('dec_deg'),
                spectral_type=target_metadata.get(target_key, {}).get(
                    'spectral_type'
                ),
                stars=stars, disks=disks,
            ))
        package_observations = observations_by_package.get(
            package['package_path'], {}
        )
        matched_observations = tuple(sorted(
            set(package_observations) & requested
        ))
        matched.update(matched_observations)
        systems.append(SampleSystem(
            package_path=package['package_path'],
            candidate_id=package['default_candidate_id'],
            result_path=package['result_path'],
            display_name=package['display_name'],
            evidence=package['evidence'],
            data_current=bool(package['data_current']),
            matched_observations=matched_observations,
            targets=tuple(targets),
            url=_package_url(
                root, package['package_path'], output_directory, sed_base_url
            ),
        ))
    return SampleProjection(
        name=name, identifiers=identifiers, systems=tuple(systems),
        unmatched=tuple(sorted(requested - matched)), root=root,
    )


def write_sample_products(index, sample, identifiers, output_directory, *,
                          table=True, plots=True, sed_base_url=None):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    projection = project_sample(
        index, sample, identifiers, output_directory=output_directory,
        sed_base_url=sed_base_url,
    )
    outputs = {}
    if table:
        outputs['votable'] = write_sample_votable(
            projection, output_directory / '{}.xml'.format(sample)
        )
        outputs['photometry'] = write_photometry_csv(
            projection,
            output_directory / '{}_photometry.csv'.format(sample),
        )
    if plots:
        outputs['plots'] = write_sample_plots(projection, output_directory)
        outputs['flux_size'] = write_flux_size_plot(
            projection, output_directory
        )
    if table:
        outputs['table'] = write_sample_table(
            projection, output_directory,
            has_photometry=True, has_plots=plots, has_flux_size=plots,
        )
    return projection, outputs


def write_sample_table(projection, output_directory, *, has_photometry=False,
                       has_plots=False, has_flux_size=False):
    output_directory = Path(output_directory)
    table = _render_template(
        'joint_sample_table.html', projection=projection,
        rows=sample_table_rows(projection),
        has_photometry=has_photometry, has_plots=has_plots,
        has_flux_size=has_flux_size,
    )
    path = output_directory / 'index.html'
    path.write_text(table, encoding='utf-8')
    return path


def sample_table_rows(projection):
    """Return compact, sortable science rows for a sample projection.

    Disk values remain grouped in one target row, matching the established
    sample-table presentation.  Numeric sorting uses the first displayed
    (warmest) disk component when more than one is present.
    """

    rows = []
    for system in projection.systems:
        for target in system.targets:
            star = target.stars[0] if target.stars else None
            luminosity = _median(
                None if star is None else star.luminosity_lsun
            )
            disks = tuple(sorted(({
                'component': disk,
                'temperature_k': _median(disk.temperature_k),
                'log_fractional_luminosity': _log10(
                    _median(disk.fractional_luminosity)
                ),
                'blackbody_radius_au': _median(disk.blackbody_radius_au),
            } for disk in target.disks), key=_disk_sort_key))
            first_disk = disks[0] if disks else {}
            rows.append({
                'system': system,
                'target': target,
                'star': star,
                'disks': disks,
                'temperature_k': _median(
                    None if star is None else star.temperature_k
                ),
                'log_luminosity': _log10(luminosity),
                'distance_pc': target.distance_pc,
                'disk_temperature_sort': first_disk.get('temperature_k'),
                'log_fractional_luminosity_sort': first_disk.get(
                    'log_fractional_luminosity'
                ),
                'blackbody_radius_sort': first_disk.get(
                    'blackbody_radius_au'
                ),
            })
    return tuple(rows)


def write_sample_votable(projection, path):
    table = sample_votable(projection)
    table.write(path, format='votable', overwrite=True)
    return Path(path)


def sample_votable(projection):
    """Build the hierarchy-aware sample result table in memory."""

    rows = []
    for system in projection.systems:
        for target in system.targets:
            star = target.stars[0] if target.stars else None
            disks = target.disks or (None,)
            for disk in disks:
                rows.append((
                    system.package_path, system.candidate_id,
                    target.target_key, target.main_id or '', target.sdbid or '',
                    target.direct_member, '' if star is None else star.component_id,
                    '' if star is None else star.model_name,
                    _median(None if star is None else star.temperature_k),
                    _median(None if star is None else star.luminosity_lsun),
                    '' if disk is None else disk.component_id,
                    '' if disk is None else disk.model_name,
                    _median(None if disk is None else disk.temperature_k),
                    _median(None if disk is None else disk.fractional_luminosity),
                    _median(None if disk is None else disk.blackbody_radius_au),
                    system.url,
                ))
    names = (
        'package', 'candidate_id', 'target_key', 'main_id', 'sdbid',
        'direct_sample_member', 'stellar_component_id', 'stellar_model',
        'stellar_temperature_k', 'stellar_luminosity_lsun',
        'disk_component_id', 'disk_model', 'disk_temperature_k',
        'disk_fractional_luminosity', 'disk_blackbody_radius_au', 'url',
    )
    dtypes = (
        'U256', 'U256', 'U256', 'U256', 'U256', bool,
        'U256', 'U256', float, float, 'U256', 'U256', float, float,
        float, 'U1024',
    )
    table = Table(rows=rows, names=names, dtype=dtypes)
    table.meta['updated'] = _created()
    return table


def write_sample_plots(projection, output_directory):
    script, div = sample_plot_components(projection)
    document = _render_template(
        'joint_sample_plot.html', projection=projection,
        plot_script=script, plot_div=div,
        bokeh_js=CDN.render_js(), bokeh_css=CDN.render_css(),
    )
    path = Path(output_directory) / 'hr.html'
    path.write_text(document, encoding='utf-8')
    return path


def sample_plot_components(projection):
    """Return embeddable Bokeh components for the HR and disk plots."""

    hr_data = {
        'temperature': [], 'luminosity': [], 'target': [], 'system': [],
        'model': [], 'url': [], 'colour': [], 'marker': [],
    }
    disk_data = {
        'temperature': [], 'fraction': [], 'target': [], 'system': [],
        'model': [], 'url': [], 'colour': [], 'marker': [],
    }
    for system in projection.systems:
        for target in system.targets:
            colour = '#2f6fad' if target.direct_member else '#a7b0bb'
            marker = 'circle' if target.direct_member else 'triangle'
            for star in target.stars:
                temperature = _median(star.temperature_k)
                luminosity = _median(star.luminosity_lsun)
                if _positive(temperature, luminosity):
                    _append_plot_row(
                        hr_data, temperature, luminosity, target, system,
                        star.model_name, colour, marker,
                    )
            for disk in target.disks:
                temperature = _median(disk.temperature_k)
                fraction = _median(disk.fractional_luminosity)
                if _positive(temperature, fraction):
                    _append_plot_row(
                        disk_data, temperature, fraction, target, system,
                        disk.model_name, colour, marker,
                    )
    hr = _sample_figure(
        hr_data, 'HR diagram', 'Effective temperature / K',
        'Stellar luminosity / Solar', reverse_x=True,
    )
    disk = _sample_figure(
        disk_data, 'Owned disk components', 'Disk temperature / K',
        'Disk fractional luminosity', reverse_x=False,
    )
    layout = gridplot([[hr, disk]], sizing_mode='stretch_width')
    return components(layout)


def _sample_figure(data, title, x_label, y_label, reverse_x):
    if not data['temperature']:
        return figure(title='{} — no fitted values'.format(title), height=500)
    x = np.asarray(data['temperature'], dtype=float)
    y = np.asarray(data['luminosity' if 'luminosity' in data else 'fraction'])
    x_range = (float(np.max(x) * 1.15), float(np.min(x) / 1.15)) \
        if reverse_x else (float(np.min(x) / 1.5), float(np.max(x) * 1.5))
    plot = figure(
        title='{} ({} components)'.format(title, len(x)),
        x_axis_label=x_label, y_axis_label=y_label,
        x_axis_type='log', y_axis_type='log', x_range=x_range,
        height=600, sizing_mode='stretch_width',
        tools='pan,wheel_zoom,box_zoom,box_select,tap,reset,save',
        active_scroll='wheel_zoom',
    )
    source = ColumnDataSource(data)
    renderer = plot.scatter(
        'temperature', 'luminosity' if 'luminosity' in data else 'fraction',
        source=source, marker='marker', size=10, color='colour', alpha=0.8,
    )
    plot.add_tools(HoverTool(
        renderers=[renderer],
        tooltips=[
            ('target', '@target'), ('system', '@system'), ('model', '@model'),
        ],
    ))
    tap = plot.select(type=TapTool)
    tap.callback = OpenURL(url='@url')
    return plot


def _append_plot_row(data, x, y, target, system, model, colour, marker):
    data['temperature'].append(x)
    y_name = 'luminosity' if 'luminosity' in data else 'fraction'
    data[y_name].append(y)
    data['target'].append(target.label)
    data['system'].append(system.label)
    data['model'].append(model)
    data['url'].append(system.url)
    data['colour'].append(colour)
    data['marker'].append(marker)


def _star(derived, fitted):
    return SampleStar(
        component_id=derived['component_id'],
        model_name=derived.get('model_name') or fitted.get('model_name', ''),
        temperature_k=fitted.get('parameters', {}).get('Teff'),
        luminosity_lsun=derived.get('luminosity_lsun'),
        radius_rsun=derived.get('radius_rsun'),
    )


def _package_target_metadata(package):
    """Read lightweight catalogue metadata from raw-photometry headers."""

    metadata = {}
    for path in Path(package).glob('*-rawphot.txt'):
        keywords = _rawphot_header(path)
        key = _target_key(keywords)
        if key is None:
            continue
        metadata[key] = {
            'ra_deg': _optional_float(keywords.get('raj2000')),
            'dec_deg': _optional_float(keywords.get('dej2000')),
            'spectral_type': _optional_text(keywords.get('sp_type')),
        }
    return metadata


def _rawphot_header(path):
    """Parse only the leading IPAC keywords, avoiding the photometry body."""

    values = {}
    with Path(path).open(encoding='utf-8') as handle:
        for line in handle:
            if not line.startswith('\\'):
                break
            key, separator, value = line[1:].partition('=')
            if not separator:
                continue
            value = value.strip()
            try:
                values[key.strip()] = ast.literal_eval(value)
            except (SyntaxError, ValueError):
                values[key.strip()] = value.strip("'\"")
    return values


def _target_key(keywords):
    main_id = _optional_text(keywords.get('main_id'))
    if main_id is not None:
        return 'simbad:' + ' '.join(main_id.split())
    sdbid = _optional_text(keywords.get('id'))
    return None if sdbid is None else 'sdb:' + sdbid


def _optional_text(value):
    if value is None:
        return None
    value = str(value).strip()
    return None if value in {'', 'None'} else value


def _optional_float(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _log10(value):
    return None if value is None or value <= 0 else math.log10(value)


def _disk_sort_key(disk):
    temperature = disk['temperature_k']
    return (temperature is None, 0.0 if temperature is None else -temperature)


def _package_url(root, package_path, output_directory, sed_base_url):
    if sed_base_url:
        return '{}/{}/'.format(
            sed_base_url.rstrip('/'), package_path.strip('/')
        )
    if output_directory is None:
        return '{}/index.html'.format(package_path.strip('/'))
    return os.path.relpath(
        root / package_path / 'index.html', Path(output_directory)
    )


def _render_template(name, **values):
    env = jinja2.Environment(
        autoescape=True,
        loader=jinja2.PackageLoader('sdf', package_path='www/templates'),
    )
    env.filters['summary'] = _summary_text
    return env.get_template(name).render(
        creation_time=_created(), **values,
    )


def _summary_text(value, precision=2):
    if value is None:
        return '—'
    median = _median(value)
    if median is None:
        return '—'
    minus = value.get('minus')
    plus = value.get('plus')
    if minus is None or plus is None:
        return ('{:.' + str(precision) + 'g}').format(median)
    template = '{:.' + str(precision) + 'g} −{:.' + str(precision) \
        + 'g}/+{:.' + str(precision) + 'g}'
    return template.format(median, minus, plus)


def summary_text(value, precision=2):
    """Format a retained median and asymmetric uncertainty for HTML tables."""

    return _summary_text(value, precision=precision)


def _median(value):
    if not isinstance(value, dict):
        return None
    median = value.get('median')
    return None if median is None else float(median)


def _positive(*values):
    return all(value is not None and np.isfinite(value) and value > 0
               for value in values)


def _created():
    return datetime.now(timezone.utc).strftime('%d/%m/%y %X UTC')
