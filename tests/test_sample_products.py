import csv
from dataclasses import replace

import numpy as np

from .context import sdf
from .test_coordinator import FakeSampler, _fit_dataset
from .test_result_index import _retained_package
import sdf.sample_products
import sdf.sample_diagnostics


def _retained_binary_package(tmp_path):
    package = tmp_path / 'binary-package'
    package.mkdir()
    source = package / 'AB-rawphot.txt'
    source.write_text('test input\n', encoding='utf-8')
    cached_model = sdf.model.PhotModel.read_model('kurucz-0.0')
    filter_name = next(
        name for name, bases in zip(
            cached_model.filters, cached_model.colour_bases
        ) if len(bases) == 0
    )
    targets = tuple(
        sdf.dataset.TargetIdentity(
            key='simbad:' + name, main_id=name, sdbid='sdb-' + name,
            keywords={'plx_value': 10.0, 'plx_err': 0.0},
        )
        for name in ('A', 'B')
    )
    photometry = sdf.photometry.Photometry(
        filters=np.array([filter_name]), fnujy=np.array([1.0]),
        e_fnujy=np.array([0.1]), upperlim=np.array([False]),
        ignore=np.array([False]), bibcode=np.array(['test']),
    )
    dataset = sdf.dataset.FitDataset(
        targets=targets,
        observations=(sdf.dataset.Observation(
            observation_id='sdb-AB', source=str(source.resolve()),
            target_key='simbad:AB',
            contributor_keys=('simbad:A', 'simbad:B'),
            photometry=photometry,
        ),),
    )
    candidate = sdf.model_selection.ModelCandidate.from_owned_models((
        ('simbad:A', 'kurucz-0.0'),
        ('simbad:A', 'bb_disk_r'),
        ('simbad:B', 'kurucz-0.0'),
    ))
    coordinator = sdf.coordinator.FitCoordinator(
        package, sampler=FakeSampler()
    )
    result = coordinator.fit_candidate(dataset, candidate)
    coordinator.library.write_selection(sdf.candidate_library.FitSelection(
        default=result.candidate_id, displayed=(result.candidate_id,)
    ))
    return package, dataset


def test_projection_finds_system_by_observation_and_marks_context(tmp_path):
    package = tmp_path / 'package'
    package.mkdir()
    dataset = _fit_dataset(package)
    observation = replace(
        dataset.observations[0], observation_id='sdb-composite-ab'
    )
    dataset = replace(dataset, observations=(observation,))
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('kurucz-0.0',)
    )
    coordinator = sdf.coordinator.FitCoordinator(
        package, sampler=FakeSampler()
    )
    result = coordinator.fit_candidate(dataset, candidate)
    coordinator.library.write_selection(sdf.candidate_library.FitSelection(
        default=result.candidate_id, displayed=(result.candidate_id,)
    ))
    index = sdf.result_index.ResultIndex(tmp_path / 'sdf-results.sqlite')
    index.rebuild(tmp_path)
    observation_id = dataset.observations[0].observation_id

    projection = sdf.sample_products.project_sample(
        index, 'test-sample', (observation_id, 'missing'),
        output_directory=tmp_path / 'samples' / 'test-sample',
    )

    assert len(projection.systems) == 1
    assert projection.systems[0].matched_observations == (observation_id,)
    assert projection.systems[0].targets[0].direct_member is False
    assert projection.unmatched == ('missing',)
    assert projection.systems[0].url.endswith('package/index.html')


def test_projection_marks_direct_physical_target_member(tmp_path):
    package, dataset, _ = _retained_package(tmp_path)
    (package / 'identities.yml').write_text(
        'version: 1\ntargets:\n  {}:\n'
        '    aliases:\n      - HD 123\n      - HIP 107000\n'.format(
            dataset.targets[0].sdbid
        ),
        encoding='utf-8',
    )
    index = sdf.result_index.ResultIndex(tmp_path / 'sdf-results.sqlite')
    index.rebuild(tmp_path)

    projection = sdf.sample_products.project_sample(
        index, 'test-sample', (dataset.targets[0].sdbid,),
    )

    target = projection.systems[0].targets[0]
    assert target.direct_member
    assert target.aliases == ('HD 123', 'HIP 107000')
    assert 'HIP 107000' in sdf.sample_products.sample_table_rows(
        projection
    )[0]['filter_text']
    assert target.stars[0].temperature_k['median'] > 0


def test_writers_create_hierarchy_table_votable_and_plots(tmp_path):
    _, dataset, _ = _retained_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'sdf-results.sqlite')
    index.rebuild(tmp_path)
    output = tmp_path / 'sample-output'

    projection, outputs = sdf.sample_products.write_sample_products(
        index, 'test-sample', (dataset.targets[0].sdbid,), output,
    )

    assert projection.target_count == 1
    assert set(outputs) == {
        'table', 'votable', 'photometry', 'plots', 'flux_size',
    }
    html = outputs['table'].read_text(encoding='utf-8')
    assert 'Filter targets' in html
    assert 'T<sub>eff</sub>' in html
    assert 'sample member' not in html
    assert '>sdb-test<' not in html
    assert outputs['votable'].is_file()
    assert outputs['photometry'].is_file()
    assert 'Bokeh' in outputs['plots'].read_text(encoding='utf-8')
    assert 'Bokeh' in outputs['flux_size'].read_text(encoding='utf-8')


def test_sample_table_rows_are_science_columns_and_group_disks(tmp_path):
    _retained_binary_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'sdf-results.sqlite')
    index.rebuild(tmp_path)
    projection = sdf.sample_products.project_sample(
        index, 'binary-sample', ('sdb-AB',)
    )

    rows = sdf.sample_products.sample_table_rows(projection)

    assert len(rows) == 2
    assert rows[0]['target'].main_id == 'A'
    assert rows[0]['temperature_k'] > 0
    assert rows[0]['log_luminosity'] is not None
    assert rows[0]['disks'][0]['temperature_k'] > 0
    assert rows[1]['target'].main_id == 'B'
    assert rows[1]['disks'] == ()
    assert rows[1]['disk_temperature_sort'] is None


def test_system_plot_lines_only_connect_components_from_the_same_system():
    data = {
        'temperature': [6000.0, 4500.0, 3000.0],
        'luminosity': [2.0, 0.4, 0.1],
        'system_key': ['one', 'one', 'two'],
    }

    lines = sdf.sample_products._system_line_data(data)

    assert lines == {
        'xs': [[6000.0, 4500.0]],
        'ys': [[2.0, 0.4]],
    }


def test_system_lines_do_not_capture_point_navigation():
    data = {
        'temperature': [6000.0, 4500.0],
        'luminosity': [2.0, 0.4],
        'target': ['A', 'B'],
        'system': ['AB', 'AB'],
        'system_key': ['one', 'one'],
        'model': ['star', 'star'],
        'url': ['/system/one/', '/system/one/'],
        'colour': ['blue', 'blue'],
        'marker': ['circle', 'circle'],
    }

    plot = sdf.sample_products._sample_figure(
        data, 'HR diagram', 'temperature', 'luminosity', True,
    )
    tap = plot.select(type=sdf.sample_products.TapTool)

    assert len(plot.renderers) == 2
    assert tap.renderers == [plot.renderers[1]]
    assert tap.renderers[0].data_source.data['url'] == [
        '/system/one/', '/system/one/',
    ]


def test_composite_member_includes_both_contributors_and_owned_disk(tmp_path):
    _retained_binary_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'sdf-results.sqlite')
    index.rebuild(tmp_path)

    projection = sdf.sample_products.project_sample(
        index, 'binary-sample', ('sdb-AB',)
    )

    system = projection.systems[0]
    assert system.matched_observations == ('sdb-AB',)
    assert tuple(target.main_id for target in system.targets) == ('A', 'B')
    assert not any(target.direct_member for target in system.targets)
    assert len(system.targets[0].disks) == 1
    assert system.targets[1].disks == ()


def test_photometry_csv_is_long_form_and_keeps_component_ownership(tmp_path):
    _retained_binary_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'sdf-results.sqlite')
    index.rebuild(tmp_path)
    projection = sdf.sample_products.project_sample(
        index, 'binary-sample', ('sdb-AB',)
    )
    path = tmp_path / 'photometry.csv'

    sdf.sample_diagnostics.write_photometry_csv(projection, path)

    with path.open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 3
    assert {row['role'] for row in rows} == {'stellar', 'disk'}
    assert {row['owner_key'] for row in rows} == {'simbad:A', 'simbad:B'}
    assert {row['contributors'] for row in rows} == {'simbad:A;simbad:B'}
    assert all(float(row['model_total_median_jy']) > 0 for row in rows)


def test_flux_size_projection_uses_owned_stellar_and_disk_spectra(tmp_path):
    _retained_binary_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'sdf-results.sqlite')
    index.rebuild(tmp_path)
    projection = sdf.sample_products.project_sample(
        index, 'binary-sample', ('sdb-AB',)
    )

    rows = sdf.sample_diagnostics.flux_size_rows(projection)

    assert {row['band'] for row in rows} == {
        band.name for band in sdf.sample_diagnostics.FLUX_SIZE_BANDS
    }
    assert {row['target'] for row in rows} == {'A'}
    assert all(row['disk_flux_mjy'] > 0 for row in rows)
    assert all(row['stellar_flux_mjy'] > 0 for row in rows)
    assert all(row['radius_arcsec'] > 0 for row in rows)
