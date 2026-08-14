from dataclasses import replace

import numpy as np

from .context import sdf
from .test_coordinator import FakeSampler, _fit_dataset
from .test_result_index import _retained_package
import sdf.sample_products


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
    _, dataset, _ = _retained_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'sdf-results.sqlite')
    index.rebuild(tmp_path)

    projection = sdf.sample_products.project_sample(
        index, 'test-sample', (dataset.targets[0].sdbid,),
    )

    target = projection.systems[0].targets[0]
    assert target.direct_member
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
    assert set(outputs) == {'table', 'votable', 'plots'}
    html = outputs['table'].read_text(encoding='utf-8')
    assert 'Physical target' in html
    assert 'sample member' in html
    assert outputs['votable'].is_file()
    assert 'Bokeh' in outputs['plots'].read_text(encoding='utf-8')


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
