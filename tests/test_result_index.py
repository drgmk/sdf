import sqlite3

import pytest

from .context import sdf
from .test_coordinator import FakeSampler, _fit_dataset
import sdf.candidate_library
import sdf.coordinator
import sdf.model_selection
import sdf.result_index


def _retained_package(tmp_path):
    package = tmp_path / 'package'
    package.mkdir()
    dataset = _fit_dataset(package)
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
    return package, dataset, result


def test_rebuild_indexes_candidate_targets_components_and_selection(tmp_path):
    package, dataset, result = _retained_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')

    summary = index.rebuild(tmp_path)

    assert summary == sdf.result_index.IndexSummary(1, 1, 1, 1, 1, 5)
    rows = index.candidates_for_target(dataset.targets[0].key)
    assert len(rows) == 1
    assert rows[0]['candidate_id'] == result.candidate_id
    assert rows[0]['package_path'] == 'package'
    assert rows[0]['data_current'] == 1
    assert rows[0]['is_default'] == 1
    assert rows[0]['display_position'] == 0
    default = index.default_targets(('sdb-test',))
    assert len(default) == 1
    assert default[0]['target_key'] == dataset.targets[0].key
    assert default[0]['candidate_id'] == result.candidate_id
    parameters = index.default_parameters(('Test',))
    assert {
        (row['quantity_scope'], row['name']) for row in parameters
    } == {
        ('fitted', 'Teff'), ('fitted', 'logg'), ('fitted', 'norm'),
        ('derived', 'luminosity_at_1pc_lsun'),
        ('derived', 'radius_at_1pc_m'),
    }
    assert index.verify(tmp_path) == ()


def test_rebuild_indexes_exported_package_before_it_is_fit(tmp_path):
    from .test_dataset import _write_rawphot

    package = tmp_path / 'unfit'
    package.mkdir()
    _write_rawphot(
        package / 'target-rawphot.txt',
        sdbid='sdb-unfit', main_id='HD 456',
    )
    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')

    summary = index.rebuild(tmp_path)
    rows = index.catalog_packages()

    assert summary == sdf.result_index.IndexSummary(1, 1, 1, 0, 0, 0)
    assert rows == ({
        'package_path': 'unfit',
        'input_fingerprint': rows[0]['input_fingerprint'],
        'default_candidate_id': None,
        'updated_at': rows[0]['updated_at'],
        'candidate_count': 0,
        'has_current_candidate': 0,
        'default_current': 0,
        'target_labels': ('HD 456',),
        'sdbids': ('sdb-unfit',),
        'status': 'unfit',
    },)
    assert index.catalog_packages(query='456') == rows
    assert index.catalog_packages(identifiers=('sdb-unfit',)) == rows


def test_rebuild_skips_export_package_without_usable_observations(tmp_path):
    from .test_dataset import _write_empty_rawphot

    package = tmp_path / 'incomplete'
    package.mkdir()
    _write_empty_rawphot(
        package / 'target-rawphot.txt',
        sdbid='sdb-incomplete', main_id='HD 15003',
    )
    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')

    summary = index.rebuild(tmp_path)

    assert summary == sdf.result_index.IndexSummary(0, 0, 0, 0, 0, 0)
    assert index.catalog_packages() == ()
    assert index.verify(tmp_path) == ()


def test_catalog_search_ignores_identifier_spacing_and_ranks_exact_first(
        tmp_path):
    from .test_dataset import _write_rawphot

    exact = tmp_path / 'exact'
    prefix = tmp_path / 'prefix'
    exact.mkdir()
    prefix.mkdir()
    _write_rawphot(
        exact / 'target-rawphot.txt',
        sdbid='sdb-exact', main_id='HD     38',
    )
    _write_rawphot(
        prefix / 'target-rawphot.txt',
        sdbid='sdb-prefix', main_id='HD 38138',
    )
    (exact / 'identities.yml').write_text(
        'version: 1\n'
        'targets:\n'
        '  sdb-exact:\n'
        '    aliases:\n'
        '      - HIP 473\n',
        encoding='utf-8',
    )
    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')
    index.rebuild(tmp_path)

    assert [row['package_path'] for row in
            index.catalog_packages(query='hd38')] == ['exact', 'prefix']
    assert [row['package_path'] for row in
            index.catalog_packages(query='HIP473')] == ['exact']


def test_identity_alias_change_reindexes_without_staling_fit(tmp_path):
    package, _, _ = _retained_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')
    index.rebuild(tmp_path)
    assert index.catalog_packages()[0]['status'] == 'current'

    (package / 'identities.yml').write_text(
        'version: 1\n'
        'targets:\n'
        '  sdb-test:\n'
        '    aliases:\n'
        '      - HIP 473\n',
        encoding='utf-8',
    )
    index.update_package(package)

    rows = index.catalog_packages(query='hip473')
    assert len(rows) == 1
    assert rows[0]['status'] == 'current'


def test_catalog_status_becomes_stale_after_input_change(tmp_path):
    package, dataset, _ = _retained_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')
    index.rebuild(tmp_path)
    assert index.catalog_packages()[0]['status'] == 'current'

    with open(dataset.observations[0].source, 'a', encoding='utf-8') as handle:
        handle.write('changed input\n')
    index.update_package(package)

    assert index.catalog_packages()[0]['status'] == 'stale'


def test_default_package_can_be_found_by_observation_scope(tmp_path):
    package, dataset, result = _retained_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')
    index.rebuild(tmp_path)
    observation_id = dataset.observations[0].observation_id

    packages = index.default_packages((observation_id,))
    observations = index.observations_for_packages(('package',))

    assert packages[0]['package_path'] == 'package'
    assert packages[0]['default_candidate_id'] == result.candidate_id
    assert observations == ({
        'package_path': 'package',
        'observation_id': observation_id,
        'target_key': dataset.targets[0].key,
        'position': 0,
    },)


def test_rebuild_marks_candidate_data_stale_after_input_change(tmp_path):
    package, dataset, _ = _retained_package(tmp_path)
    with open(dataset.observations[0].source, 'a', encoding='utf-8') as handle:
        handle.write('changed input\n')
    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')

    index.rebuild(tmp_path)

    rows = index.candidates_for_target(dataset.targets[0].key)
    assert rows[0]['data_current'] == 0


def test_default_queries_accept_sample_sized_identifier_lists(tmp_path):
    _, dataset, _ = _retained_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')
    index.rebuild(tmp_path)
    identifiers = tuple('missing-{}'.format(i) for i in range(1200)) + (
        'sdb-test',
    )

    assert index.default_targets(identifiers)[0]['target_key'] == \
        dataset.targets[0].key
    assert index.default_parameters(identifiers)


def test_verify_reports_artifact_removed_after_rebuild(tmp_path):
    _, _, result = _retained_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')
    index.rebuild(tmp_path)
    (result.path / 'result.json').unlink()

    problems = index.verify(tmp_path)

    assert problems == (
        'missing candidate artifact: {}'.format(
            (result.path / 'result.json').relative_to(tmp_path).as_posix()
        ),
    )


def test_rebuilt_index_passes_integrity_check(tmp_path):
    _retained_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')
    index.rebuild(tmp_path)

    with sqlite3.connect(index.path) as connection:
        assert connection.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'


def test_old_index_schema_requests_rebuild(tmp_path):
    _retained_package(tmp_path)
    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')
    index.rebuild(tmp_path)
    with sqlite3.connect(index.path) as connection:
        connection.execute(
            "UPDATE index_metadata SET value = '1' "
            "WHERE key = 'schema_version'"
        )

    with pytest.raises(ValueError, match='rebuild the index'):
        index.default_packages()
    assert any('rebuild the index' in problem for problem in index.verify())


def test_verify_reports_missing_index(tmp_path):
    index = sdf.result_index.ResultIndex(tmp_path / 'missing.sqlite')

    assert index.verify() == (
        'result index does not exist: {}'.format(index.path),
    )


def test_verify_reports_new_exported_package_after_rebuild(tmp_path):
    from .test_dataset import _write_rawphot

    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')
    index.rebuild(tmp_path)
    package = tmp_path / 'later'
    package.mkdir()
    _write_rawphot(
        package / 'later-rawphot.txt', sdbid='later', main_id='HD 789'
    )

    assert index.verify(tmp_path) == (
        'exported package absent from index: later',
    )


def test_coordinator_updates_nearest_existing_default_index(tmp_path):
    index = sdf.result_index.ResultIndex(tmp_path / 'sdf-results.sqlite')
    assert index.rebuild(tmp_path).candidates == 0
    package = tmp_path / 'nested' / 'package'
    package.mkdir(parents=True)
    dataset = _fit_dataset(package)
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('kurucz-0.0',)
    )

    coordinator = sdf.coordinator.FitCoordinator(
        package, sampler=FakeSampler()
    )
    result = coordinator.fit_candidate(dataset, candidate)

    assert coordinator.result_index.path == index.path
    assert index.summary().candidates == 1
    assert index.candidates_for_target(dataset.targets[0].key)[0][
        'candidate_id'
    ] == result.candidate_id


def test_no_result_index_is_created_implicitly(tmp_path):
    package, _, _ = _retained_package(tmp_path)

    assert sdf.result_index.ResultIndex.find_for_package(package) is None
    assert not (tmp_path / 'sdf-results.sqlite').exists()
