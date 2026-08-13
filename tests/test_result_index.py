import sqlite3

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

    assert summary == sdf.result_index.IndexSummary(1, 1, 1, 1, 5)
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


def test_verify_reports_missing_index(tmp_path):
    index = sdf.result_index.ResultIndex(tmp_path / 'missing.sqlite')

    assert index.verify() == (
        'result index does not exist: {}'.format(index.path),
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
