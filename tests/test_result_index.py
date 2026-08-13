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

    assert summary == sdf.result_index.IndexSummary(1, 1, 1, 1)
    rows = index.candidates_for_target(dataset.targets[0].key)
    assert len(rows) == 1
    assert rows[0]['candidate_id'] == result.candidate_id
    assert rows[0]['package_path'] == 'package'
    assert rows[0]['data_current'] == 1
    assert rows[0]['is_default'] == 1
    assert rows[0]['display_position'] == 0
    assert index.verify(tmp_path) == ()


def test_rebuild_marks_candidate_data_stale_after_input_change(tmp_path):
    package, dataset, _ = _retained_package(tmp_path)
    with open(dataset.observations[0].source, 'a', encoding='utf-8') as handle:
        handle.write('changed input\n')
    index = sdf.result_index.ResultIndex(tmp_path / 'results.sqlite')

    index.rebuild(tmp_path)

    rows = index.candidates_for_target(dataset.targets[0].key)
    assert rows[0]['data_current'] == 0


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
