from .context import sdf
from .test_coordinator import FakeSampler, _fit_dataset
import sdf.coordinator
import sdf.fit_management
import sdf.model_selection


def _two_candidates(tmp_path):
    tmp_path.mkdir(parents=True, exist_ok=True)
    dataset = _fit_dataset(tmp_path)
    first = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('kurucz-0.0',)
    )
    second = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('phoenix_sol',)
    )
    coordinator = sdf.coordinator.FitCoordinator(
        tmp_path, sampler=FakeSampler()
    )
    first_result = coordinator.fit_candidate(dataset, first)
    # Reuse the simple test observation with a different real model family.
    second_result = coordinator.fit_candidate(dataset, second)
    return dataset, first_result, second_result


def test_available_fits_report_default_and_display_order(tmp_path):
    _, first, second = _two_candidates(tmp_path)
    manager = sdf.fit_management.FitManager(tmp_path)

    selection = manager.select(
        default=second.path.name,
        displayed=(first.path.name, second.path.name),
    )
    available = manager.available()

    assert selection.displayed == (second.candidate_id, first.candidate_id)
    by_id = {fit.candidate_id: fit for fit in available}
    assert by_id[second.candidate_id].is_default
    assert by_id[second.candidate_id].display_position == 0
    assert by_id[first.candidate_id].display_position == 1


def test_select_resolves_unambiguous_storage_prefix(tmp_path):
    _, first, second = _two_candidates(tmp_path)
    manager = sdf.fit_management.FitManager(tmp_path)
    prefix = second.path.name[:second.path.name.index('--') + 4]

    selection = manager.select(default=prefix, displayed=(first.path.name,))

    assert selection.default == second.candidate_id
    assert selection.displayed == (second.candidate_id, first.candidate_id)


def test_select_updates_existing_result_index(tmp_path):
    index = sdf.result_index.ResultIndex(tmp_path / 'sdf-results.sqlite')
    index.rebuild(tmp_path)
    dataset, first, second = _two_candidates(tmp_path / 'package')
    manager = sdf.fit_management.FitManager(tmp_path / 'package')

    manager.select(default=first.path.name, displayed=(second.path.name,))

    rows = index.candidates_for_target(dataset.targets[0].key)
    assert rows[0]['candidate_id'] == first.candidate_id
    assert rows[0]['is_default'] == 1
    assert rows[0]['display_position'] == 0
