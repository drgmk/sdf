from .context import sdf
from .test_coordinator import FakeSampler, _fit_dataset
import pytest
import sdf.coordinator
import sdf.fit_management
import sdf.model_selection
import sdf.result_index


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


def test_exact_candidate_resolves_target_identity_and_disk_ownership(tmp_path):
    dataset = _fit_dataset(tmp_path)

    candidate, policies = sdf.fit_management.exact_candidate(
        dataset,
        stellar=(('Test', 'kurucz-0.0'),),
        disks=(
            ('sdb-test', 'modbb_disk_r'),
            ('simbad:Test', 'modbb_disk_r'),
        ),
    )

    assert candidate.model_names == (
        'kurucz-0.0', 'modbb_disk_r', 'modbb_disk_r'
    )
    assert {component.owner_key for component in candidate.components} == {
        'simbad:Test'
    }
    assert policies[0].search_mode == 'exact'
    assert policies[0].max_disks == 2
    assert policies[0].origin == 'explicit'


def test_exact_candidate_requires_one_stellar_model_per_target(tmp_path):
    dataset = _fit_dataset(tmp_path)

    with pytest.raises(ValueError, match='stellar model missing for: Test'):
        sdf.fit_management.exact_candidate(dataset, stellar=())


def test_refit_forces_sampler_for_current_candidate(tmp_path):
    dataset, first, _ = _two_candidates(tmp_path)
    sampler = FakeSampler()
    manager = sdf.fit_management.FitManager(tmp_path)

    result = manager.refit(
        first.path.name, sampler=sampler, dataset=dataset
    )

    assert result.candidate_id == first.candidate_id
    assert len(sampler.calls) == 1


def test_fit_exact_retains_explicit_policy(tmp_path):
    dataset = _fit_dataset(tmp_path)
    sampler = FakeSampler()
    manager = sdf.fit_management.FitManager(tmp_path)

    result = manager.fit_exact(
        stellar=(('Test', 'kurucz-0.0'),),
        sampler=sampler,
        dataset=dataset,
    )

    assert len(sampler.calls) == 1
    assert result.metadata['policy']['targets'][0]['search_mode'] == 'exact'
    assert result.metadata['policy']['targets'][0]['origin'] == 'explicit'


def test_preview_does_not_change_or_expand_public_selection(tmp_path):
    dataset, first, second = _two_candidates(tmp_path)
    manager = sdf.fit_management.FitManager(tmp_path)
    manager.select(default=first.path.name, displayed=(first.path.name,))

    public = manager.render(tmp_path / 'public.html', dataset=dataset)
    preview = manager.render(
        tmp_path / 'preview.html', fits=(second.path.name,), dataset=dataset
    )

    public_html = public.read_text(encoding='utf-8')
    preview_html = preview.read_text(encoding='utf-8')
    assert first.metadata['display_name'] in public_html
    assert second.metadata['display_name'] not in public_html
    assert second.metadata['display_name'] in preview_html
    assert manager.library.read_selection().displayed == (first.candidate_id,)


def test_fit_manager_can_render_operator_extension(tmp_path):
    dataset, first, _ = _two_candidates(tmp_path)
    manager = sdf.fit_management.FitManager(tmp_path)
    manager.select(default=first.path.name, displayed=(first.path.name,))

    output = manager.render(
        tmp_path / 'operator.html', dataset=dataset, operator=True,
    )

    assert 'id="fit-controls"' in output.read_text(encoding='utf-8')
