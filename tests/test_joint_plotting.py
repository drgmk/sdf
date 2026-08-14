from .context import sdf
from .test_coordinator import FakeSampler, _fit_dataset
import sdf.coordinator
import sdf.joint_plotting
import sdf.model_selection


def test_joint_results_page_contains_candidate_and_scoped_plot(tmp_path):
    dataset = _fit_dataset(tmp_path)
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('kurucz-0.0',)
    )
    coordinator = sdf.coordinator.FitCoordinator(
        tmp_path, sampler=FakeSampler()
    )
    coordinator.fit_candidate(dataset, candidate)

    output = sdf.joint_plotting.write_joint_results_page(
        tmp_path, dataset=dataset
    )
    html = output.read_text(encoding='utf-8')

    assert output == tmp_path / 'index.html'
    assert 'Bokeh' in html
    assert 'Test: kurucz-0.0' in html
    assert 'Evidence 0.00' in html
    assert 'Residual' in html
    assert 'current input' in html
    assert 'Physical contributors' in html
    assert 'Observed scopes' in html
    assert 'result.json' in html
    assert 'posterior-' in html
    assert 'fit-controls' not in html


def test_operator_results_page_extends_public_page_with_fit_drawer(tmp_path):
    dataset = _fit_dataset(tmp_path)
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('kurucz-0.0',)
    )
    sdf.coordinator.FitCoordinator(
        tmp_path, sampler=FakeSampler()
    ).fit_candidate(dataset, candidate)

    output = sdf.joint_plotting.write_joint_results_page(
        tmp_path, dataset=dataset, operator=True,
    )
    html = output.read_text(encoding='utf-8')

    assert output == tmp_path / 'operator.html'
    assert 'class="sed-page operator-mode"' in html
    assert 'id="fit-controls"' in html
    assert 'id="open-fit-controls"' in html
    assert 'data-target-key="simbad:Test"' in html
    assert 'data-fit-reference=' in html


def test_scope_colours_reuse_owner_colour_only_for_resolved_observations():
    scopes = (
        {'observation_id': 'A', 'contributors': ('a',)},
        {'observation_id': 'AB', 'contributors': ('a', 'b')},
        {'observation_id': 'B', 'contributors': ('b',)},
    )

    colours = sdf.joint_plotting._scope_colours(
        scopes, {'a': 'blue', 'b': 'orange'},
    )

    assert colours['A'] == 'blue'
    assert colours['B'] == 'orange'
    assert colours['AB'] not in {'blue', 'orange'}
