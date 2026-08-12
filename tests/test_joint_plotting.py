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
