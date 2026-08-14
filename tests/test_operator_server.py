from .context import sdf
from .test_fit_management import _two_candidates
import sdf.operator_server
import sdf.operator_state


def _service(tmp_path):
    dataset, first, second = _two_candidates(tmp_path / 'package')
    manager = sdf.fit_management.FitManager(tmp_path / 'package')
    manager.select(default=first.path.name, displayed=(first.path.name,))
    state = sdf.operator_state.OperatorState(tmp_path / 'operator.sqlite')
    service = sdf.operator_server.OperatorService(
        tmp_path / 'package', state, actor='Grant', token='test-token',
        dataset=dataset,
    )
    return service, dataset, first, second


def test_operator_page_contains_live_policy_and_all_retained_fits(tmp_path):
    service, _, first, second = _service(tmp_path)

    html = service.render_page()

    assert 'Save and queue fit' in html
    assert 'test-token' in html
    assert first.path.name in html
    assert second.path.name in html
    assert 'Use automatic' not in html


def test_policy_save_prefills_a_reason_and_queues_without_fitting(tmp_path):
    service, dataset, _, _ = _service(tmp_path)
    target = dataset.targets[0]

    result = service.save_policy({
        'target_key': target.key,
        'scope': 'global',
        'stellar_model': 'kurucz-0.0',
        'disk_model': 'modbb_disk_r',
        'max_disks': 1,
        'search_mode': 'evidence-guided',
        'actor': 'Grant',
        'reason': '',
    })

    record = service.state.policy(target.key)
    assert result['status'] == 'queued'
    assert record.policy.stellar_models == ('kurucz-0.0',)
    assert record.reason.startswith('Use kurucz-0.0 for Test')
    assert service.state.jobs()[0].status == 'queued'


def test_selection_actions_do_not_run_a_fit(tmp_path):
    service, _, first, second = _service(tmp_path)

    shown = service.update_selection({
        'reference': second.path.name, 'action': 'display', 'displayed': True,
    })
    selected = service.update_selection({
        'reference': second.path.name, 'action': 'default',
    })

    assert shown['displayed'] == [first.candidate_id, second.candidate_id]
    assert selected['default'] == second.candidate_id
    assert service.state.jobs() == ()
