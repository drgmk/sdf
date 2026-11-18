from .context import sdf
import sdf.operator_jobs
import sdf.operator_state


def _queued_state(tmp_path):
    state = sdf.operator_state.OperatorState(tmp_path / 'operator.sqlite')
    package = tmp_path / 'package'
    package.mkdir()
    state.save_policy_and_queue(
        sdf.model_policy.TargetModelPolicy(
            target_key='simbad:Test', stellar_models=('phoenix_sol',),
        ),
        package=package, actor='Grant', reason='Try the normal-star model',
    )
    return state, package


def test_worker_runs_fit_outside_request_and_completes_job(
        tmp_path, monkeypatch):
    state, package = _queued_state(tmp_path)
    calls = []

    monkeypatch.setattr(
        sdf.operator_jobs, 'fit_package',
        lambda path, **kwargs: calls.append(('fit', path, kwargs)),
    )

    class FakeManager:
        def __init__(self, path):
            calls.append(('manager', path))

        def render(self):
            calls.append(('render',))

    monkeypatch.setattr(sdf.operator_jobs, 'FitManager', FakeManager)

    job = sdf.operator_jobs.run_next_job(state, worker_id='worker-1')

    assert job.package_path == str(package.resolve())
    assert calls[0] == (
        'fit', str(package.resolve()), {'sampler': None, 'policy_store': state},
    )
    assert state.jobs()[0].status == 'completed'
    assert state.jobs()[0].worker_id == 'worker-1'


def test_worker_records_failure_without_losing_job(tmp_path, monkeypatch):
    state, _ = _queued_state(tmp_path)

    def fail(*args, **kwargs):
        raise RuntimeError('sampler stopped')

    monkeypatch.setattr(sdf.operator_jobs, 'fit_package', fail)

    sdf.operator_jobs.run_next_job(state)

    job = state.jobs()[0]
    assert job.status == 'failed'
    assert job.error == 'RuntimeError: sampler stopped'
