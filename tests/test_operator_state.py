from .context import sdf
import sdf.model_policy
import sdf.operator_state


def _policy(target_key='simbad:Test', stellar='phoenix_sol'):
    return sdf.model_policy.TargetModelPolicy(
        target_key=target_key,
        stellar_models=(stellar,),
        disk_models=('modbb_disk_r',),
        max_disks=2,
        origin='operator',
    )


def test_policy_is_persistent_and_package_scope_wins(tmp_path):
    state = sdf.operator_state.OperatorState(tmp_path / 'operator.sqlite')
    package = tmp_path / 'package'
    package.mkdir()

    state.save_policy_and_queue(
        _policy(stellar='phoenix_sol'), package=package, main_id='Test',
        scope='global', actor='Grant', reason='Use the normal-star grid',
    )
    state.save_policy_and_queue(
        _policy(stellar='koester_wd'), package=package, main_id='Test',
        scope='package', actor='Grant', reason='This component is a WD',
    )

    record = state.policy('simbad:Test', package=package)
    assert record.policy.stellar_models == ('koester_wd',)
    assert record.reason == 'This component is a WD'
    assert state.policy('simbad:Test').policy.stellar_models == (
        'phoenix_sol',
    )
    assert state.last_actor() == 'Grant'


def test_repeated_policy_edits_coalesce_one_queued_package_job(tmp_path):
    state = sdf.operator_state.OperatorState(tmp_path / 'operator.sqlite')
    package = tmp_path / 'package'
    package.mkdir()

    first = state.save_policy_and_queue(
        _policy('simbad:A'), package=package, actor='Grant',
        reason='Use phoenix_sol',
    )
    second = state.save_policy_and_queue(
        _policy('simbad:B'), package=package, actor='Grant',
        reason='Also update B',
    )

    jobs = state.jobs(package=package)
    assert first == second
    assert len(jobs) == 1
    assert jobs[0].payload['target_keys'] == ['simbad:A', 'simbad:B']
    assert jobs[0].reason == 'Also update B'


def test_fit_job_claim_and_completion_are_explicit(tmp_path):
    state = sdf.operator_state.OperatorState(tmp_path / 'operator.sqlite')
    package = tmp_path / 'package'
    package.mkdir()
    job_id = state.save_policy_and_queue(
        _policy(), package=package, actor='Grant', reason='Try this model',
    )

    claimed = state.claim_next()
    assert claimed.id == job_id
    assert claimed.status == 'running'
    assert state.claim_next() is None

    state.finish(job_id)
    assert state.jobs()[0].status == 'completed'


def test_reset_removes_only_requested_scope_and_queues_refit(tmp_path):
    state = sdf.operator_state.OperatorState(tmp_path / 'operator.sqlite')
    package = tmp_path / 'package'
    package.mkdir()
    state.save_policy_and_queue(
        _policy(), package=package, scope='global', actor='Grant',
        reason='Set global',
    )
    state.save_policy_and_queue(
        _policy(stellar='koester_wd'), package=package, scope='package',
        actor='Grant', reason='Set package exception',
    )
    state.claim_next()

    state.reset_policy_and_queue(
        'simbad:Test', package=package, scope='package', actor='Grant',
        reason='Return to global policy',
    )

    assert state.policy('simbad:Test', package=package).policy.stellar_models \
        == ('phoenix_sol',)
    assert state.jobs(status='queued')[0].payload == {
        'target_keys': ['simbad:Test'],
    }


def test_seed_recovers_retained_policy_without_queue_or_overwrite(tmp_path):
    state = sdf.operator_state.OperatorState(tmp_path / 'operator.sqlite')
    package = tmp_path / 'package'
    package.mkdir()

    assert state.seed_policy(_policy(), package=package, main_id='Test')
    assert not state.seed_policy(
        _policy(stellar='koester_wd'), package=package, main_id='Test'
    )

    record = state.policy('simbad:Test', package=package)
    assert record.policy.stellar_models == ('phoenix_sol',)
    assert record.scope == str(package.resolve())
    assert state.jobs() == ()


def test_interrupted_running_job_is_recovered_for_restart(tmp_path):
    state = sdf.operator_state.OperatorState(tmp_path / 'operator.sqlite')
    package = tmp_path / 'package'
    package.mkdir()
    state.queue_package(
        package, target_keys=('simbad:Test',), actor='sdf watcher',
        reason='Fit changed data',
    )
    state.claim_next()

    assert state.recover_running_jobs() == 1

    job = state.jobs()[0]
    assert job.status == 'queued'
    assert job.started_at is None
    assert job.error == 'worker stopped before completion'
