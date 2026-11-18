import sqlite3

import pytest

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


def test_explicit_default_prevents_retained_fit_reseeding(tmp_path):
    state = sdf.operator_state.OperatorState(tmp_path / 'operator.sqlite')
    package = tmp_path / 'package'
    package.mkdir()
    state.seed_policy(_policy(), package=package, main_id='Test')

    state.reset_policy_and_queue(
        'simbad:Test', package=package, scope='package', actor='Grant',
        reason='Use the default config',
    )

    assert state.policy('simbad:Test', package=package) is None
    assert not state.seed_policy(_policy(), package=package, main_id='Test')
    assert state.policy('simbad:Test', package=package) is None


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


def test_workers_claim_different_packages_concurrently(tmp_path):
    state = sdf.operator_state.OperatorState(tmp_path / 'operator.sqlite')
    first = tmp_path / 'first'
    second = tmp_path / 'second'
    first.mkdir()
    second.mkdir()
    for package in (first, second):
        state.queue_package(
            package, target_keys=('simbad:{}'.format(package.name),),
            actor='sdf watcher', reason='Fit new data',
        )

    one = state.claim_next(worker_id='worker-1')
    two = state.claim_next(worker_id='worker-2')

    assert {one.package_path, two.package_path} == {
        str(first.resolve()), str(second.resolve()),
    }
    assert {one.worker_id, two.worker_id} == {'worker-1', 'worker-2'}


def test_second_job_for_running_package_waits_for_first(tmp_path):
    state = sdf.operator_state.OperatorState(tmp_path / 'operator.sqlite')
    first = tmp_path / 'first'
    second = tmp_path / 'second'
    first.mkdir()
    second.mkdir()
    state.queue_package(
        first, target_keys=('simbad:A',), actor='Grant', reason='First fit',
    )
    running = state.claim_next(worker_id='worker-1')
    state.save_policy_and_queue(
        _policy('simbad:A'), package=first, actor='Grant',
        reason='Changed policy during fit',
    )
    state.queue_package(
        second, target_keys=('simbad:B',), actor='Grant',
        reason='Independent fit',
    )

    independent = state.claim_next(worker_id='worker-2')
    assert independent.package_path == str(second.resolve())
    assert state.claim_next(worker_id='worker-3') is None

    state.finish(running.id)
    waiting = state.claim_next(worker_id='worker-3')
    assert waiting.package_path == str(first.resolve())


def test_recovery_can_be_limited_to_failed_worker(tmp_path):
    state = sdf.operator_state.OperatorState(tmp_path / 'operator.sqlite')
    for name in ('first', 'second'):
        package = tmp_path / name
        package.mkdir()
        state.queue_package(
            package, target_keys=('simbad:{}'.format(name),),
            actor='sdf watcher', reason='Fit new data',
        )
    state.claim_next(worker_id='worker-1')
    state.claim_next(worker_id='worker-2')

    assert state.recover_running_jobs(worker_id='worker-1') == 1
    assert len(state.jobs(status='queued')) == 1
    running = state.jobs(status='running')
    assert len(running) == 1
    assert running[0].worker_id == 'worker-2'


def test_version_one_operator_database_is_migrated_in_place(tmp_path):
    path = tmp_path / 'operator.sqlite'
    with sqlite3.connect(path) as connection:
        connection.executescript('''
            CREATE TABLE operator_metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );
            INSERT INTO operator_metadata VALUES ('schema_version', '1');
            CREATE TABLE fit_jobs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                package_path TEXT NOT NULL,
                action TEXT NOT NULL,
                payload TEXT NOT NULL,
                status TEXT NOT NULL,
                actor TEXT NOT NULL,
                reason TEXT NOT NULL,
                created_at TEXT NOT NULL,
                started_at TEXT,
                completed_at TEXT,
                error TEXT
            );
        ''')

    sdf.operator_state.OperatorState(path).initialize()

    with sqlite3.connect(path) as connection:
        version = connection.execute(
            "SELECT value FROM operator_metadata WHERE key = 'schema_version'"
        ).fetchone()[0]
        columns = {
            row[1] for row in connection.execute('PRAGMA table_info(fit_jobs)')
        }
    assert version == '7'
    assert 'worker_id' in columns


def test_version_two_operator_database_records_default_choices(tmp_path):
    path = tmp_path / 'operator.sqlite'
    with sqlite3.connect(path) as connection:
        connection.executescript('''
            CREATE TABLE operator_metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );
            INSERT INTO operator_metadata VALUES ('schema_version', '2');
        ''')

    sdf.operator_state.OperatorState(path).initialize()

    with sqlite3.connect(path) as connection:
        version = connection.execute(
            "SELECT value FROM operator_metadata WHERE key = 'schema_version'"
        ).fetchone()[0]
        tables = {
            row[0] for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            )
        }
    assert version == '7'
    assert 'fit_policy_default' in tables


def test_version_five_operator_database_adds_sample_edit_fields(tmp_path):
    path = tmp_path / 'operator.sqlite'
    with sqlite3.connect(path) as connection:
        connection.executescript('''
            CREATE TABLE operator_metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );
            INSERT INTO operator_metadata VALUES ('schema_version', '5');
            CREATE TABLE sample_acquisition_jobs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                sample_name TEXT NOT NULL,
                note TEXT,
                visibility TEXT NOT NULL,
                targets TEXT NOT NULL,
                status TEXT NOT NULL,
                actor TEXT NOT NULL,
                created_at TEXT NOT NULL,
                started_at TEXT,
                completed_at TEXT,
                result TEXT,
                error TEXT,
                worker_id TEXT
            );
            INSERT INTO sample_acquisition_jobs (
                sample_name, note, visibility, targets, status, actor,
                created_at
            ) VALUES ('Existing', NULL, 'private', '["Vega"]', 'queued',
                      'alice', '2026-01-01T00:00:00+00:00');
        ''')

    state = sdf.operator_state.OperatorState(path).initialize()
    job = state.sample_acquisition(1)

    assert job.operation == 'add'
    assert job.owner_username is None
    with sqlite3.connect(path) as connection:
        version = connection.execute(
            "SELECT value FROM operator_metadata WHERE key = 'schema_version'"
        ).fetchone()[0]
    assert version == '7'


def test_acquisition_queue_coalesces_claims_completes_and_recovers(tmp_path):
    state = sdf.operator_state.OperatorState(
        tmp_path / 'operator.sqlite'
    ).initialize()

    first = state.queue_acquisition('  HD   123  ', actor='alice')
    duplicate = state.queue_acquisition('hd 123', actor='alice')
    another_user = state.queue_acquisition('HD 123', actor='bob')

    assert duplicate == first
    assert another_user != first
    claimed = state.claim_next_acquisition(worker_id='import-1')
    assert claimed.query == 'HD 123'
    assert claimed.actor == 'alice'
    assert claimed.status == 'running'
    assert state.claim_next_acquisition().actor == 'bob'

    assert state.recover_running_acquisitions() == 2
    claimed = state.claim_next_acquisition(worker_id='import-2')
    state.finish_acquisition(
        claimed.id,
        result={'sdbids': ['sdb-test'], 'packages': ['sdb-test']},
    )

    completed = state.acquisition(first)
    assert completed.status == 'completed'
    assert completed.result['packages'] == ['sdb-test']
    assert state.acquisitions(actor='alice') == (completed,)


def test_sample_acquisition_queue_validates_claims_and_completes(tmp_path):
    state = sdf.operator_state.OperatorState(
        tmp_path / 'operator.sqlite'
    ).initialize()

    job_id = state.queue_sample_acquisition(
        'Nearby stars', [' HD   123 ', 'Vega', 'HD 123'], actor='alice',
        note='A test sample', visibility='private',
    )
    job = state.claim_next_sample_acquisition(worker_id='import-1')

    assert job.id == job_id
    assert job.sample_name == 'Nearby stars'
    assert job.targets == ('HD 123', 'Vega')
    assert job.note == 'A test sample'
    assert job.operation == 'add'
    assert job.status == 'running'
    assert job.progress_stage == 'starting'
    assert job.progress_current == 0
    assert job.progress_total == 2
    assert state.update_sample_acquisition_progress(job.id, {
        'stage': 'providers',
        'description': 'Updating catalogue providers',
        'current': 1,
        'total': 2,
    })
    progress = state.sample_acquisition(job.id)
    assert progress.progress_stage == 'providers'
    assert progress.progress_message == 'Updating catalogue providers'
    assert progress.progress_current == 1
    state.update_sample_acquisition_progress(job.id, {
        'stage': 'allwise',
        'description': 'Starting AllWISE chunks',
        'current': None,
        'total': None,
    })
    progress = state.sample_acquisition(job.id)
    assert progress.progress_current == 1
    assert progress.progress_total == 2
    assert state.recover_running_sample_acquisitions() == 1

    job = state.claim_next_sample_acquisition(worker_id='import-2')
    state.finish_sample_acquisition(
        job.id, result={'sample': 'Nearby stars', 'packages': ['one']}
    )
    completed = state.sample_acquisition(job_id)
    assert completed.status == 'completed'
    assert completed.result['packages'] == ['one']

    remove_id = state.queue_sample_acquisition(
        'Nearby stars', ['sdb-one'], actor='alice', operation='remove',
        owner_username='alice',
    )
    removal = state.sample_acquisition(remove_id)
    assert removal.operation == 'remove'
    assert removal.owner_username == 'alice'

    with pytest.raises(ValueError, match='1000'):
        state.queue_sample_acquisition(
            'Too large', [str(value) for value in range(1001)], actor='alice'
        )
