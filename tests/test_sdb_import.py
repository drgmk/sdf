import json
from types import SimpleNamespace
from threading import Event

from .context import sdf
import sdf.operator_state
import sdf.sdb_import


def test_sdb_importer_ensures_then_exports_exact_target(tmp_path):
    export_root = tmp_path / 'export'
    export_root.mkdir()
    manifest = export_root / 'export-7-manifest.json'
    calls = []

    def runner(command, **kwargs):
        calls.append((command, kwargs))
        if 'add' in command:
            payload = {
                'created_count': 1, 'existing_count': 0,
                'providers': ['simbad', 'gaia_dr3'],
                'items': [{
                    'status': 'created', 'sdbid': 'sdb-test',
                }],
                'update_summary': {'failed': 1},
            }
            return SimpleNamespace(
                returncode=1, stdout=json.dumps(payload),
                stderr='one optional provider failed',
            )
        manifest.write_text(json.dumps({
            'packages': [{
                'directory': 'sdb-test', 'status': 'completed',
            }],
        }), encoding='utf-8')
        payload = {
            'manifest': str(manifest), 'failed': 0,
            'exported': 1, 'skipped': 0,
        }
        return SimpleNamespace(
            returncode=0, stdout=json.dumps(payload), stderr='',
        )

    config = sdf.sdb_import.SdbImportConfig.create(
        tmp_path / 'sdb.sqlite', executable='/opt/sdb',
        reference_database=tmp_path / 'reference.sqlite',
        cache_database=tmp_path / 'cache.sqlite',
        config=tmp_path / 'sdb.toml',
    )
    result = sdf.sdb_import.SdbImporter(
        config, export_root, runner=runner
    ).acquire('HD 123', actor='alice')

    assert result['sdbids'] == ['sdb-test']
    assert result['packages'] == ['sdb-test']
    assert result['warnings']
    assert calls[0][0][-3:] == ['add', 'HD 123', '--ensure']
    assert calls[1][0][-6:] == [
        'export', 'sdb-test', '--output-dir', str(export_root),
        '--workers', '1',
    ]
    assert calls[0][1]['env']['SDB_ACTOR'] == 'alice'
    assert all(call[1]['timeout'] == 1800.0 for call in calls)


def test_acquisition_worker_records_result_and_registration(tmp_path):
    state = sdf.operator_state.OperatorState(
        tmp_path / 'operator.sqlite'
    ).initialize()
    job_id = state.queue_acquisition('Vega', actor='alice')
    stopped = Event()

    class Importer:
        def acquire(self, query, *, actor):
            assert (query, actor) == ('Vega', 'alice')
            return {'query': query, 'packages': ['vega']}

    def register(result, *, actor):
        stopped.set()
        return dict(result, fit_job_ids=[42], actor=actor)

    sdf.sdb_import.run_acquisition_worker(
        state, Importer(), on_export=register,
        poll_seconds=0.01, stop_event=stopped,
    )

    job = state.acquisition(job_id)
    assert job.status == 'completed'
    assert job.result['fit_job_ids'] == [42]
    assert job.result['actor'] == 'alice'


def test_sdb_importer_acquires_and_exports_sample(tmp_path):
    export_root = tmp_path / 'export'
    export_root.mkdir()
    manifest = export_root / 'export-8-manifest.json'
    calls = []

    def runner(command, **kwargs):
        calls.append(command)
        if 'acquire' in command:
            target_file = command[command.index('acquire') + 2]
            assert open(target_file, encoding='utf-8').read() == 'HD 123\nVega\n'
            return SimpleNamespace(
                returncode=1,
                stdout=json.dumps({
                    'requested_count': 2, 'succeeded_count': 1,
                    'failed_count': 1,
                }), stderr='one target failed',
            )
        manifest.write_text(json.dumps({
            'packages': [{
                'directory': 'sdb-test', 'status': 'completed',
            }],
        }), encoding='utf-8')
        return SimpleNamespace(
            returncode=0,
            stdout=json.dumps({
                'manifest': str(manifest), 'failed': 0,
                'exported': 1, 'skipped': 0,
            }), stderr='',
        )

    importer = sdf.sdb_import.SdbImporter(
        sdf.sdb_import.SdbImportConfig.create(tmp_path / 'sdb.sqlite'),
        export_root, runner=runner,
    )
    result = importer.acquire_sample(SimpleNamespace(
        sample_name='Nearby', targets=('HD 123', 'Vega'), actor='alice',
        owner_username='alice', visibility='private', note='Nearby systems',
    ))

    assert result['sample'] == 'Nearby'
    assert result['packages'] == ['sdb-test']
    assert result['failed_count'] == 1
    assert result['warnings']
    assert calls[0][calls[0].index('sample'):] [:2] == ['sample', 'acquire']
    assert calls[0][calls[0].index('--owner') + 1] == 'alice'
    assert calls[1][-7:] == [
        'export', '--sample', 'Nearby', '--output-dir', str(export_root),
        '--workers', '1',
    ]


def test_sdb_importer_checks_sample_targets_without_importing(tmp_path):
    calls = []

    def runner(command, **kwargs):
        calls.append(command)
        target_file = command[-1]
        assert open(target_file, encoding='utf-8').read() == '49 Cet\nNew\n'
        return SimpleNamespace(
            returncode=0,
            stdout=json.dumps({
                'requested_count': 2, 'existing_count': 1,
                'new_count': 1, 'ambiguous_count': 0,
                'items': [],
            }), stderr='',
        )

    importer = sdf.sdb_import.SdbImporter(
        sdf.sdb_import.SdbImportConfig.create(tmp_path / 'sdb.sqlite'),
        tmp_path, runner=runner,
    )
    result = importer.check_sample_targets(('49 Cet', 'New', '49 Cet'))

    assert result['existing_count'] == 1
    assert calls[0][-3:-1] == ['sample', 'check']


def test_sdb_importer_updates_sample_collaborator(tmp_path):
    calls = []

    def runner(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(
            returncode=0,
            stdout=json.dumps({
                'sample_id': 7, 'username': 'bob@example.com',
            }),
            stderr='',
        )

    importer = sdf.sdb_import.SdbImporter(
        sdf.sdb_import.SdbImportConfig.create(tmp_path / 'sdb.sqlite'),
        tmp_path, runner=runner,
    )
    result = importer.update_sample_collaborator(
        'Nearby', 'bob@example.com', operation='share', actor='alice@example.com',
    )

    assert result['username'] == 'bob@example.com'
    assert calls[0][0][-4:] == [
        'sample', 'share', 'Nearby', 'bob@example.com',
    ]
    assert calls[0][1]['env']['SDB_ACTOR'] == 'alice@example.com'


def test_sdb_importer_removes_members_and_exports_sample(tmp_path):
    export_root = tmp_path / 'export'
    export_root.mkdir()
    manifest = export_root / 'export-9-manifest.json'
    calls = []

    def runner(command, **kwargs):
        calls.append(command)
        if 'remove' in command:
            return SimpleNamespace(
                returncode=0,
                stdout=json.dumps({'action_id': len(calls)}), stderr='',
            )
        if 'delete-if-empty' in command:
            return SimpleNamespace(
                returncode=0,
                stdout=json.dumps({
                    'deleted': False, 'member_count': 3, 'deletion': None,
                }), stderr='',
            )
        manifest.write_text(json.dumps({'packages': []}), encoding='utf-8')
        return SimpleNamespace(
            returncode=0,
            stdout=json.dumps({
                'manifest': str(manifest), 'failed': 0,
                'exported': 0, 'skipped': 0,
            }), stderr='',
        )

    importer = sdf.sdb_import.SdbImporter(
        sdf.sdb_import.SdbImportConfig.create(tmp_path / 'sdb.sqlite'),
        export_root, runner=runner,
    )
    result = importer.acquire_sample(SimpleNamespace(
        sample_name='Nearby', targets=('sdb-one', 'sdb-two'), actor='alice',
        operation='remove', visibility='private', note=None,
    ))

    assert result['operation'] == 'remove'
    assert result['succeeded_count'] == 2
    assert result['failed_count'] == 0
    assert result['packages'] == []
    assert calls[0][calls[0].index('sample'):] [:3] == [
        'sample', 'remove', 'Nearby',
    ]
    assert calls[-1][-7:] == [
        'export', '--sample', 'Nearby', '--output-dir', str(export_root),
        '--workers', '1',
    ]
    assert any('delete-if-empty' in call for call in calls)


def test_sdb_importer_deletes_sample_after_final_member(tmp_path):
    export_root = tmp_path / 'export'
    export_root.mkdir()
    calls = []

    def runner(command, **kwargs):
        calls.append(command)
        if 'remove' in command:
            return SimpleNamespace(
                returncode=0,
                stdout=json.dumps({'action_id': 1}), stderr='',
            )
        if 'delete-if-empty' in command:
            return SimpleNamespace(
                returncode=0,
                stdout=json.dumps({
                    'deleted': True, 'member_count': 0,
                    'deletion': {'name': 'Solo'},
                }), stderr='',
            )
        assert 'catalog' in command
        sample_catalog = export_root / 'samples.json'
        sample_catalog.write_text('{}', encoding='utf-8')
        return SimpleNamespace(
            returncode=0,
            stdout=json.dumps({
                'path': str(sample_catalog), 'sample_count': 0,
            }), stderr='',
        )

    importer = sdf.sdb_import.SdbImporter(
        sdf.sdb_import.SdbImportConfig.create(tmp_path / 'sdb.sqlite'),
        export_root, runner=runner,
    )
    result = importer.acquire_sample(SimpleNamespace(
        sample_name='Solo', targets=('sdb-one',), actor='alice',
        operation='remove', visibility='private', note=None,
    ))

    assert result['sample_deleted'] is True
    assert result['remaining_member_count'] == 0
    assert result['sample_catalog'] == 'samples.json'
    assert not any('export' in call for call in calls)
    assert calls[-1][calls[-1].index('sample'):] == [
        'sample', 'catalog', '--output-dir', str(export_root),
    ]


def test_sdb_importer_retains_all_failed_sample_result_without_export(tmp_path):
    calls = []

    def runner(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(
            returncode=1,
            stdout=json.dumps({
                'requested_count': 1, 'succeeded_count': 0,
                'failed_count': 1,
                'import_result': {
                    'items': [{
                        'requested_name': 'New', 'status': 'failed',
                        'error': 'SIMBAD unavailable',
                    }],
                    'update_summary': None,
                },
            }), stderr='SIMBAD unavailable',
        )

    importer = sdf.sdb_import.SdbImporter(
        sdf.sdb_import.SdbImportConfig.create(tmp_path / 'sdb.sqlite'),
        tmp_path, runner=runner,
    )
    result = importer.acquire_sample(SimpleNamespace(
        sample_name='New sample', targets=('New',), actor='alice',
        owner_username='alice', visibility='private', note=None,
    ))

    assert result['succeeded_count'] == 0
    assert result['failed_count'] == 1
    assert result['packages'] == []
    assert len(calls) == 1


def test_acquisition_worker_drains_sample_queue(tmp_path):
    state = sdf.operator_state.OperatorState(
        tmp_path / 'operator.sqlite'
    ).initialize()
    job_id = state.queue_sample_acquisition(
        'Nearby', ['Vega'], actor='alice'
    )
    stopped = Event()

    class Importer:
        def acquire_sample(self, job):
            assert job.sample_name == 'Nearby'
            return {'sample': job.sample_name, 'packages': ['vega']}

    def register(result, *, actor):
        stopped.set()
        return dict(result, fit_job_ids=[43], actor=actor)

    sdf.sdb_import.run_acquisition_worker(
        state, Importer(), on_sample_export=register,
        poll_seconds=0.01, stop_event=stopped,
    )

    job = state.sample_acquisition(job_id)
    assert job.status == 'completed'
    assert job.result['fit_job_ids'] == [43]
    assert job.result['actor'] == 'alice'
