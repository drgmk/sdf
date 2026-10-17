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
