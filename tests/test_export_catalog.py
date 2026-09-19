import json

import sdf.export_catalog


def _manifest(path, run_id, sample, identifiers, status='completed',
              samples=None, kind='sample'):
    payload = {
        'schema': 'sdb-fit-package-export',
        'schema_version': 2,
        'run_id': run_id,
        'status': status,
        'completed_at': '2026-08-15T00:00:00+00:00',
        'selection': {
            'kind': kind, 'sample': sample, 'value': sample,
            'selected_sdbids': identifiers,
        },
    }
    if samples is not None:
        payload['samples'] = samples
    path.write_text(json.dumps(payload), encoding='utf-8')


def test_samples_use_newest_completed_manifest(tmp_path):
    _manifest(tmp_path / 'export-1-manifest.json', 1, 'nearby', ['old'])
    _manifest(tmp_path / 'export-2-manifest.json', 2, 'nearby', ['new'])
    _manifest(
        tmp_path / 'export-3-manifest.json', 3, 'nearby', ['partial'],
        status='partial',
    )

    samples = sdf.export_catalog.ExportCatalog(tmp_path).samples()

    assert len(samples) == 1
    assert samples[0].name == 'nearby'
    assert samples[0].identifiers == ('new',)
    assert samples[0].run_id == 2


def test_samples_use_latest_complete_root_snapshot_for_any_selector(tmp_path):
    _manifest(tmp_path / 'export-1-manifest.json', 1, 'legacy', ['old'])
    _manifest(
        tmp_path / 'export-2-manifest.json', 2, None, ['a', 'b'],
        kind='all', samples=[
            {
                'name': 'reasons', 'sdbids': ['a', 'b'],
                'sample_date': '2026-08-01', 'note': 'Disk sample',
            },
            {'name': 'follow-up', 'sdbids': ['b']},
        ],
    )
    _manifest(
        tmp_path / 'export-3-manifest.json', 3, None, ['partial'],
        status='partial', kind='all',
        samples=[{'name': 'ignored', 'sdbids': ['partial']}],
    )

    samples = sdf.export_catalog.ExportCatalog(tmp_path).samples()

    assert [sample.name for sample in samples] == ['follow-up', 'reasons']
    assert samples[1].identifiers == ('a', 'b')
    assert samples[1].sample_date == '2026-08-01'
    assert samples[1].note == 'Disk sample'
    assert samples[1].run_id == 2
