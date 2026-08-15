import json

import sdf.export_catalog


def _manifest(path, run_id, sample, identifiers, status='completed'):
    path.write_text(json.dumps({
        'schema': 'sdb-fit-package-export',
        'schema_version': 2,
        'run_id': run_id,
        'status': status,
        'completed_at': '2026-08-15T00:00:00+00:00',
        'selection': {
            'kind': 'sample', 'sample': sample, 'value': sample,
            'selected_sdbids': identifiers,
        },
    }), encoding='utf-8')


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
