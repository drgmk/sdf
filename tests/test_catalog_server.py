import json
from threading import Thread
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import pytest

from .test_dataset import _write_rawphot
import sdf.catalog_server


def _catalog(tmp_path):
    package = tmp_path / 'package'
    package.mkdir()
    _write_rawphot(
        package / 'target-rawphot.txt',
        sdbid='sdb-test', main_id='HD 123',
    )
    (tmp_path / 'export-1-manifest.json').write_text(json.dumps({
        'schema': 'sdb-fit-package-export', 'schema_version': 2,
        'run_id': 1, 'status': 'completed',
        'selection': {
            'kind': 'sample', 'sample': 'test', 'value': 'test',
            'selected_sdbids': ['sdb-test'],
        },
    }), encoding='utf-8')
    return sdf.catalog_server.CatalogService(
        tmp_path, actor='Grant', token='test-token', rebuild=True,
    )


def test_catalog_renders_unfit_search_sample_and_queue(tmp_path):
    service = _catalog(tmp_path)

    assert 'HD 123' in service.render_dashboard()
    assert 'HD 123' in service.render_dashboard('123')
    assert 'test' in service.render_samples()
    assert 'HD 123' in service.render_sample('test')
    assert 'Queue automatic fit' in service.render_package_state('package')

    result = service.queue_package('package', actor='Grant')

    assert result['status'] == 'queued'
    assert service.state.jobs()[0].package_path == str(tmp_path / 'package')


def test_catalog_http_routes_and_rejects_bad_token(tmp_path):
    _catalog(tmp_path)
    try:
        server = sdf.catalog_server.make_catalog_server(
            tmp_path, actor='Grant', token='test-token', port=0,
        )
    except PermissionError:
        pytest.skip('test environment does not permit loopback sockets')
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    root = 'http://127.0.0.1:{}'.format(server.server_address[1])
    try:
        assert 'HD 123' in urlopen(root + '/').read().decode()
        response = urlopen(root + '/target/sdb-test')
        assert response.url.endswith('/system/package/')
        assert 'Queue automatic fit' in response.read().decode()
        request = Request(
            root + '/system/package/api/queue',
            data=json.dumps({'actor': 'Grant'}).encode(),
            headers={
                'Content-Type': 'application/json',
                'X-SDF-Token': 'wrong',
            },
            method='POST',
        )
        try:
            urlopen(request)
            assert False, 'bad token should fail'
        except HTTPError as error:
            assert error.code == 403
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)
