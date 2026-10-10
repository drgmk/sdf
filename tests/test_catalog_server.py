import json
from types import SimpleNamespace
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
    principal = service.auth.principal({})

    assert 'HD 123' in service.render_dashboard()
    assert 'HD 123' in service.render_public_search('123')
    assert 'test' in service.render_samples(principal)
    assert 'HD 123' in service.render_sample('test', principal)
    assert 'Queue automatic fit' in service.render_package_state(
        'package', admin=True, principal=principal,
    )
    assert 'Queue automatic fit' not in service.render_public_system('package')

    result = service.queue_package('package', actor='Grant')

    assert result['status'] == 'queued'
    assert service.state.jobs()[0].package_path == str(tmp_path / 'package')


def test_catalog_http_routes_and_rejects_bad_token(tmp_path):
    _catalog(tmp_path)
    (tmp_path / 'package' / 'operator.html').write_text(
        'legacy page containing test-token', encoding='utf-8'
    )
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
        assert 'Search fitted systems' in urlopen(root + '/').read().decode()
        assert 'HD 123' in urlopen(root + '/search?q=123').read().decode()
        sample_page = urlopen(root + '/samples/test/').read().decode()
        assert 'Physical target' in sample_page
        assert 'Awaiting a current fit' in sample_page
        assert 'Bokeh' in urlopen(
            root + '/samples/test/plots'
        ).read().decode()
        assert b'VOTABLE' in urlopen(
            root + '/samples/test/sample.xml'
        ).read()
        assert urlopen(
            root + '/samples/test/photometry.csv'
        ).read().startswith(b'package,candidate_id')
        response = urlopen(root + '/target/sdb-test')
        assert response.url.endswith('/system/package/')
        assert 'Queue automatic fit' not in response.read().decode()
        assert 'Queue automatic fit' in urlopen(
            root + '/admin/system/package/'
        ).read().decode()
        for asset in ('operator.html', 'operator%2Ehtml', 'OPERATOR.html'):
            try:
                urlopen(root + '/system/package/' + asset)
                assert False, 'legacy operator page should not be public'
            except HTTPError as error:
                assert error.code == 404
        request = Request(
            root + '/admin/system/package/api/queue',
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


def test_proxy_auth_separates_public_sample_and_admin_access(tmp_path):
    _catalog(tmp_path)
    (tmp_path / 'export-2-manifest.json').write_text(json.dumps({
        'schema': 'sdb-fit-package-export', 'schema_version': 2,
        'run_id': 2, 'status': 'completed',
        'samples': [
            {'name': 'test', 'sdbids': ['sdb-test']},
            {'name': 'private', 'sdbids': ['sdb-test']},
            {'name': 'other', 'sdbids': ['sdb-test']},
        ],
    }), encoding='utf-8')
    groups = tmp_path / '.htgroup'
    groups.write_text(
        'admin: grant\nprivate: alice\nother: bob\n', encoding='utf-8'
    )
    public = tmp_path / 'public-samples'
    public.write_text('test\n', encoding='utf-8')
    try:
        server = sdf.catalog_server.make_catalog_server(
            tmp_path, actor='Local', token='test-token', port=0,
            auth_user_header='X-SDF-User', auth_group_file=groups,
            public_sample_file=public,
        )
    except PermissionError:
        pytest.skip('test environment does not permit loopback sockets')
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    root = 'http://127.0.0.1:{}'.format(server.server_address[1])

    def read(path, user=None, cookie=None):
        headers = {} if user is None else {'X-SDF-User': user}
        if cookie is not None:
            headers['Cookie'] = cookie
        return urlopen(Request(root + path, headers=headers)).read().decode()

    def status(path, user=None):
        try:
            read(path, user)
        except HTTPError as error:
            return error.code
        return 200

    def post_status(user=None):
        headers = {
            'Content-Type': 'application/json',
            'X-SDF-Token': 'test-token',
        }
        if user is not None:
            headers['X-SDF-User'] = user
        request = Request(
            root + '/admin/system/package/api/queue',
            data=json.dumps({'reason': 'Access check'}).encode(),
            headers=headers, method='POST',
        )
        try:
            urlopen(request)
        except HTTPError as error:
            return error.code
        return 200

    try:
        assert 'Search fitted systems' in read('/')
        assert 'No public fit' in read('/system/package/')
        public_page = read('/samples')
        assert 'test' in public_page
        assert 'private' not in public_page
        assert 'Sign in to see shared samples' in public_page
        assert 'HD 123' in read('/samples/test')
        assert status('/samples/private') == 401
        assert status('/admin') == 401
        assert post_status() == 401

        alice_page = read('/samples', 'alice')
        assert 'test' in alice_page
        assert 'private' in alice_page
        assert 'other' not in alice_page
        assert 'HD 123' in read('/samples/test', 'alice')
        assert 'HD 123' in read('/samples/private', 'alice')
        assert status('/admin', 'alice') == 403
        assert post_status('alice') == 403

        assert 'private' not in read('/samples', 'bob')
        assert status('/samples/private', 'bob') == 403

        session = server.catalog_service.auth.issue_session('alice')
        session_page = read('/samples', cookie='sdf_session=' + session)
        assert 'private' in session_page
        assert 'alice' in session_page

        assert 'SDF catalog' in read('/admin', 'grant')
        admin_page = read('/admin/system/package/', 'grant')
        assert 'Queue automatic fit' in admin_page
        assert 'Signed in as <strong>grant</strong>' in admin_page

        request = Request(
            root + '/admin/system/package/api/queue',
            data=json.dumps({
                'actor': 'spoofed', 'reason': 'Authenticated fit',
            }).encode(),
            headers={
                'Content-Type': 'application/json',
                'X-SDF-Token': 'test-token',
                'X-SDF-User': 'grant',
            },
            method='POST',
        )
        assert json.loads(urlopen(request).read())['status'] == 'queued'
        assert server.catalog_service.state.jobs()[0].actor == 'grant'
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


def test_catalog_runtime_starts_requested_isolated_worker_pool(
        tmp_path, monkeypatch):
    processes = []
    threads = []

    class FakeProcess:
        def __init__(self, *, target, args, name, daemon):
            self.target = target
            self.args = args
            self.name = name
            self.daemon = daemon
            self.pid = 1000 + len(processes)
            processes.append(self)

        def start(self):
            pass

        def is_alive(self):
            return True

        def join(self, timeout=None):
            pass

    class FakeThread:
        def __init__(self, *, target, name, daemon):
            self.target = target
            self.name = name
            self.daemon = daemon
            threads.append(self)

        def start(self):
            pass

        def join(self, timeout=None):
            pass

    state = sdf.catalog_server.OperatorState(
        tmp_path / 'operator.sqlite'
    ).initialize()
    monkeypatch.setattr(sdf.catalog_server.multiprocessing, 'Process', FakeProcess)
    monkeypatch.setattr(sdf.catalog_server, 'Thread', FakeThread)
    runtime = sdf.catalog_server.CatalogRuntime(
        SimpleNamespace(state=state), workers=3, watch=False,
    )

    runtime.start()

    assert [process.name for process in processes] == [
        'sdf-fit-worker-1', 'sdf-fit-worker-2', 'sdf-fit-worker-3',
    ]
    assert len({process.args[-1] for process in processes}) == 3
    assert all(process.args[-1] in runtime.worker_ids.values()
               for process in processes)
    assert [thread.name for thread in threads] == ['sdf-worker-supervisor']
