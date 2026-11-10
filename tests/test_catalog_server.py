import json
import re
import sqlite3
from types import SimpleNamespace
from threading import Thread
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import pytest

from .test_dataset import _write_rawphot
import sdf.catalog_server
import sdf.sdb_import


def test_sigterm_handler_uses_normal_keyboard_interrupt_cleanup():
    with pytest.raises(KeyboardInterrupt):
        sdf.catalog_server._raise_keyboard_interrupt(None, None)


def _catalog(tmp_path):
    package = tmp_path / 'package'
    package.mkdir()
    _write_rawphot(
        package / 'target-rawphot.txt',
        sdbid='sdb-test', main_id='HD 123',
    )
    (package / 'identities.yml').write_text(
        'version: 1\ntargets:\n  sdb-test:\n'
        '    aliases:\n      - HD 123\n      - HIP 107000\n',
        encoding='utf-8',
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
    sample_page = service.render_sample('test', principal)
    assert 'HD 123' in sample_page
    assert 'Edit sample' in sample_page
    assert 'Remove targets' in sample_page
    assert re.search(
        r'<label data-member-label="[^"]*hip 107000', sample_page,
    )
    assert '.sample-member-list label[hidden]{display:none}' in (
        sample_page
    )
    assert 'Queue automatic fit' in service.render_package_state(
        'package', admin=True, principal=principal,
    )
    assert 'Queue automatic fit' not in service.render_public_system('package')

    result = service.queue_package('package', actor='Grant')

    assert result['status'] == 'queued'
    assert service.state.jobs()[0].package_path == str(tmp_path / 'package')


def test_public_system_render_is_revision_cached_and_uses_shared_resources(
    tmp_path, monkeypatch,
):
    service = _catalog(tmp_path)
    package_calls = []

    def package(package_path):
        package_calls.append(package_path)
        return {'package_path': 'package', 'status': 'fitted'}

    monkeypatch.setattr(service, 'package', package)
    calls = []

    def render(manager, output, **options):
        calls.append(options)
        path = sdf.catalog_server.Path(output)
        path.write_text('render {}'.format(len(calls)), encoding='utf-8')
        return path

    monkeypatch.setattr(sdf.catalog_server.FitManager, 'render', render)

    assert service.render_public_system('package') == 'render 1'
    assert service.render_public_system('package') == 'render 1'
    assert len(calls) == 1
    assert package_calls == ['package']
    assert calls[0]['inline'] is False

    identities = tmp_path / 'package' / 'identities.yml'
    identities.write_text(
        identities.read_text(encoding='utf-8') + '\n# updated\n',
        encoding='utf-8',
    )

    assert service.render_public_system('package') == 'render 2'
    assert len(calls) == 2


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
        assert 'class="sample-page"' in sample_page
        assert 'Filter targets' in sample_page
        assert 'T<sub>eff</sub>' in sample_page
        assert 'VOTable' in sample_page
        assert 'Photometry' in sample_page
        assert 'class="sample-download-links"' in sample_page
        assert 'Awaiting a current fit' in sample_page
        plot_page = urlopen(root + '/samples/test/plots').read().decode()
        assert 'Bokeh' in plot_page
        assert 'Circles are direct sample members' not in plot_page
        assert '<h1>test</h1>' in plot_page
        assert 'VOTable' in plot_page
        assert 'Photometry' in plot_page
        flux_size_page = urlopen(
            root + '/samples/test/flux-size'
        ).read().decode()
        assert 'VOTable' in flux_size_page
        assert 'Photometry' in flux_size_page
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


def test_signed_in_search_can_queue_and_view_sdb_import(tmp_path):
    _catalog(tmp_path)
    config = sdf.sdb_import.SdbImportConfig.create(
        tmp_path / 'sdb.sqlite'
    )
    try:
        server = sdf.catalog_server.make_catalog_server(
            tmp_path, actor='Grant', token='test-token', port=0,
            sdb_import_config=config,
        )
    except PermissionError:
        pytest.skip('test environment does not permit loopback sockets')
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    root = 'http://127.0.0.1:{}'.format(server.server_address[1])
    try:
        missing = urlopen(root + '/search?q=Vega').read().decode()
        assert 'Import Vega from SIMBAD' in missing
        request = Request(
            root + '/imports',
            data=json.dumps({'query': 'Vega'}).encode(),
            headers={
                'Content-Type': 'application/json',
                'X-SDF-Token': 'test-token',
            },
            method='POST',
        )
        response = urlopen(request)
        assert response.status == 202
        job_id = json.loads(response.read())['job_id']
        page = urlopen(root + '/imports/{}'.format(job_id)).read().decode()
        assert 'Waiting to import' in page
        assert 'Vega' in page

        sample_form = urlopen(root + '/samples/new').read().decode()
        assert 'Create a sample' in sample_form
        assert 'Target list file' in sample_form
        assert 'Check target list' in sample_form
        request = Request(
            root + '/sample-imports',
            data=json.dumps({
                'sample_name': 'Nearby', 'targets': ['HD 123', 'Vega'],
                'note': 'Nearby systems', 'visibility': 'private',
            }).encode(),
            headers={
                'Content-Type': 'application/json',
                'X-SDF-Token': 'test-token',
            }, method='POST',
        )
        response = urlopen(request)
        assert response.status == 202
        sample_job_id = json.loads(response.read())['job_id']
        page = urlopen(
            root + '/sample-imports/{}'.format(sample_job_id)
        ).read().decode()
        assert 'Waiting to import' in page
        assert 'Nearby' in page
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


def test_completed_sample_import_renders_target_results(tmp_path):
    service = _catalog(tmp_path)
    principal = service.auth.principal({})
    job_id = service.state.queue_sample_acquisition(
        'Nearby', ('49 Ceti',), actor=principal.username,
    )
    assert service.state.claim_next_sample_acquisition() is not None
    service.state.finish_sample_acquisition(job_id, result={
        'succeeded_count': 1,
        'failed_count': 0,
        'fit_job_ids': [],
        'manifests': ['export-manifest.json'],
        'import_result': {
            'import_result': {
                'items': [{
                    'requested_name': '49 Ceti',
                    'status': 'existing',
                    'error': None,
                }],
            },
        },
    })

    page = service.render_sample_acquisition(job_id, principal)

    assert 'Sample imported' in page
    assert '49 Ceti' in page
    assert 'Existing target' in page
    assert '/samples/Nearby/' in page


def test_running_sample_import_renders_progress(tmp_path):
    service = _catalog(tmp_path)
    principal = service.auth.principal({})
    job_id = service.state.queue_sample_acquisition(
        'Nearby', ('49 Ceti', 'Vega'), actor=principal.username,
    )
    assert service.state.claim_next_sample_acquisition() is not None
    service.state.update_sample_acquisition_progress(job_id, {
        'stage': 'targets',
        'description': 'Resolving sample identifiers',
        'current': 1,
        'total': 2,
    })

    page = service.render_sample_acquisition(job_id, principal)

    assert 'Resolving sample identifiers' in page
    assert '1 of 2' in page
    assert '<progress value="1" max="2">' in page


def test_completed_sample_removal_renders_membership_results(tmp_path):
    service = _catalog(tmp_path)
    principal = service.auth.principal({})
    job_id = service.state.queue_sample_acquisition(
        'Nearby', ('sdb-test',), actor=principal.username,
        operation='remove',
    )
    assert service.state.claim_next_sample_acquisition() is not None
    service.state.finish_sample_acquisition(job_id, result={
        'succeeded_count': 1,
        'failed_count': 0,
        'fit_job_ids': [],
        'manifests': ['export-manifest.json'],
        'member_results': [{'target': 'sdb-test', 'action_id': 4}],
    })

    page = service.render_sample_acquisition(job_id, principal)

    assert 'Sample updated' in page
    assert 'Membership results' in page
    assert 'sdb-test' in page
    assert 'Removed' in page
    assert '/samples/Nearby/' in page


def test_registering_acquisition_export_indexes_and_queues_fit(tmp_path):
    service = _catalog(tmp_path)
    imported = tmp_path / 'new-package'
    imported.mkdir()
    _write_rawphot(
        imported / 'new-rawphot.txt', sdbid='sdb-new', main_id='Vega'
    )

    result = service.register_acquisition_export({
        'query': 'Vega', 'packages': ['new-package'],
    }, actor='alice')

    assert len(result['fit_job_ids']) == 1
    assert service.package('new-package')['status'] == 'unfit'
    assert service.state.jobs()[0].actor == 'alice'


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
    access_database = tmp_path / 'sdb.sqlite'
    with sqlite3.connect(access_database) as connection:
        connection.executescript('''
            CREATE TABLE users (
                username TEXT PRIMARY KEY,
                is_admin INTEGER NOT NULL,
                disabled INTEGER NOT NULL
            );
            CREATE TABLE samples (
                id INTEGER PRIMARY KEY,
                name TEXT UNIQUE NOT NULL,
                owner_username TEXT,
                visibility TEXT NOT NULL
            );
            CREATE TABLE sample_grants (
                sample_id INTEGER NOT NULL,
                username TEXT NOT NULL,
                PRIMARY KEY (sample_id, username)
            );
            INSERT INTO users VALUES ('grant', 1, 0);
            INSERT INTO users VALUES ('alice', 0, 0);
            INSERT INTO users VALUES ('bob', 0, 0);
            INSERT INTO samples VALUES (1, 'test', NULL, 'public');
            INSERT INTO samples VALUES (2, 'private', 'alice', 'private');
            INSERT INTO samples VALUES (3, 'other', 'bob', 'private');
        ''')
    try:
        server = sdf.catalog_server.make_catalog_server(
            tmp_path, actor='Local', token='test-token', port=0,
            auth_user_header='X-SDF-User',
            auth_database=access_database,
            review_base_url='/review',
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
        assert 'Sign in to see private samples' in public_page
        assert 'HD 123' in read('/samples/test')
        assert status('/samples/private') == 401
        assert status('/admin') == 401
        assert status('/admin/activity') == 401
        assert post_status() == 401

        alice_page = read('/samples', 'alice')
        assert 'test' in alice_page
        assert 'private' in alice_page
        assert 'other' not in alice_page
        assert 'HD 123' in read('/samples/test', 'alice')
        assert 'HD 123' in read('/samples/private', 'alice')
        assert status('/admin', 'alice') == 403
        assert status('/admin/activity', 'alice') == 403
        assert post_status('alice') == 403

        assert 'private' not in read('/samples', 'bob')
        assert status('/samples/private', 'bob') == 403

        collaborator_calls = []

        class CollaboratorImporter:
            def update_sample_collaborator(
                    self, sample_name, username, *, operation, actor):
                collaborator_calls.append(
                    (sample_name, username, operation, actor)
                )
                return {
                    'sample': sample_name, 'username': username,
                    'operation': operation,
                }

        server.catalog_service.importer = CollaboratorImporter()
        with sqlite3.connect(access_database) as connection:
            connection.execute(
                "INSERT INTO sample_grants VALUES (2, 'bob')"
            )
        shared_page = read('/samples/private', 'bob')
        assert 'Edit sample' in shared_page
        assert '>Sharing</button>' in shared_page
        assert '<strong>Co-owners</strong>' in shared_page
        assert 'alice' in shared_page
        assert 'data-unshare-user="bob"' in shared_page

        request = Request(
            root + '/sample-collaborators',
            data=json.dumps({
                'sample_name': 'private', 'operation': 'share',
                'username': 'grant',
            }).encode(),
            headers={
                'Content-Type': 'application/json',
                'X-SDF-Token': 'test-token',
                'X-SDF-User': 'bob',
            },
            method='POST',
        )
        assert json.loads(urlopen(request).read())['username'] == 'grant'
        assert collaborator_calls == [
            ('private', 'grant', 'share', 'bob'),
        ]

        session = server.catalog_service.auth.issue_session('alice')
        session_page = read('/samples', cookie='sdf_session=' + session)
        assert 'private' in session_page
        assert 'alice' in session_page

        dashboard = read('/admin', 'grant')
        assert 'SED catalogue' in dashboard
        assert 'href="/review/">Data review</a>' in dashboard
        assert 'href="/docs">Docs</a>' in dashboard
        assert 'href="/">Home</a>' in dashboard
        assert 'href="/admin/activity">Activity</a>' in dashboard
        activity = read('/admin/activity?days=7', 'grant')
        assert '<h1>Site activity</h1>' in activity
        assert 'What is here' in activity
        assert 'Recent activity' in activity
        admin_page = read('/admin/system/package/', 'grant')
        assert 'Queue automatic fit' in admin_page
        assert 'Signed in as <strong>grant</strong>' in admin_page
        assert 'Queue automatic fit' in read(
            '/admin/target/HD%20123', 'grant'
        )
        fit_status = json.loads(read(
            '/admin/target/HD%20123/api/status', 'grant'
        ))
        assert fit_status['package'] == 'package'
        assert fit_status['job'] is None

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


def test_catalog_runtime_allows_only_one_worker_service(tmp_path, monkeypatch):
    class FakeProcess:
        def __init__(self, *, target, args, name, daemon):
            self.pid = 1

        def start(self):
            pass

        def is_alive(self):
            return True

    class FakeThread:
        def __init__(self, *, target, name, daemon):
            pass

        def start(self):
            pass

    state = sdf.catalog_server.OperatorState(
        tmp_path / 'operator.sqlite'
    ).initialize()
    monkeypatch.setattr(
        sdf.catalog_server.multiprocessing, 'Process', FakeProcess
    )
    monkeypatch.setattr(sdf.catalog_server, 'Thread', FakeThread)
    first = sdf.catalog_server.CatalogRuntime(
        SimpleNamespace(state=state), workers=1, watch=False,
    )
    second = sdf.catalog_server.CatalogRuntime(
        SimpleNamespace(state=state), workers=1, watch=False,
    )

    first.start()
    with pytest.raises(RuntimeError, match='already active'):
        second.start()
    first.worker_lock.release()


def test_catalog_runtime_does_not_lock_around_sdb_import(
        tmp_path, monkeypatch):
    threads = []

    class FakeThread:
        def __init__(self, *, target, args=(), kwargs=None, name, daemon):
            self.target = target
            self.args = args
            self.kwargs = kwargs or {}
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
    service = SimpleNamespace(
        state=state,
        importer=object(),
        register_acquisition_export=lambda result, actor: result,
        register_sample_acquisition_export=lambda result, actor: result,
        export_sync_lock=object(),
    )
    monkeypatch.setattr(sdf.catalog_server, 'Thread', FakeThread)
    runtime = sdf.catalog_server.CatalogRuntime(
        service, worker=False, watch=False,
    )

    runtime.start()

    assert len(threads) == 1
    assert threads[0].target is sdf.sdb_import.run_acquisition_worker
    assert 'work_lock' not in threads[0].kwargs
    assert 'on_sample_export' in threads[0].kwargs
