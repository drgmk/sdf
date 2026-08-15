"""Catalog-wide local HTTP service over an SDB export root."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict
from functools import partial
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import mimetypes
import multiprocessing
from pathlib import Path, PurePosixPath
import secrets
from threading import Event, Thread
from urllib.parse import parse_qs, quote, unquote, urlparse

import jinja2

from .catalog_sync import CatalogSynchronizer
from .dataset import FitDataset
from .export_catalog import ExportCatalog
from .operator_server import OperatorService
from .operator_jobs import run_worker
from .operator_state import OperatorState
from .result_index import ResultIndex


class CatalogService:
    """Application operations shared by the catalog HTTP boundary."""

    def __init__(self, root, *, index_path=None, state_path=None, actor=None,
                 token=None, rebuild=False):
        self.root = Path(root).expanduser().resolve()
        if not self.root.is_dir():
            raise ValueError('catalog root is not a directory: {}'.format(
                self.root
            ))
        self.index = ResultIndex(
            index_path or self.root / 'sdf-results.sqlite'
        )
        if rebuild or not self.index.path.is_file():
            self.index.rebuild(self.root)
        else:
            try:
                if self.index.root() != self.root:
                    raise ValueError('result index belongs to another root')
                self.index.summary()
                if self.index.verify(self.root):
                    raise ValueError('result index is out of sync')
            except (OSError, ValueError):
                self.index.rebuild(self.root)
        self.state = OperatorState(
            state_path or self.root / 'sdf-operator.sqlite'
        ).initialize()
        self.actor = actor or self.state.last_actor() or ''
        self.token = token or secrets.token_urlsafe(24)
        self.exports = ExportCatalog(self.root)
        self._operators = {}

    def packages(self, *, query=None, identifiers=None):
        rows = self.index.catalog_packages(
            query=query, identifiers=identifiers
        )
        jobs = {}
        for job in reversed(self.state.jobs()):
            try:
                package = Path(job.package_path).relative_to(self.root).as_posix()
            except ValueError:
                continue
            jobs[package] = job
        return tuple(dict(row, job=jobs.get(row['package_path'])) for row in rows)

    def summary(self):
        packages = self.packages()
        return {
            'index': asdict(self.index.summary()),
            'status': dict(Counter(row['status'] for row in packages)),
            'jobs': dict(Counter(job.status for job in self.state.jobs())),
        }

    def package(self, package_path):
        clean = _safe_relative(package_path)
        rows = self.packages()
        for row in rows:
            if row['package_path'] == clean:
                return row
        raise KeyError('exported package not found: {}'.format(clean))

    def package_directory(self, package_path):
        row = self.package(package_path)
        return self.root / row['package_path']

    def operator(self, package_path):
        directory = self.package_directory(package_path)
        key = str(directory)
        service = self._operators.get(key)
        if service is None:
            service = OperatorService(
                directory, self.state, actor=self.actor, token=self.token,
            )
            self._operators[key] = service
        return service

    def queue_package(self, package_path, *, actor=None, reason=None,
                      remember_actor=True):
        directory = self.package_directory(package_path)
        dataset = FitDataset.from_path(directory, nospec=True)
        if dataset is None:
            raise ValueError('{} has no usable observations'.format(directory))
        actor = str(actor or self.actor).strip()
        if not actor:
            raise ValueError('actor is required')
        labels = [
            target.main_id or target.sdbid or target.key
            for target in dataset.targets
        ]
        reason = str(reason or '').strip() or 'Fit {}'.format(
            ' + '.join(labels)
        )
        job_id = self.state.queue_package(
            directory, actor=actor, reason=reason,
            target_keys=(target.key for target in dataset.targets),
            remember_actor=remember_actor,
        )
        if remember_actor:
            self.actor = actor
        return {'job_id': job_id, 'status': 'queued'}

    def render_dashboard(self, query=None, *, status=None, page=1,
                         per_page=100):
        packages = self.packages(query=query)
        if status:
            packages = tuple(
                row for row in packages if row['status'] == status
            )
        total = len(packages)
        page_count = max(1, (total + per_page - 1) // per_page)
        page = min(max(1, int(page)), page_count)
        start = (page - 1) * per_page
        return _render(
            'catalog_dashboard.html', title='SDF catalog',
            packages=packages[start:start + per_page], query=query or '',
            summary=self.summary(), samples=self.exports.samples(),
            status_filter=status or '', total=total, page=page,
            page_count=page_count,
        )

    def render_samples(self):
        return _render(
            'catalog_samples.html', title='Samples',
            samples=self.exports.samples(),
        )

    def render_sample(self, name):
        sample = self.exports.sample(name)
        packages = self.packages(identifiers=sample.identifiers)
        matched = {
            identifier for row in packages for identifier in row['sdbids']
        }
        return _render(
            'catalog_sample.html', title=sample.name, sample=sample,
            packages=packages,
            unmatched=tuple(
                value for value in sample.identifiers if value not in matched
            ),
        )

    def render_jobs(self):
        jobs = []
        for job in self.state.jobs():
            value = asdict(job)
            try:
                value['package'] = Path(job.package_path).relative_to(
                    self.root
                ).as_posix()
            except ValueError:
                value['package'] = job.package_path
            jobs.append(value)
        return _render('catalog_jobs.html', title='Fit jobs', jobs=jobs)

    def render_package_state(self, package_path):
        row = self.package(package_path)
        return _render(
            'catalog_package.html',
            title=' / '.join(row['target_labels']) or row['package_path'],
            package=row, token=self.token, actor=self.actor,
        )


def make_catalog_server(root, *, index_path=None, state_path=None, actor=None,
                        token=None, rebuild=False, host='127.0.0.1', port=0):
    if host not in {'127.0.0.1', 'localhost', '::1'}:
        raise ValueError('catalog service must bind to a loopback address')
    service = CatalogService(
        root, index_path=index_path, state_path=state_path, actor=actor,
        token=token, rebuild=rebuild,
    )
    handler = partial(_CatalogHandler, service=service)
    server = ThreadingHTTPServer((host, port), handler)
    server.catalog_service = service
    return server


def serve_catalog(root, *, index_path=None, state_path=None, actor=None,
                  rebuild=False, host='127.0.0.1', port=8765,
                  worker=True, watch=True, autofit_existing=False,
                  poll_seconds=5.0):
    server = make_catalog_server(
        root, index_path=index_path, state_path=state_path, actor=actor,
        rebuild=rebuild, host=host, port=port,
    )
    address, actual_port = server.server_address[:2]
    print('SDF catalog: http://{}:{}/'.format(address, actual_port))
    print('Result index: {}'.format(server.catalog_service.index.path))
    print('Operator state: {}'.format(server.catalog_service.state.path))
    runtime = CatalogRuntime(
        server.catalog_service, worker=worker, watch=watch,
        autofit_existing=autofit_existing, poll_seconds=poll_seconds,
    )
    runtime.start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        runtime.stop()


class CatalogRuntime:
    """Supervise the optional watcher thread and isolated fit worker."""

    def __init__(self, service, *, worker=True, watch=True,
                 autofit_existing=False, poll_seconds=5.0):
        self.service = service
        self.worker_enabled = worker
        self.watch_enabled = watch
        self.autofit_existing = autofit_existing
        self.poll_seconds = poll_seconds
        self.stop_event = multiprocessing.Event()
        self.watch_stop = Event()
        self.worker_process = None
        self.watch_thread = None
        self.supervisor_thread = None

    def start(self):
        if self.worker_enabled:
            self._start_worker()
            self.supervisor_thread = Thread(
                target=self._supervise_worker,
                name='sdf-worker-supervisor', daemon=True,
            )
            self.supervisor_thread.start()
        if self.watch_enabled:
            synchronizer = CatalogSynchronizer(self.service)
            self.watch_thread = Thread(
                target=synchronizer.run,
                args=(self.watch_stop,),
                kwargs={
                    'poll_seconds': self.poll_seconds,
                    'autofit_existing': self.autofit_existing,
                },
                name='sdf-export-watcher', daemon=True,
            )
            self.watch_thread.start()
            print('Export watcher: every {:.1f}s{}'.format(
                self.poll_seconds,
                ' (queueing existing work)' if self.autofit_existing else '',
            ))

    def stop(self):
        self.watch_stop.set()
        self.stop_event.set()
        if self.watch_thread is not None:
            self.watch_thread.join(timeout=self.poll_seconds + 1)
        if self.supervisor_thread is not None:
            self.supervisor_thread.join(timeout=self.poll_seconds + 1)
        if self.worker_process is not None:
            self.worker_process.join(timeout=self.poll_seconds + 1)
            if self.worker_process.is_alive():
                self.worker_process.terminate()
                self.worker_process.join(timeout=2)

    def _start_worker(self):
        self.worker_process = multiprocessing.Process(
            target=_worker_main,
            args=(str(self.service.state.path), self.stop_event,
                  self.poll_seconds),
            name='sdf-fit-worker', daemon=False,
        )
        self.worker_process.start()
        print('Fit worker: pid {}'.format(self.worker_process.pid))

    def _supervise_worker(self):
        while not self.watch_stop.wait(self.poll_seconds):
            if self.stop_event.is_set():
                return
            if self.worker_process is not None \
                    and not self.worker_process.is_alive():
                print('Fit worker stopped; restarting')
                self._start_worker()


def _worker_main(state_path, stop_event, poll_seconds):
    run_worker(
        OperatorState(state_path).initialize(),
        poll_seconds=poll_seconds, stop_event=stop_event,
    )


class _CatalogHandler(BaseHTTPRequestHandler):
    def __init__(self, *args, service, **kwargs):
        self.service = service
        super().__init__(*args, **kwargs)

    def do_GET(self):
        parsed = urlparse(self.path)
        try:
            if parsed.path == '/':
                values = parse_qs(parsed.query)
                query = values.get('q', [''])[0].strip()
                status = values.get('status', [''])[0].strip()
                page = int(values.get('page', ['1'])[0])
                self._html(self.service.render_dashboard(
                    query or None, status=status or None, page=page,
                ))
                return
            if parsed.path == '/search':
                values = parse_qs(parsed.query)
                query = values.get('q', [''])[0].strip()
                status = values.get('status', [''])[0].strip()
                page = int(values.get('page', ['1'])[0])
                self._html(self.service.render_dashboard(
                    query or None, status=status or None, page=page,
                ))
                return
            if parsed.path == '/samples':
                self._html(self.service.render_samples())
                return
            if parsed.path.startswith('/samples/'):
                name = unquote(parsed.path.removeprefix('/samples/')).strip('/')
                self._html(self.service.render_sample(name))
                return
            if parsed.path == '/jobs':
                self._html(self.service.render_jobs())
                return
            if parsed.path.startswith('/target/'):
                identifier = unquote(parsed.path.removeprefix('/target/'))
                rows = self.service.packages(identifiers=(identifier,))
                if len(rows) == 1:
                    self.send_response(HTTPStatus.SEE_OTHER)
                    self.send_header(
                        'Location', _system_url(rows[0]['package_path'])
                    )
                    self.end_headers()
                    return
                self._html(self.service.render_dashboard(identifier),
                           status=HTTPStatus.NOT_FOUND if not rows else HTTPStatus.OK)
                return
            if parsed.path.startswith('/system/'):
                package, asset = _system_parts(parsed.path)
                row = self.service.package(package)
                if not asset:
                    if row['status'] in {'unfit', 'needs selection'}:
                        self._html(self.service.render_package_state(package))
                    else:
                        prefix = _system_url(package).rstrip('/') + '/api'
                        self._html(self.service.operator(package).render_page(
                            operator_api=prefix
                        ))
                    return
                self._static(self.service.package_directory(package), asset)
                return
            if parsed.path == '/favicon.ico':
                self._send(HTTPStatus.NO_CONTENT, b'', 'image/x-icon')
                return
            self._error(KeyError('page not found'), HTTPStatus.NOT_FOUND)
        except (KeyError, ValueError, OSError) as error:
            self._error(error, HTTPStatus.NOT_FOUND)

    def do_POST(self):
        parsed = urlparse(self.path)
        try:
            package, asset = _system_parts(parsed.path)
            if not asset.startswith('api/'):
                raise KeyError('unknown catalog endpoint')
            if self.headers.get('X-SDF-Token') != self.service.token:
                raise PermissionError('invalid operator token')
            payload = self._payload()
            action = asset.removeprefix('api/')
            if action == 'queue':
                result = self.service.queue_package(
                    package, actor=payload.get('actor'),
                    reason=payload.get('reason'),
                )
            else:
                operator = self.service.operator(package)
                operations = {
                    'policy': operator.save_policy,
                    'system-policy': operator.save_system,
                    'policy/reset': operator.reset_policy,
                    'selection': operator.update_selection,
                }
                if action not in operations:
                    raise KeyError('unknown operator endpoint')
                result = operations[action](payload)
            self._json(result)
        except PermissionError as error:
            self._error(error, HTTPStatus.FORBIDDEN, json_response=True)
        except (KeyError, ValueError, OSError, json.JSONDecodeError) as error:
            self._error(error, HTTPStatus.BAD_REQUEST, json_response=True)

    def _payload(self):
        length = int(self.headers.get('Content-Length', '0'))
        if length <= 0 or length > 65536:
            raise ValueError('invalid request size')
        return json.loads(self.rfile.read(length))

    def _static(self, directory, asset):
        relative = _safe_relative(asset)
        path = (directory / relative).resolve()
        if directory.resolve() not in path.parents or not path.is_file():
            raise KeyError('package asset not found')
        content_type = mimetypes.guess_type(path.name)[0] or \
            'application/octet-stream'
        self._send(HTTPStatus.OK, path.read_bytes(), content_type)

    def _html(self, value, status=HTTPStatus.OK):
        self._send(status, value.encode('utf-8'), 'text/html; charset=utf-8')

    def _json(self, value, status=HTTPStatus.OK):
        self._send(
            status, json.dumps(value).encode('utf-8'),
            'application/json; charset=utf-8',
        )

    def _error(self, error, status, json_response=False):
        if json_response:
            self._json({'error': str(error)}, status=status)
            return
        self._html(_render(
            'catalog_message.html', title='SDF catalog',
            heading=HTTPStatus(status).phrase, message=str(error),
        ), status=status)

    def _send(self, status, body, content_type):
        self.send_response(status)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(body)


def _render(template_name, **values):
    environment = jinja2.Environment(
        autoescape=True,
        loader=jinja2.PackageLoader('sdf', package_path='www/templates'),
    )
    return environment.get_template(template_name).render(**values)


def _safe_relative(value):
    path = PurePosixPath(unquote(str(value)))
    if path.is_absolute() or '..' in path.parts or not path.parts:
        raise ValueError('invalid catalog path')
    return path.as_posix()


def _system_url(package):
    return '/system/{}/'.format(quote(package, safe='/'))


def _system_parts(path):
    relative = path.removeprefix('/system/').lstrip('/')
    if not relative:
        raise KeyError('system package is required')
    parts = relative.split('/')
    # Export packages are currently one directory below the root. Preserve a
    # clear boundary here; nested package roots can later use an encoded key.
    package = _safe_relative(parts[0])
    asset = '/'.join(parts[1:]).strip('/')
    return package, asset
