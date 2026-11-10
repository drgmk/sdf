"""Catalog-wide local HTTP service over an sdb export root."""

from __future__ import annotations

from collections import Counter, OrderedDict
from dataclasses import asdict
from functools import partial
import getpass
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import io
import json
import mimetypes
import multiprocessing
import os
from pathlib import Path, PurePosixPath
import secrets
import signal
import tempfile
from threading import Event, RLock, Thread
from urllib.parse import parse_qs, quote, unquote, urlparse

from bokeh.resources import CDN
import jinja2

from .catalog_sync import CatalogSynchronizer
from .catalog_auth import (
    AuthenticationRequired, AuthorizationDenied, CatalogAuthorizer,
)
from .dataset import FitDataset
from .export_catalog import ExportCatalog
from .external_services import ExternalServices
from .fit_management import FitManager
from .identity_manifest import read_identity_manifest
from .locking import InterProcessLock
from .operator_server import OperatorService
from .operator_jobs import run_worker
from .operator_state import OperatorState
from .result_index import ResultIndex
from .sample_diagnostics import (
    PHOTOMETRY_FIELDS, flux_size_plot_components, photometry_rows,
)
from .sample_products import (
    project_sample, sample_plot_components, sample_table_rows,
    sample_votable, summary_text,
)
from .sdb_import import SdbImporter, run_acquisition_worker
from .site_activity import activity_report


class CatalogService:
    """Application operations shared by the catalog HTTP boundary."""

    def __init__(self, root, *, index_path=None, state_path=None, actor=None,
                 token=None, rebuild=False, auth_user_header=None,
                 auth_database=None, auth_session_secret=None,
                 auth_cookie_secure=False,
                 sdb_import_config=None, review_base_url=None,
                 access_log=None):
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
        self.actor = actor or self.state.last_actor() or getpass.getuser()
        self.token = token or secrets.token_urlsafe(24)
        self.exports = ExportCatalog(self.root)
        self.auth = CatalogAuthorizer(
            user_header=auth_user_header, access_database=auth_database,
            local_actor=self.actor, session_secret=auth_session_secret,
        )
        self.external_services = ExternalServices.from_sdb_config(
            None if sdb_import_config is None else sdb_import_config.config
        )
        self.auth_cookie_secure = bool(auth_cookie_secure)
        self.review_base_url = (
            str(review_base_url).rstrip('/') if review_base_url else None
        )
        self.access_log = (
            None if access_log is None
            else Path(access_log).expanduser().resolve()
        )
        self.importer = (
            None if sdb_import_config is None
            else SdbImporter(sdb_import_config, self.root)
        )
        self.export_sync_lock = RLock()
        self._operators = {}
        self._sample_cache = {}
        self._public_system_cache = OrderedDict()
        self._public_system_cache_lock = RLock()
        # Served documents contain plot data and are a few MB each even with
        # shared Bokeh resources. Keep a useful recent working set without
        # turning a long-running catalog process into an unbounded HTML cache.
        self._public_system_cache_size = 24

    def render_template(self, template_name, **values):
        values.setdefault('review_base_url', self.review_base_url)
        return _render(template_name, **values)

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

    def target_fit_status(self, identifier, *, since=None):
        rows = self.packages(identifiers=(str(identifier),))
        if len(rows) != 1:
            raise KeyError(
                'fitting package is not uniquely identified: {}'.format(
                    identifier
                )
            )
        row = rows[0]
        jobs = self.state.jobs(package=self.root / row['package_path'])
        if since:
            jobs = tuple(job for job in jobs if job.created_at >= str(since))
        return {
            'package': row['package_path'],
            'package_status': row['status'],
            'job': None if not jobs else asdict(jobs[0]),
            'system_url': _system_url(row['package_path'], admin=True),
        }

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
        dataset = FitDataset.from_path(directory)
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
                         per_page=100, principal=None):
        packages = self.packages(query=query)
        if status:
            packages = tuple(
                row for row in packages if row['status'] == status
            )
        total = len(packages)
        page_count = max(1, (total + per_page - 1) // per_page)
        page = min(max(1, int(page)), page_count)
        start = (page - 1) * per_page
        return self.render_template(
            'catalog_dashboard.html', title='SED catalogue',
            packages=packages[start:start + per_page], query=query or '',
            summary=self.summary(), samples=self.exports.samples(),
            status_filter=status or '', total=total, page=page,
            page_count=page_count, admin=True, principal=principal,
            search_action='/admin/search', page_base='/admin',
            system_base='/admin/system',
            import_enabled=self.importer is not None, token=self.token,
            auth_enabled=self.auth.proxy_enabled,
        )

    def render_activity(self, principal, *, days=30):
        """Render the admin-only inventory, use, and action overview."""

        days = int(days)
        if days not in {0, 7, 30, 90}:
            raise ValueError('activity period must be 7, 30, 90, or all')
        sdb_path = (
            None if self.auth.access is None else self.auth.access.path
        )
        report = activity_report(
            index_summary=asdict(self.index.summary()),
            packages=self.packages(), state_path=self.state.path,
            sdb_path=sdb_path, access_log=self.access_log, days=days,
        )
        return self.render_template(
            'catalog_activity.html', title='Site activity', report=report,
            principal=principal, admin=True, query='',
            search_action='/admin/search', request_path='/admin/activity',
            auth_enabled=self.auth.proxy_enabled,
        )

    def render_public_home(self, principal=None):
        return self.render_template(
            'catalog_public_home.html', title='sdb', query='',
            search_action='/search', principal=principal,
            auth_enabled=self.auth.proxy_enabled, request_path='/',
        )

    def render_public_search(self, query, *, page=1, per_page=100,
                             principal=None):
        query = str(query or '').strip()
        if not query:
            return self.render_public_home(principal)
        packages = self.packages(query=query)
        total = len(packages)
        page_count = max(1, (total + per_page - 1) // per_page)
        page = min(max(1, int(page)), page_count)
        start = (page - 1) * per_page
        return self.render_template(
            'catalog_public_search.html', title='Search results', query=query,
            packages=packages[start:start + per_page], total=total,
            page=page, page_count=page_count, search_action='/search',
            page_base='/search', system_base='/system',
            principal=principal, auth_enabled=self.auth.proxy_enabled,
            request_path='/search?q={}'.format(quote(query, safe='')),
            import_enabled=self.importer is not None, token=self.token,
        )

    def queue_acquisition(self, query, *, actor):
        if self.importer is None:
            raise ValueError('target import is not configured')
        job_id = self.state.queue_acquisition(query, actor=actor)
        return {'job_id': job_id, 'status': 'queued'}

    def acquisition(self, job_id, principal):
        job = self.state.acquisition(job_id)
        is_admin = principal.is_admin
        if job.actor != principal.username and not is_admin:
            raise AuthorizationDenied(
                'this target import belongs to another user'
            )
        return job

    def render_acquisition(self, job_id, principal):
        job = self.acquisition(job_id, principal)
        package_rows = []
        if job.result is not None:
            available = {
                row['package_path']: row for row in self.packages()
            }
            package_rows = [
                available[name] for name in job.result.get('packages', ())
                if name in available
            ]
        fit_active = any(
            row.get('job') is not None
            and row['job'].status in {'queued', 'running'}
            for row in package_rows
        )
        return self.render_template(
            'catalog_import_job.html',
            title='Import {}'.format(job.query), job=job,
            packages=package_rows,
            auto_refresh=(job.status in {'queued', 'running'} or fit_active),
            principal=principal, search_action='/search',
            system_base='/system', auth_enabled=self.auth.proxy_enabled,
            request_path='/imports/{}'.format(job.id),
            admin_access=principal.is_admin,
        )

    def queue_sample_acquisition(self, sample_name, targets, *, actor,
                                 note=None, visibility='private',
                                 operation='add', owner_username=None):
        if self.importer is None:
            raise ValueError('sample import is not configured')
        job_id = self.state.queue_sample_acquisition(
            sample_name, targets, actor=actor, note=note,
            visibility=visibility, operation=operation,
            owner_username=owner_username,
        )
        return {'job_id': job_id, 'status': 'queued'}

    def check_sample_targets(self, targets):
        if self.importer is None:
            raise ValueError('sample import is not configured')
        return self.importer.check_sample_targets(targets)

    def update_sample_collaborator(self, sample_name, username, *, operation,
                                   actor):
        if self.importer is None:
            raise ValueError('sample collaboration is not configured')
        return self.importer.update_sample_collaborator(
            sample_name, username, operation=operation, actor=actor,
        )

    def sample_acquisition(self, job_id, principal):
        job = self.state.sample_acquisition(job_id)
        if job.actor != principal.username and not principal.is_admin:
            raise AuthorizationDenied(
                'this sample import belongs to another user'
            )
        return job

    def render_new_sample(self, principal):
        return self.render_template(
            'catalog_sample_new.html', title='Create sample',
            principal=principal, search_action='/search', token=self.token,
            auth_enabled=self.auth.proxy_enabled, request_path='/samples/new',
            import_enabled=self.importer is not None,
        )

    def render_sample_acquisition(self, job_id, principal):
        job = self.sample_acquisition(job_id, principal)
        return self.render_template(
            'catalog_sample_import_job.html',
            title='Import {}'.format(job.sample_name), job=job,
            auto_refresh=job.status in {'queued', 'running'},
            principal=principal, search_action='/search',
            auth_enabled=self.auth.proxy_enabled,
            request_path='/sample-imports/{}'.format(job.id),
            admin_access=principal.is_admin,
        )

    def register_acquisition_export(self, result, *, actor):
        """Index complete exported packages and queue their automatic fits."""

        return self._register_acquisition_export(
            result, actor=actor,
            reason='Fit newly imported {}'.format(result['query']),
        )

    def register_sample_acquisition_export(self, result, *, actor):
        """Index one completed sample export and queue its needed fits."""

        registered = self._register_acquisition_export(
            result, actor=actor,
            reason='Fit newly imported sample {}'.format(result['sample']),
        )
        # ExportCatalog reads manifest snapshots directly, but discard any
        # projections created before this manifest appeared.
        self._sample_cache.clear()
        return registered

    def _register_acquisition_export(self, result, *, actor, reason):

        package_names = tuple(dict.fromkeys(
            _safe_relative(value) for value in result.get('packages', ())
        ))
        fit_job_ids = []
        with self.export_sync_lock:
            for name in package_names:
                package = (self.root / name).resolve()
                if self.root not in package.parents or not package.is_dir():
                    raise ValueError(
                        'sdb exported an invalid package path: {}'.format(name)
                    )
                self.index.update_package(package)
            active = {
                job.package_path for job in self.state.jobs()
                if job.status in {'queued', 'running'}
            }
            for name in package_names:
                row = self.package(name)
                directory = self.root / name
                if row['status'] not in {'unfit', 'stale'}:
                    continue
                if str(directory) in active:
                    continue
                queued = self.queue_package(
                    name, actor=actor,
                    reason=reason,
                    remember_actor=False,
                )
                fit_job_ids.append(queued['job_id'])
                active.add(str(directory))
        return dict(result, fit_job_ids=fit_job_ids)

    def render_samples(self, principal):
        samples = []
        for sample in self.exports.samples():
            if not self.auth.can_access_sample(principal, sample.name):
                continue
            policy = (
                None if not self.auth.proxy_enabled
                else self.auth.access.sample(
                    sample.name,
                    principal.username if principal.authenticated else None,
                )
            )
            if policy is None or policy.visibility == 'public':
                access = 'public'
            elif policy.owner_username == principal.username:
                access = 'yours'
            elif policy.granted:
                access = 'shared'
            else:
                access = 'private'
            samples.append({'sample': sample, 'access': access})
        return self.render_template(
            'catalog_samples.html', title='Samples',
            samples=tuple(samples), principal=principal,
            search_action='/search',
            admin_access=principal.is_admin,
            auth_enabled=self.auth.proxy_enabled, request_path='/samples',
            import_enabled=self.importer is not None,
        )

    def render_sample(self, name, principal):
        sample = self.exports.sample(name)
        policy = (
            None if not self.auth.proxy_enabled
            else self.auth.access.sample(
                sample.name,
                principal.username if principal.authenticated else None,
            )
        )
        can_share = (
            self.importer is not None
            and self.auth.can_manage_sample_collaborators(
                principal, sample.name,
            )
        )
        collaborators = (
            self.auth.access.sample_collaborators(sample.name)
            if can_share else ()
        )
        projection = self.sample_projection(sample)
        packages = self.packages(identifiers=sample.identifiers)
        projected = {system.package_path for system in projection.systems}
        labels = {}
        search_terms = {}
        for package in packages:
            path = self.root / package['package_path'] / 'identities.yml'
            if not path.is_file():
                continue
            for identifier, aliases in read_identity_manifest(path).items():
                search_terms[identifier] = tuple(dict.fromkeys((
                    identifier, *aliases,
                )))
        for system in projection.systems:
            for identifier, aliases in system.identity_aliases:
                search_terms[identifier] = tuple(dict.fromkeys((
                    identifier, *aliases,
                )))
            for target in system.targets:
                aliases = set(target.search_terms)
                for identifier in set(sample.identifiers) & aliases:
                    labels[identifier] = target.label
                    search_terms[identifier] = target.search_terms
        members = tuple({
            'identifier': identifier,
            'label': labels.get(identifier, identifier),
            'filter_text': ' '.join(search_terms.get(
                identifier, (labels.get(identifier, identifier), identifier),
            )),
        } for identifier in sample.identifiers)
        return self.render_template(
            'catalog_sample.html', title=sample.name, sample=sample,
            projection=projection, rows=sample_table_rows(projection),
            pending=tuple(row for row in packages
                          if row['package_path'] not in projected),
            sample_members=members,
            can_edit=self.auth.can_edit_sample(principal, sample.name),
            can_share=can_share, collaborators=collaborators,
            sample_owner=(
                policy.owner_username if policy is not None
                else sample.owner_username
            ),
            token=self.token,
            principal=principal, search_action='/search',
            system_base='/system',
            admin_access=principal.is_admin,
            auth_enabled=self.auth.proxy_enabled,
            request_path='/samples/{}/'.format(quote(sample.name, safe='')),
        )

    def render_sample_plots(self, name, principal):
        sample = self.exports.sample(name)
        projection = self.sample_projection(sample)
        script, div = sample_plot_components(projection)
        return self.render_template(
            'catalog_sample_plot.html', title='{} plots'.format(sample.name),
            sample=sample, projection=projection,
            plot_script=script, plot_div=div,
            bokeh_js=CDN.render_js(), bokeh_css=CDN.render_css(),
            principal=principal, search_action='/search',
            admin_access=principal.is_admin,
            auth_enabled=self.auth.proxy_enabled,
            request_path='/samples/{}/plots'.format(quote(sample.name, safe='')),
        )

    def render_sample_flux_size(self, name, principal):
        sample = self.exports.sample(name)
        projection = self.sample_projection(sample)
        script, div = flux_size_plot_components(projection)
        return self.render_template(
            'catalog_sample_flux_size.html',
            title='{} flux and size'.format(sample.name), sample=sample,
            projection=projection, plot_script=script, plot_div=div,
            bokeh_js=CDN.render_js(), bokeh_css=CDN.render_css(),
            principal=principal, search_action='/search',
            admin_access=principal.is_admin,
            auth_enabled=self.auth.proxy_enabled,
            request_path='/samples/{}/flux-size'.format(
                quote(sample.name, safe='')
            ),
        )

    def sample_projection(self, sample):
        stamp = (
            sample.run_id,
            sample.manifest_path.stat().st_mtime_ns,
            self.index.path.stat().st_mtime_ns,
        )
        cached = self._sample_cache.get(sample.name)
        if cached is not None and cached[0] == stamp:
            return cached[1]
        projection = project_sample(
            self.index, sample.name, sample.identifiers,
            sed_base_url='/system',
        )
        self._sample_cache[sample.name] = (stamp, projection)
        return projection

    def sample_votable_bytes(self, name):
        sample = self.exports.sample(name)
        stream = io.BytesIO()
        sample_votable(self.sample_projection(sample)).write(
            stream, format='votable'
        )
        return stream.getvalue()

    def sample_photometry_csv(self, name):
        import csv
        sample = self.exports.sample(name)
        stream = io.StringIO(newline='')
        writer = csv.DictWriter(stream, fieldnames=PHOTOMETRY_FIELDS)
        writer.writeheader()
        writer.writerows(photometry_rows(self.sample_projection(sample)))
        return stream.getvalue().encode('utf-8')

    def render_jobs(self, principal=None):
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
        return self.render_template(
            'catalog_jobs.html', title='Jobs', jobs=jobs,
            acquisitions=self.state.acquisitions(), admin=True,
            principal=principal, search_action='/admin/search',
            system_base='/admin/system',
        )

    def render_package_state(self, package_path, *, actor=None, admin=False,
                             principal=None):
        row = self.package(package_path)
        return self.render_template(
            'catalog_package.html' if admin else 'catalog_public_package.html',
            title=' / '.join(row['target_labels']) or row['package_path'],
            package=row, token=self.token, actor=actor or self.actor,
            actor_locked=self.auth.proxy_enabled, admin=admin,
            principal=principal,
            search_action='/admin/search' if admin else '/search',
        )

    def render_public_system(self, package_path, *, principal=None):
        principal = principal or self.auth.principal({})
        row = self.package(package_path)
        if row['status'] in {'unfit', 'needs selection'}:
            return self.render_package_state(
                package_path, principal=principal,
            )
        directory = self.package_directory(package_path)
        cache_key = (
            row['package_path'], principal.username, bool(principal.is_admin),
        )
        stamp = _public_system_stamp(directory)
        with self._public_system_cache_lock:
            cached = self._public_system_cache.get(cache_key)
            if cached is not None and cached[0] == stamp:
                self._public_system_cache.move_to_end(cache_key)
                return cached[1]
        with tempfile.NamedTemporaryFile(
            prefix='sdf-public-', suffix='.html', delete=False,
        ) as temporary:
            name = temporary.name
        try:
            output = FitManager(
                directory
            ).render(
                name,
                navigation=self.sed_navigation(
                    package_path, principal=principal, admin=False,
                ),
                inline=False,
            )
            document = output.read_text(encoding='utf-8')
            with self._public_system_cache_lock:
                self._public_system_cache[cache_key] = (stamp, document)
                self._public_system_cache.move_to_end(cache_key)
                while (
                    len(self._public_system_cache)
                    > self._public_system_cache_size
                ):
                    self._public_system_cache.popitem(last=False)
            return document
        finally:
            try:
                Path(name).unlink()
            except FileNotFoundError:
                pass

    def sed_navigation(self, package_path, *, principal, admin):
        public_url = _system_url(package_path)
        admin_url = _system_url(package_path, admin=True)
        return {
            'principal': principal,
            'admin_access': principal.is_admin,
            'auth_enabled': self.auth.proxy_enabled,
            'request_path': admin_url if admin else public_url,
            'search_action': '/admin/search' if admin else '/search',
            'admin_url': None if admin else admin_url,
            'review_base_url': self.review_base_url,
            'review_url': (
                None if self.review_base_url is None else
                self.review_base_url + '/target/' + quote(
                    PurePosixPath(package_path).name, safe=''
                )
            ),
            'external_route_root': '',
        }


def make_catalog_server(root, *, index_path=None, state_path=None, actor=None,
                        token=None, rebuild=False, host='127.0.0.1', port=0,
                        auth_user_header=None, auth_database=None,
                        auth_session_secret=None,
                        auth_cookie_secure=False, sdb_import_config=None,
                        review_base_url=None, access_log=None):
    if host not in {'127.0.0.1', 'localhost', '::1'}:
        raise ValueError('catalog service must bind to a loopback address')
    service = CatalogService(
        root, index_path=index_path, state_path=state_path, actor=actor,
        token=token, rebuild=rebuild, auth_user_header=auth_user_header,
        auth_database=auth_database,
        auth_session_secret=auth_session_secret,
        auth_cookie_secure=auth_cookie_secure,
        sdb_import_config=sdb_import_config,
        review_base_url=review_base_url,
        access_log=access_log,
    )
    handler = partial(_CatalogHandler, service=service)
    server = ThreadingHTTPServer((host, port), handler)
    server.catalog_service = service
    return server


def serve_catalog(root, *, index_path=None, state_path=None, actor=None,
                  rebuild=False, host='127.0.0.1', port=8765,
                  worker=True, workers=1, watch=True, autofit_existing=False,
                  poll_seconds=5.0, auth_user_header=None,
                  auth_database=None, auth_session_secret=None,
                  auth_cookie_secure=False,
                  sdb_import_config=None, review_base_url=None,
                  access_log=None):
    server = make_catalog_server(
        root, index_path=index_path, state_path=state_path, actor=actor,
        rebuild=rebuild, host=host, port=port,
        auth_user_header=auth_user_header,
        auth_database=auth_database,
        auth_session_secret=auth_session_secret,
        auth_cookie_secure=auth_cookie_secure,
        sdb_import_config=sdb_import_config,
        review_base_url=review_base_url,
        access_log=access_log,
    )
    address, actual_port = server.server_address[:2]
    print('sdf-fits service: http://{}:{}/'.format(address, actual_port))
    print('Result index: {}'.format(server.catalog_service.index.path))
    print('Administrative state: {}'.format(server.catalog_service.state.path))
    runtime = CatalogRuntime(
        server.catalog_service, worker=worker, workers=workers, watch=watch,
        autofit_existing=autofit_existing, poll_seconds=poll_seconds,
    )
    previous_sigterm = signal.getsignal(signal.SIGTERM)
    signal.signal(signal.SIGTERM, _raise_keyboard_interrupt)
    try:
        runtime.start()
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        signal.signal(signal.SIGTERM, previous_sigterm)
        server.server_close()
        runtime.stop()


def _raise_keyboard_interrupt(signum, frame):
    """Route service termination through the normal worker cleanup path."""

    raise KeyboardInterrupt


class CatalogRuntime:
    """Supervise the optional watcher thread and isolated fit workers."""

    def __init__(self, service, *, worker=True, workers=1, watch=True,
                 autofit_existing=False, poll_seconds=5.0):
        self.service = service
        self.worker_enabled = worker
        self.worker_count = int(workers) if worker else 0
        if self.worker_count < 0 or (worker and self.worker_count == 0):
            raise ValueError('workers must be a positive integer')
        self.watch_enabled = watch
        self.autofit_existing = autofit_existing
        self.poll_seconds = poll_seconds
        self.stop_event = multiprocessing.Event()
        self.watch_stop = Event()
        self.worker_processes = {}
        self.worker_ids = {}
        self.worker_lock = None
        self.watch_thread = None
        self.supervisor_thread = None
        self.import_thread = None

    def start(self):
        if self.worker_enabled:
            self.worker_lock = InterProcessLock(
                Path(str(self.service.state.path) + '.workers.lock'),
                description='fit worker service',
            )
            self.worker_lock.acquire()
            try:
                # The process-wide lock makes global recovery safe. Later
                # pool-slot restarts recover only work owned by that slot.
                recovered = self.service.state.recover_running_jobs()
                if recovered:
                    print('Recovered {} interrupted fit job(s)'.format(
                        recovered
                    ))
                for slot in range(self.worker_count):
                    self._start_worker(slot)
                self.supervisor_thread = Thread(
                    target=self._supervise_worker,
                    name='sdf-worker-supervisor', daemon=True,
                )
                self.supervisor_thread.start()
            except Exception:
                self.worker_lock.release()
                self.worker_lock = None
                raise
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
        if getattr(self.service, 'importer', None) is not None:
            recovered = self.service.state.recover_running_acquisitions()
            recovered_samples = (
                self.service.state.recover_running_sample_acquisitions()
            )
            if recovered:
                print('Recovered {} interrupted target import(s)'.format(
                    recovered
                ))
            if recovered_samples:
                print('Recovered {} interrupted sample import(s)'.format(
                    recovered_samples
                ))
            self.import_thread = Thread(
                target=run_acquisition_worker,
                args=(self.service.state, self.service.importer),
                kwargs={
                    'on_export': self.service.register_acquisition_export,
                    'on_sample_export': (
                        self.service.register_sample_acquisition_export
                    ),
                    'poll_seconds': self.poll_seconds,
                    'stop_event': self.watch_stop,
                    'worker_id': 'sdf-sdb-import',
                },
                name='sdf-sdb-import', daemon=True,
            )
            self.import_thread.start()
            print('sdb target import: enabled')

    def stop(self):
        self.watch_stop.set()
        self.stop_event.set()
        if self.watch_thread is not None:
            self.watch_thread.join(timeout=self.poll_seconds + 1)
        if self.import_thread is not None:
            self.import_thread.join(timeout=self.poll_seconds + 1)
        if self.supervisor_thread is not None:
            self.supervisor_thread.join(timeout=self.poll_seconds + 1)
        for process in self.worker_processes.values():
            process.join(timeout=self.poll_seconds + 1)
        for process in self.worker_processes.values():
            if process.is_alive():
                process.terminate()
                process.join(timeout=2)
        if self.worker_lock is not None:
            self.worker_lock.release()
            self.worker_lock = None

    def _start_worker(self, slot):
        worker_id = '{}-{}'.format(id(self), slot + 1)
        process = multiprocessing.Process(
            target=_worker_main,
            args=(str(self.service.state.path), self.stop_event,
                  self.poll_seconds, worker_id),
            name='sdf-fit-worker-{}'.format(slot + 1), daemon=False,
        )
        process.start()
        self.worker_processes[slot] = process
        self.worker_ids[slot] = worker_id
        print('Fit worker {}: pid {}'.format(slot + 1, process.pid))

    def _supervise_worker(self):
        while not self.watch_stop.wait(self.poll_seconds):
            if self.stop_event.is_set():
                return
            for slot, process in tuple(self.worker_processes.items()):
                if process.is_alive():
                    continue
                worker_id = self.worker_ids[slot]
                recovered = self.service.state.recover_running_jobs(
                    worker_id=worker_id
                )
                print('Fit worker {} stopped; recovered {} job(s); restarting'.
                      format(slot + 1, recovered))
                self._start_worker(slot)


def _worker_main(state_path, stop_event, poll_seconds, worker_id):
    # The parent service owns shutdown and signals workers through stop_event.
    # Ignoring terminal interrupts here avoids one traceback per pool process.
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    run_worker(
        OperatorState(state_path).initialize(),
        poll_seconds=poll_seconds, stop_event=stop_event,
        recover=False, worker_id=worker_id,
    )


class _CatalogHandler(BaseHTTPRequestHandler):
    def __init__(self, *args, service, **kwargs):
        self.service = service
        super().__init__(*args, **kwargs)

    def do_GET(self):
        parsed = urlparse(self.path)
        try:
            principal = self.service.auth.principal(self.headers)
            if parsed.path == '/external/simbad/id':
                values = parse_qs(parsed.query)
                name = values.get('name', [''])[0].strip()
                if not name:
                    raise ValueError('SIMBAD name is required')
                self._redirect_external(
                    self.service.external_services.simbad_identifier_url(name)
                )
                return
            if parsed.path == '/external/vizier/catalog':
                values = parse_qs(parsed.query)
                source = values.get('source', [''])[0].strip()
                if not source:
                    raise ValueError('VizieR source is required')
                self._redirect_external(
                    self.service.external_services.vizier_catalog_url(source)
                )
                return
            if parsed.path == '/external/vizier/row':
                values = parse_qs(parsed.query)
                source = values.get('source', [''])[0].strip()
                column = values.get('column', [''])[0].strip()
                value = values.get('value', [''])[0].strip()
                if not source or not column or not value:
                    raise ValueError(
                        'VizieR source, column, and value are required'
                    )
                self._redirect_external(
                    self.service.external_services.vizier_row_url(
                        source, column, value,
                    )
                )
                return
            if parsed.path == '/':
                self._html(self.service.render_public_home(principal))
                return
            if parsed.path == '/search':
                values = parse_qs(parsed.query)
                query = values.get('q', [''])[0].strip()
                page = int(values.get('page', ['1'])[0])
                self._html(self.service.render_public_search(
                    query, page=page, principal=principal
                ))
                return
            if parsed.path.startswith('/imports/'):
                self.service.auth.require_authenticated(principal)
                job_id = parsed.path.removeprefix('/imports/').strip('/')
                if not job_id.isdigit():
                    raise KeyError('target import not found')
                self._html(self.service.render_acquisition(
                    int(job_id), principal
                ))
                return
            if parsed.path.startswith('/sample-imports/'):
                self.service.auth.require_authenticated(principal)
                job_id = parsed.path.removeprefix(
                    '/sample-imports/'
                ).strip('/')
                if not job_id.isdigit():
                    raise KeyError('sample import not found')
                self._html(self.service.render_sample_acquisition(
                    int(job_id), principal
                ))
                return
            if parsed.path == '/auth/login':
                self.service.auth.require_authenticated(principal)
                values = parse_qs(parsed.query)
                destination = _safe_next(values.get('next', ['/samples'])[0])
                token = self.service.auth.issue_session(principal.username)
                cookie = (
                    'sdf_session={}; Path=/; HttpOnly; SameSite=Lax; '
                    'Max-Age={}'.format(
                        token, self.service.auth.session_seconds
                    )
                )
                if self.service.auth_cookie_secure:
                    cookie += '; Secure'
                self._redirect(destination, headers=(('Set-Cookie', cookie),))
                return
            if parsed.path == '/auth/logout':
                self._redirect('/samples', headers=((
                    'Set-Cookie',
                    'sdf_session=; Path=/; HttpOnly; SameSite=Lax; Max-Age=0',
                ),))
                return
            if parsed.path in {'/samples', '/samples/'}:
                self._html(self.service.render_samples(principal))
                return
            if parsed.path == '/samples/new':
                self.service.auth.require_authenticated(principal)
                self._html(self.service.render_new_sample(principal))
                return
            if parsed.path.startswith('/samples/'):
                name, asset = _sample_parts(parsed.path)
                self.service.auth.require_sample(principal, name)
                if not asset:
                    self._html(self.service.render_sample(name, principal))
                    return
                if asset == 'plots':
                    self._html(self.service.render_sample_plots(
                        name, principal
                    ))
                    return
                if asset == 'flux-size':
                    self._html(self.service.render_sample_flux_size(
                        name, principal
                    ))
                    return
                if asset == 'sample.xml':
                    self._send(
                        HTTPStatus.OK, self.service.sample_votable_bytes(name),
                        'application/x-votable+xml; charset=utf-8',
                        headers=((
                            'Content-Disposition',
                            'attachment; filename="sample.xml"',
                        ),),
                    )
                    return
                if asset == 'photometry.csv':
                    self._send(
                        HTTPStatus.OK,
                        self.service.sample_photometry_csv(name),
                        'text/csv; charset=utf-8',
                        headers=((
                            'Content-Disposition',
                            'attachment; filename="photometry.csv"',
                        ),),
                    )
                    return
                raise KeyError('unknown sample page')
            if parsed.path in {'/admin', '/admin/'}:
                self.service.auth.require_admin(principal)
                values = parse_qs(parsed.query)
                query = values.get('q', [''])[0].strip()
                status = values.get('status', [''])[0].strip()
                page = int(values.get('page', ['1'])[0])
                self._html(self.service.render_dashboard(
                    query or None, status=status or None, page=page,
                    principal=principal,
                ))
                return
            if parsed.path == '/admin/search':
                self.service.auth.require_admin(principal)
                values = parse_qs(parsed.query)
                query = values.get('q', [''])[0].strip()
                status = values.get('status', [''])[0].strip()
                page = int(values.get('page', ['1'])[0])
                self._html(self.service.render_dashboard(
                    query or None, status=status or None, page=page,
                    principal=principal,
                ))
                return
            if parsed.path == '/admin/jobs':
                self.service.auth.require_admin(principal)
                self._html(self.service.render_jobs(principal))
                return
            if parsed.path == '/admin/activity':
                self.service.auth.require_admin(principal)
                values = parse_qs(parsed.query)
                days = int(values.get('days', ['30'])[0])
                self._html(self.service.render_activity(
                    principal, days=days,
                ))
                return
            if (parsed.path.startswith('/admin/target/')
                    and parsed.path.endswith('/api/status')):
                self.service.auth.require_admin(principal)
                identifier = unquote(
                    parsed.path.removeprefix('/admin/target/').removesuffix(
                        '/api/status'
                    ).rstrip('/')
                )
                values = parse_qs(parsed.query)
                since = values.get('since', [None])[0]
                self._json(self.service.target_fit_status(
                    identifier, since=since,
                ))
                return
            if parsed.path.startswith('/admin/target/'):
                self.service.auth.require_admin(principal)
                identifier = unquote(
                    parsed.path.removeprefix('/admin/target/')
                )
                rows = self.service.packages(identifiers=(identifier,))
                if len(rows) == 1:
                    self.send_response(HTTPStatus.SEE_OTHER)
                    self.send_header(
                        'Location', _system_url(
                            rows[0]['package_path'], admin=True
                        )
                    )
                    self.end_headers()
                    return
                raise KeyError(
                    'fitting package is not uniquely identified: {}'.format(
                        identifier
                    )
                )
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
                self._html(self.service.render_public_search(
                    identifier, principal=principal
                ),
                           status=HTTPStatus.NOT_FOUND if not rows else HTTPStatus.OK)
                return
            if parsed.path.startswith('/system/'):
                package, asset = _system_parts(parsed.path, '/system/')
                if not asset:
                    self._html(self.service.render_public_system(
                        package, principal=principal,
                    ))
                    return
                self._static(
                    self.service.package_directory(package), asset,
                    public=True,
                )
                return
            if parsed.path.startswith('/admin/system/'):
                self.service.auth.require_admin(principal)
                package, asset = _system_parts(
                    parsed.path, '/admin/system/'
                )
                row = self.service.package(package)
                if not asset:
                    if row['status'] in {'unfit', 'needs selection'}:
                        self._html(self.service.render_package_state(
                            package, actor=principal.username, admin=True,
                            principal=principal,
                        ))
                    else:
                        prefix = _system_url(
                            package, admin=True
                        ).rstrip('/') + '/api'
                        self._html(self.service.operator(package).render_page(
                            operator_api=prefix, actor=principal.username,
                            actor_locked=self.service.auth.proxy_enabled,
                            navigation=self.service.sed_navigation(
                                package, principal=principal, admin=True,
                            ),
                        ))
                    return
                self._static(self.service.package_directory(package), asset)
                return
            if parsed.path == '/favicon.ico':
                self._send(HTTPStatus.NO_CONTENT, b'', 'image/x-icon')
                return
            self._error(KeyError('page not found'), HTTPStatus.NOT_FOUND)
        except AuthenticationRequired as error:
            if self.service.auth.proxy_enabled:
                self._redirect('/auth/login?next={}'.format(
                    quote(self.path, safe='')
                ))
            else:
                self._error(error, HTTPStatus.UNAUTHORIZED)
        except AuthorizationDenied as error:
            self._error(error, HTTPStatus.FORBIDDEN)
        except (KeyError, ValueError, OSError) as error:
            self._error(error, HTTPStatus.NOT_FOUND)

    def do_POST(self):
        parsed = urlparse(self.path)
        try:
            principal = self.service.auth.principal(self.headers)
            if parsed.path == '/imports':
                self.service.auth.require_authenticated(principal)
                if self.headers.get('X-SDF-Token') != self.service.token:
                    raise PermissionError('invalid request token')
                payload = self._payload()
                result = self.service.queue_acquisition(
                    payload.get('query'), actor=principal.username,
                )
                self._json(result, status=HTTPStatus.ACCEPTED)
                return
            if parsed.path == '/sample-imports':
                self.service.auth.require_authenticated(principal)
                if self.headers.get('X-SDF-Token') != self.service.token:
                    raise PermissionError('invalid request token')
                payload = self._payload()
                visibility = payload.get('visibility', 'private')
                sample_name = payload.get('sample_name')
                owner_username = principal.username
                if self.service.auth.proxy_enabled:
                    policy = self.service.auth.access.sample(
                        sample_name, principal.username,
                    )
                    if policy is not None:
                        self.service.auth.require_sample_edit(
                            principal, sample_name,
                        )
                        owner_username = policy.owner_username
                else:
                    try:
                        owner_username = self.service.exports.sample(
                            sample_name
                        ).owner_username
                    except KeyError:
                        owner_username = principal.username
                if visibility == 'public' and not principal.is_admin:
                    raise AuthorizationDenied(
                        'only administrators can publish a sample'
                    )
                result = self.service.queue_sample_acquisition(
                    sample_name, payload.get('targets', ()),
                    actor=principal.username, note=payload.get('note'),
                    visibility=visibility, owner_username=owner_username,
                )
                self._json(result, status=HTTPStatus.ACCEPTED)
                return
            if parsed.path == '/sample-edits':
                self.service.auth.require_authenticated(principal)
                if self.headers.get('X-SDF-Token') != self.service.token:
                    raise PermissionError('invalid request token')
                payload = self._payload()
                sample_name = payload.get('sample_name')
                self.service.auth.require_sample_edit(principal, sample_name)
                policy = (
                    None if not self.service.auth.proxy_enabled
                    else self.service.auth.access.sample(
                        sample_name, principal.username,
                    )
                )
                exported_owner = None
                if policy is None:
                    exported_owner = self.service.exports.sample(
                        sample_name
                    ).owner_username
                operation = str(payload.get('operation') or '').strip()
                if operation not in {'add', 'remove'}:
                    raise ValueError('sample operation must be add or remove')
                result = self.service.queue_sample_acquisition(
                    sample_name, payload.get('targets', ()),
                    actor=principal.username, operation=operation,
                    owner_username=(
                        policy.owner_username if policy is not None
                        else exported_owner
                    ),
                )
                self._json(result, status=HTTPStatus.ACCEPTED)
                return
            if parsed.path == '/sample-collaborators':
                self.service.auth.require_authenticated(principal)
                if self.headers.get('X-SDF-Token') != self.service.token:
                    raise PermissionError('invalid request token')
                payload = self._payload()
                sample_name = payload.get('sample_name')
                if not self.service.auth.can_manage_sample_collaborators(
                        principal, sample_name):
                    raise AuthorizationDenied(
                        'managing collaborators for this sample is not permitted'
                    )
                operation = str(payload.get('operation') or '').strip()
                if operation not in {'share', 'unshare'}:
                    raise ValueError(
                        'sample collaborator operation must be share or unshare'
                    )
                result = self.service.update_sample_collaborator(
                    sample_name, payload.get('username'),
                    operation=operation, actor=principal.username,
                )
                self._json(result)
                return
            if parsed.path == '/sample-checks':
                self.service.auth.require_authenticated(principal)
                if self.headers.get('X-SDF-Token') != self.service.token:
                    raise PermissionError('invalid request token')
                payload = self._payload()
                self._json(self.service.check_sample_targets(
                    payload.get('targets', ())
                ))
                return
            self.service.auth.require_admin(principal)
            package, asset = _system_parts(
                parsed.path, '/admin/system/'
            )
            if not asset.startswith('api/'):
                raise KeyError('unknown catalog endpoint')
            if self.headers.get('X-SDF-Token') != self.service.token:
                raise PermissionError('invalid operator token')
            payload = self._payload()
            payload['actor'] = principal.username
            action = asset.removeprefix('api/')
            if action == 'queue':
                result = self.service.queue_package(
                    package, actor=principal.username,
                    reason=payload.get('reason'),
                    remember_actor=not self.service.auth.proxy_enabled,
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
        except AuthenticationRequired as error:
            self._error(error, HTTPStatus.UNAUTHORIZED, json_response=True)
        except (AuthorizationDenied, PermissionError) as error:
            self._error(error, HTTPStatus.FORBIDDEN, json_response=True)
        except (KeyError, ValueError, RuntimeError, OSError,
                json.JSONDecodeError) as error:
            self._error(error, HTTPStatus.BAD_REQUEST, json_response=True)

    def _payload(self):
        length = int(self.headers.get('Content-Length', '0'))
        if length <= 0 or length > 1024 * 1024:
            raise ValueError('invalid request size')
        return json.loads(self.rfile.read(length))

    def _static(self, directory, asset, *, public=False):
        relative = _safe_relative(asset)
        if public and relative.casefold() == 'operator.html':
            raise KeyError('package asset not found')
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

    def _redirect(self, location, *, headers=()):
        self.send_response(HTTPStatus.SEE_OTHER)
        self.send_header('Location', location)
        self.send_header('Content-Length', '0')
        self.send_header('Cache-Control', 'no-store')
        for name, value in headers:
            self.send_header(name, value)
        self.end_headers()

    def _redirect_external(self, location):
        self.send_response(HTTPStatus.FOUND)
        self.send_header('Location', location)
        self.send_header('Content-Length', '0')
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()

    def _error(self, error, status, json_response=False):
        if json_response:
            self._json({'error': str(error)}, status=status)
            return
        self._html(_render(
            'catalog_message.html', title='sdb',
            heading=HTTPStatus(status).phrase, message=str(error),
        ), status=status)

    def _send(self, status, body, content_type, *, headers=()):
        self.send_response(status)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        for name, value in headers:
            self.send_header(name, value)
        self.end_headers()
        self.wfile.write(body)


def _render(template_name, **values):
    environment = jinja2.Environment(
        autoescape=True,
        loader=jinja2.PackageLoader('sdf', package_path='www/templates'),
    )
    environment.filters['summary'] = summary_text
    return environment.get_template(template_name).render(**values)


def _public_system_stamp(directory):
    """Return the cheap package revision that affects a rendered SED page.

    Sampler workspaces and files produced by page rendering itself are omitted.
    Inputs, retained results, selection, and diagnostic images remain part of
    the revision, so exports and completed fits invalidate only their package.
    """

    directory = Path(directory)
    ignored_names = {
        'index.html', 'operator.html', 'plot-data.npz', '.sdf-fit.lock',
    }
    states = []
    for root, directories, filenames in os.walk(directory):
        directories[:] = [
            name for name in directories
            if name not in {'.work', 'limits-data'}
        ]
        root = Path(root)
        for filename in filenames:
            if filename in ignored_names:
                continue
            path = root / filename
            try:
                stat = path.stat()
            except FileNotFoundError:
                # Atomic fit/export updates can replace a file between the
                # directory listing and stat. The next request sees it.
                continue
            relative = path.relative_to(directory)
            states.append((
                relative.as_posix(), stat.st_mtime_ns, stat.st_size,
            ))
    return tuple(sorted(states))


def _safe_relative(value):
    path = PurePosixPath(unquote(str(value)))
    if path.is_absolute() or '..' in path.parts or not path.parts:
        raise ValueError('invalid catalog path')
    return path.as_posix()


def _system_url(package, admin=False):
    prefix = '/admin/system' if admin else '/system'
    return '{}/{}/'.format(prefix, quote(package, safe='/'))


def _sample_parts(path):
    prefix = '/samples/'
    if not path.startswith(prefix):
        raise KeyError('unknown sample endpoint')
    relative = path.removeprefix(prefix).strip('/')
    if not relative:
        raise KeyError('sample name is required')
    parts = relative.split('/')
    name = unquote(parts[0]).strip()
    if not name or name in {'.', '..'}:
        raise ValueError('invalid sample name')
    return name, '/'.join(parts[1:]).strip('/')


def _safe_next(value):
    value = str(value or '').strip()
    parsed = urlparse(value)
    if (not (value == '/samples' or value.startswith('/samples/'))
            or value.startswith('//')
            or parsed.scheme or parsed.netloc):
        return '/samples'
    return value


def _system_parts(path, prefix='/system/'):
    if not path.startswith(prefix):
        raise KeyError('unknown catalog endpoint')
    relative = path.removeprefix(prefix).lstrip('/')
    if not relative:
        raise KeyError('system package is required')
    parts = relative.split('/')
    # Export packages are currently one directory below the root. Preserve a
    # clear boundary here; nested package roots can later use an encoded key.
    package = _safe_relative(parts[0])
    asset = '/'.join(parts[1:]).strip('/')
    return package, asset
