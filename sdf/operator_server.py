"""Loopback-only HTTP boundary for the operator extension of SED pages."""

from __future__ import annotations

from functools import partial
import getpass
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import secrets
from urllib.parse import parse_qs, urlparse

from .dataset import FitDataset
from .fit_management import FitManager
from .model_policy import ModelFamilyRegistry, TargetModelPolicy
from .operator_state import default_operator_state


class OperatorService:
    """Validated operations used by the local HTML/JSON boundary."""

    def __init__(self, package, state, *, actor=None, token=None,
                 dataset=None):
        self.package = Path(package).expanduser().resolve()
        self.state = state.initialize()
        self.actor = actor or self.state.last_actor() or getpass.getuser()
        self.token = token or secrets.token_urlsafe(24)
        self._dataset_override = dataset
        self._seed_retained_policy()

    def dataset(self):
        return self._dataset_override or FitDataset.from_package(self.package)

    def _seed_retained_policy(self):
        manager = FitManager(self.package)
        candidates = manager.library.candidates()
        if not candidates:
            return
        try:
            default_id = manager.library.read_selection().default
            candidate = manager.library.find(default_id)
        except FileNotFoundError:
            candidate = max(
                candidates,
                key=lambda value: (
                    value.posterior.evidence is not None,
                    value.posterior.evidence
                    if value.posterior.evidence is not None else float('-inf'),
                ),
            )
        if candidate is None:
            return
        targets = {
            target.key: target for target in self.dataset().targets
        }
        for values in candidate.metadata.get('policy', {}).get('targets', ()):
            target_key = values.get('target_key')
            if target_key not in targets:
                continue
            policy = TargetModelPolicy(
                target_key=target_key,
                stellar_models=tuple(values['stellar_models']),
                disk_models=tuple(values.get('disk_models', ())),
                max_disks=values.get('max_disks', 0),
                search_mode='evidence-guided',
                reddening=values.get('reddening', False),
                origin='operator',
            )
            self.state.seed_policy(
                policy, package=self.package,
                main_id=targets[target_key].main_id,
            )

    def render_page(self, fits=None):
        output = FitManager(self.package).render(
            self.package / 'operator.html', fits=fits, operator=True,
            dataset=self.dataset(),
            policy_store=self.state, operator_api='/api',
            operator_token=self.token, operator_actor=self.actor,
        )
        return output.read_text(encoding='utf-8')

    def save_policy(self, payload):
        dataset = self.dataset()
        targets = {target.key: target for target in dataset.targets}
        target_key = str(payload.get('target_key', ''))
        if target_key not in targets:
            raise ValueError('unknown physical contributor {}'.format(target_key))
        actor = str(payload.get('actor') or self.actor).strip()
        if not actor:
            raise ValueError('actor is required')
        scope = _scope(payload.get('scope'))
        stellar_model = str(payload.get('stellar_model', '')).strip()
        disk_model = str(payload.get('disk_model', '')).strip()
        try:
            max_disks = int(payload.get('max_disks'))
        except (TypeError, ValueError) as error:
            raise ValueError('maximum disks must be an integer') from error
        search_mode = str(payload.get('search_mode', 'evidence-guided'))
        registry = ModelFamilyRegistry.from_config()
        registry.require(stellar_model, 'stellar')
        disk_models = () if not disk_model else (disk_model,)
        for model_name in disk_models:
            registry.require(model_name, 'disk')
        policy = TargetModelPolicy(
            target_key=target_key, stellar_models=(stellar_model,),
            disk_models=disk_models, max_disks=max_disks,
            search_mode=search_mode, origin='operator',
        )
        target = targets[target_key]
        reason = str(payload.get('reason', '')).strip() or _policy_reason(
            target.main_id or target.sdbid or target.key, policy,
        )
        job_id = self.state.save_policy_and_queue(
            policy, package=self.package, main_id=target.main_id, scope=scope,
            actor=actor, reason=reason,
        )
        self.actor = actor
        return {'job_id': job_id, 'status': 'queued'}

    def save_system(self, payload):
        """Save changed contributor policies and queue one system fit."""

        dataset = self.dataset()
        targets = {target.key: target for target in dataset.targets}
        actor = str(payload.get('actor') or self.actor).strip()
        if not actor:
            raise ValueError('actor is required')
        changes = payload.get('changes', ())
        if not isinstance(changes, list):
            raise ValueError('system policy changes must be a list')
        policies = []
        resets = []
        seen = set()
        registry = ModelFamilyRegistry.from_config()
        for change in changes:
            if not isinstance(change, dict):
                raise ValueError('each system policy change must be an object')
            target_key = str(change.get('target_key', ''))
            if target_key not in targets:
                raise ValueError(
                    'unknown physical contributor {}'.format(target_key)
                )
            if target_key in seen:
                raise ValueError(
                    'duplicate physical contributor {}'.format(target_key)
                )
            seen.add(target_key)
            scope = _scope(change.get('scope'))
            if change.get('reset'):
                resets.append((target_key, scope))
                continue
            stellar_model = str(change.get('stellar_model', '')).strip()
            disk_model = str(change.get('disk_model', '')).strip()
            try:
                max_disks = int(change.get('max_disks'))
            except (TypeError, ValueError) as error:
                raise ValueError('maximum disks must be an integer') from error
            search_mode = str(change.get('search_mode', 'evidence-guided'))
            registry.require(stellar_model, 'stellar')
            disk_models = () if not disk_model else (disk_model,)
            for model_name in disk_models:
                registry.require(model_name, 'disk')
            policies.append((TargetModelPolicy(
                target_key=target_key, stellar_models=(stellar_model,),
                disk_models=disk_models, max_disks=max_disks,
                search_mode=search_mode, origin='operator',
            ), targets[target_key].main_id, scope))
        reason = str(payload.get('reason', '')).strip()
        if not reason:
            labels = [
                target.main_id or target.sdbid or target.key
                for target in dataset.targets
            ]
            reason = 'Update model choices and fit {}'.format(
                ' + '.join(labels)
            )
        job_id = self.state.save_system_and_queue(
            policies, resets, package=self.package, actor=actor,
            reason=reason, target_keys=targets,
        )
        self.actor = actor
        return {'job_id': job_id, 'status': 'queued'}

    def reset_policy(self, payload):
        dataset = self.dataset()
        targets = {target.key: target for target in dataset.targets}
        target_key = str(payload.get('target_key', ''))
        if target_key not in targets:
            raise ValueError('unknown physical contributor {}'.format(target_key))
        actor = str(payload.get('actor') or self.actor).strip()
        if not actor:
            raise ValueError('actor is required')
        scope = _scope(payload.get('scope'))
        label = targets[target_key].main_id or target_key
        reason = str(payload.get('reason', '')).strip()
        if not reason:
            reason = 'Return {} to inherited automatic model policy'.format(
                label
            )
        job_id = self.state.reset_policy_and_queue(
            target_key, package=self.package, scope=scope, actor=actor,
            reason=reason,
        )
        self.actor = actor
        return {'job_id': job_id, 'status': 'queued'}

    def update_selection(self, payload):
        manager = FitManager(self.package)
        reference = str(payload.get('reference', ''))
        action = payload.get('action')
        candidate_id = manager.resolve(reference)
        if action == 'default':
            selection = manager.select(default=candidate_id)
        elif action == 'display':
            current = manager.library.read_selection()
            displayed = list(current.displayed)
            should_display = bool(payload.get('displayed'))
            if should_display and candidate_id not in displayed:
                displayed.append(candidate_id)
            elif not should_display:
                if candidate_id == current.default:
                    raise ValueError('the default fit must remain displayed')
                displayed = [
                    value for value in displayed if value != candidate_id
                ]
            selection = manager.select(displayed=displayed)
        else:
            raise ValueError('unsupported selection action')
        manager.render(dataset=self.dataset())
        return {
            'default': selection.default,
            'displayed': list(selection.displayed),
        }


def make_operator_server(package, *, state_path=None, actor=None,
                         host='127.0.0.1', port=0):
    if host not in {'127.0.0.1', 'localhost', '::1'}:
        raise ValueError('operator service must bind to a loopback address')
    package = Path(package).expanduser().resolve()
    state = default_operator_state(package, state_path)
    service = OperatorService(package, state, actor=actor)
    handler = partial(_OperatorHandler, service=service, directory=str(package))
    server = ThreadingHTTPServer((host, port), handler)
    server.operator_service = service
    return server


def serve_operator(package, *, state_path=None, actor=None,
                   host='127.0.0.1', port=8765):
    server = make_operator_server(
        package, state_path=state_path, actor=actor, host=host, port=port,
    )
    address, actual_port = server.server_address[:2]
    print('SDF operator page: http://{}:{}/'.format(address, actual_port))
    print('Operator state: {}'.format(server.operator_service.state.path))
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


class _OperatorHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, service, **kwargs):
        self.service = service
        super().__init__(*args, **kwargs)

    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path in {'/', '/operator.html'}:
            fits = parse_qs(parsed.query).get('fit')
            try:
                html = self.service.render_page(fits=fits)
            except Exception as error:
                self._json_error(error, status=400)
                return
            self._send(
                200, html.encode('utf-8'), 'text/html; charset=utf-8'
            )
            return
        if parsed.path == '/favicon.ico':
            self._send(204, b'', 'image/x-icon')
            return
        super().do_GET()

    def do_POST(self):
        if self.headers.get('X-SDF-Token') != self.service.token:
            self._json_error(ValueError('invalid operator token'), status=403)
            return
        try:
            length = int(self.headers.get('Content-Length', '0'))
            if length <= 0 or length > 65536:
                raise ValueError('invalid request size')
            payload = json.loads(self.rfile.read(length))
            path = urlparse(self.path).path
            if path == '/api/policy':
                result = self.service.save_policy(payload)
            elif path == '/api/system-policy':
                result = self.service.save_system(payload)
            elif path == '/api/policy/reset':
                result = self.service.reset_policy(payload)
            elif path == '/api/selection':
                result = self.service.update_selection(payload)
            else:
                self._json_error(ValueError('unknown operator endpoint'), 404)
                return
        except Exception as error:
            self._json_error(error, status=400)
            return
        self._send(
            200, json.dumps(result).encode('utf-8'),
            'application/json; charset=utf-8',
        )

    def end_headers(self):
        self.send_header('Cache-Control', 'no-store')
        super().end_headers()

    def _json_error(self, error, status):
        self._send(
            status, json.dumps({'error': str(error)}).encode('utf-8'),
            'application/json; charset=utf-8',
        )

    def _send(self, status, body, content_type):
        self.send_response(status)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def _scope(value):
    value = str(value or 'global')
    if value not in {'global', 'package'}:
        raise ValueError('policy scope must be global or package')
    return value


def _policy_reason(label, policy):
    if policy.max_disks:
        return 'Use {} for {}; allow up to {} {} disk component(s)'.format(
            policy.stellar_models[0], label, policy.max_disks,
            policy.disk_models[0],
        )
    return 'Use {} for {}; do not fit a disk'.format(
        policy.stellar_models[0], label,
    )
