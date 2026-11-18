"""Durable SDB target acquisition behind the catalog web service."""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import subprocess
import tempfile
from threading import Event, Thread


@dataclass(frozen=True)
class SdbImportConfig:
    """Explicit, non-shell command configuration for one SDB database."""

    executable: str
    database: Path
    reference_database: Path | None = None
    cache_database: Path | None = None
    config: Path | None = None
    timeout_seconds: float = 1800.0
    sample_timeout_seconds: float = 43200.0

    @classmethod
    def create(cls, database, *, executable='sdb', reference_database=None,
               cache_database=None, config=None, timeout_seconds=1800.0,
               sample_timeout_seconds=43200.0):
        database = Path(database).expanduser().resolve()
        return cls(
            executable=str(executable), database=database,
            reference_database=_optional_path(reference_database),
            cache_database=_optional_path(cache_database),
            config=_optional_path(config),
            timeout_seconds=float(timeout_seconds),
            sample_timeout_seconds=float(sample_timeout_seconds),
        )

    def command(self):
        command = [
            self.executable, '--database', str(self.database),
            '--compact-json', '--quiet',
        ]
        if self.reference_database is not None:
            command[1:1] = [
                '--reference-database', str(self.reference_database),
            ]
        if self.cache_database is not None:
            command[1:1] = [
                '--cache-database', str(self.cache_database),
            ]
        if self.config is not None:
            command[1:1] = ['--config', str(self.config)]
        return command


class SdbImporter:
    """Identify, ensure, and export one exact SIMBAD name through SDB."""

    def __init__(self, config, export_root, *, runner=None):
        self.config = config
        self.export_root = Path(export_root).expanduser().resolve()
        self.runner = runner or subprocess.run

    def acquire(self, query, *, actor):
        query = ' '.join(str(query or '').split())
        if not query:
            raise ValueError('a SIMBAD identifier is required')
        environment = os.environ.copy()
        environment['SDB_ACTOR'] = str(actor)
        added, add_status, add_stderr = self._json_command(
            [*self.config.command(), 'add', query, '--ensure'],
            environment=environment,
        )
        items = added.get('items', ())
        sdbids = tuple(dict.fromkeys(
            str(item['sdbid']) for item in items
            if isinstance(item, dict) and item.get('sdbid')
            and item.get('status') != 'failed'
        ))
        if not sdbids:
            detail = next((
                str(item.get('error')) for item in items
                if isinstance(item, dict) and item.get('error')
            ), None)
            raise RuntimeError(detail or _command_error(
                'SDB could not identify the target', add_status, add_stderr
            ))

        packages = []
        manifests = []
        exported = 0
        skipped = 0
        for sdbid in sdbids:
            summary, status, stderr = self._json_command(
                [
                    *self.config.command(), 'export', sdbid,
                    '--output-dir', str(self.export_root), '--workers', '1',
                ],
                environment=environment,
            )
            if status != 0 or int(summary.get('failed', 0)):
                raise RuntimeError(_command_error(
                    'SDB export failed', status, stderr
                ))
            manifest = Path(summary['manifest']).expanduser().resolve()
            try:
                manifest.relative_to(self.export_root)
            except ValueError as error:
                raise RuntimeError(
                    'SDB export manifest lies outside the served root'
                ) from error
            payload = json.loads(manifest.read_text(encoding='utf-8'))
            for package in payload.get('packages', ()):
                if package.get('status') == 'completed' and package.get('directory'):
                    packages.append(str(package['directory']))
            manifests.append(manifest.name)
            exported += int(summary.get('exported', 0))
            skipped += int(summary.get('skipped', 0))

        packages = tuple(dict.fromkeys(packages))
        if not packages:
            raise RuntimeError('SDB export produced no complete fitting package')
        update = added.get('update_summary') or {}
        warnings = []
        if add_status != 0:
            warnings.append(
                '{} provider update(s) failed; the available data were exported'.
                format(update.get('failed', 'Some'))
            )
        return {
            'query': query,
            'sdbids': list(sdbids),
            'packages': list(packages),
            'manifests': manifests,
            'created': int(added.get('created_count', 0)),
            'existing': int(added.get('existing_count', 0)),
            'providers': list(added.get('providers', ())),
            'exported': exported,
            'skipped': skipped,
            'warnings': warnings,
        }

    def acquire_sample(self, job, *, progress=None):
        """Create one owned SDB sample, import its targets, and export it."""

        if getattr(job, 'operation', 'add') == 'remove':
            return self._remove_sample_members(job)

        progress = progress or (lambda event: None)
        progress({
            'stage': 'starting',
            'description': 'Starting sample import',
            'current': 0,
            'total': len(job.targets),
        })
        environment = os.environ.copy()
        environment['SDB_ACTOR'] = str(job.actor)
        with tempfile.NamedTemporaryFile(
                mode='w', encoding='utf-8', suffix='.txt') as target_file:
            target_file.write('\n'.join(job.targets) + '\n')
            target_file.flush()
            command = [
                *self.config.command(), 'sample', 'acquire', job.sample_name,
                target_file.name, '--visibility', job.visibility,
                '--actor', job.actor,
                '--reason', 'Imported from the web sample page',
            ]
            if getattr(job, 'owner_username', None):
                command.extend(('--owner', job.owner_username))
            if job.note:
                command.extend(('--note', job.note))
            acquired, status, stderr = self._json_command_with_progress(
                command, environment=environment, progress=progress,
                timeout_seconds=self.config.sample_timeout_seconds,
            )
        # Exit 1 means that useful work completed but one or more targets or
        # optional provider updates failed. Exit 2 is a command-level failure.
        if status not in {0, 1}:
            raise RuntimeError(_command_error(
                'SDB sample import failed', status, stderr
            ))

        warnings = []
        failed_count = int(acquired.get('failed_count', 0))
        update = (acquired.get('import_result') or {}).get(
            'update_summary'
        ) or {}
        if failed_count:
            warnings.append(
                '{} of {} target identifier(s) could not be resolved'.format(
                    failed_count,
                    acquired.get('requested_count', len(job.targets)),
                )
            )
        if int(update.get('failed', 0)):
            warnings.append(
                '{} provider update(s) failed; successful targets were '
                'retained'.format(update['failed'])
            )
        result = {
            'sample': job.sample_name,
            'operation': 'add',
            'requested_count': len(job.targets),
            'succeeded_count': int(acquired.get('succeeded_count', 0)),
            'failed_count': failed_count,
            'packages': [],
            'manifests': [],
            'exported': 0,
            'skipped': 0,
            'warnings': warnings,
            'import_result': acquired,
        }
        if not result['succeeded_count']:
            return result

        progress({
            'stage': 'exporting',
            'description': 'Exporting successful targets',
            'current': None,
            'total': None,
        })
        exported, export_status, export_stderr = self._json_command(
            [
                *self.config.command(), 'export', '--sample', job.sample_name,
                '--output-dir', str(self.export_root), '--workers', '1',
            ],
            environment=environment,
            timeout_seconds=self.config.sample_timeout_seconds,
        )
        if export_status != 0 or int(exported.get('failed', 0)):
            warnings.append(
                _command_error(
                    'Targets were imported but SDB export failed',
                    export_status, export_stderr,
                )
            )
            if not exported.get('manifest'):
                return result
        manifest, packages = self._manifest_packages(
            exported, allow_empty=True
        )
        result.update(
            packages=list(packages), manifests=[manifest.name],
            exported=int(exported.get('exported', 0)),
            skipped=int(exported.get('skipped', 0)),
        )
        if not packages:
            warnings.append(
                'Targets were imported, but none is ready for a fitting package'
            )
        return result

    def _json_command_with_progress(
            self, command, *, environment, progress, timeout_seconds=None):
        """Run an SDB command while forwarding its JSONL progress sidecar."""

        with tempfile.NamedTemporaryFile(
                mode='w', encoding='utf-8', suffix='.jsonl') as event_file:
            child_environment = dict(environment)
            child_environment['SDB_PROGRESS_FILE'] = event_file.name
            stop = Event()
            watcher = Thread(
                target=_watch_progress_file,
                args=(event_file.name, progress, stop),
                name='sdb-import-progress', daemon=True,
            )
            watcher.start()
            try:
                return self._json_command(
                    command, environment=child_environment,
                    timeout_seconds=timeout_seconds,
                )
            finally:
                stop.set()
                watcher.join(timeout=2)

    def _remove_sample_members(self, job):
        """Remove selected members and publish a fresh sample snapshot."""

        environment = os.environ.copy()
        environment['SDB_ACTOR'] = str(job.actor)
        removed = []
        failed = []
        for target in job.targets:
            try:
                value, status, stderr = self._json_command(
                    [
                        *self.config.command(), 'sample', 'remove',
                        job.sample_name, target, '--actor', job.actor,
                        '--reason', 'Removed from the web sample page',
                    ],
                    environment=environment,
                )
            except RuntimeError as error:
                failed.append({'target': target, 'error': str(error)})
                continue
            if status != 0:
                failed.append({
                    'target': target,
                    'error': _command_error(
                        'SDB could not remove the target', status, stderr,
                    ),
                })
                continue
            removed.append({
                'target': target,
                'action_id': value.get('action_id'),
            })

        result = {
            'sample': job.sample_name,
            'operation': 'remove',
            'requested_count': len(job.targets),
            'succeeded_count': len(removed),
            'failed_count': len(failed),
            'member_results': removed + failed,
            'packages': [],
            'manifests': [],
            'exported': 0,
            'skipped': 0,
            'warnings': [],
        }
        if failed:
            result['warnings'].append(
                '{} target removal(s) failed'.format(len(failed))
            )

        empty, empty_status, empty_stderr = self._json_command(
            [
                *self.config.command(), 'sample', 'delete-if-empty',
                job.sample_name, '--actor', job.actor,
                '--reason', 'Deleted after its final member was removed '
                'from the web sample page',
            ],
            environment=environment,
        )
        if empty_status != 0:
            raise RuntimeError(_command_error(
                'Membership changed but the empty-sample check failed',
                empty_status, empty_stderr,
            ))
        result['sample_deleted'] = bool(empty.get('deleted'))
        result['remaining_member_count'] = int(empty.get('member_count', 0))

        if result['sample_deleted']:
            catalog, catalog_status, catalog_stderr = self._json_command(
                [
                    *self.config.command(), 'sample', 'catalog',
                    '--output-dir', str(self.export_root),
                ],
                environment=environment,
            )
            if catalog_status != 0:
                raise RuntimeError(_command_error(
                    'Sample was deleted but the served sample list could not '
                    'be refreshed', catalog_status, catalog_stderr,
                ))
            result['sample_catalog'] = Path(catalog['path']).name
            return result

        exported, export_status, export_stderr = self._json_command(
            [
                *self.config.command(), 'export', '--sample', job.sample_name,
                '--output-dir', str(self.export_root), '--workers', '1',
            ],
            environment=environment,
        )
        if export_status != 0 or int(exported.get('failed', 0)):
            raise RuntimeError(_command_error(
                'Membership changed but SDB export failed',
                export_status, export_stderr,
            ))
        manifest, packages = self._manifest_packages(
            exported, allow_empty=True,
        )
        result.update(
            packages=list(packages), manifests=[manifest.name],
            exported=int(exported.get('exported', 0)),
            skipped=int(exported.get('skipped', 0)),
        )
        return result

    def check_sample_targets(self, targets):
        """Classify target identifiers using only the local SDB database."""

        targets = tuple(dict.fromkeys(
            ' '.join(str(value or '').split()) for value in targets
            if str(value or '').strip()
        ))
        if not targets:
            raise ValueError('at least one target identifier is required')
        if len(targets) > 1000:
            raise ValueError('sample checks are limited to 1000 targets')
        with tempfile.NamedTemporaryFile(
                mode='w', encoding='utf-8', suffix='.txt') as target_file:
            target_file.write('\n'.join(targets) + '\n')
            target_file.flush()
            checked, status, stderr = self._json_command(
                [
                    *self.config.command(), 'sample', 'check',
                    target_file.name,
                ],
                environment=os.environ.copy(),
            )
        if status not in {0, 1}:
            raise RuntimeError(_command_error(
                'SDB sample check failed', status, stderr
            ))
        return checked

    def update_sample_collaborator(self, sample_name, username, *, operation,
                                   actor):
        """Add or remove one collaborator through the authoritative SDB CLI."""

        operation = str(operation or '').strip().lower()
        if operation not in {'share', 'unshare'}:
            raise ValueError(
                'sample collaborator operation must be share or unshare'
            )
        sample_name = str(sample_name or '').strip()
        username = str(username or '').strip()
        if not sample_name:
            raise ValueError('sample name is required')
        if not username:
            raise ValueError('username is required')
        environment = os.environ.copy()
        environment['SDB_ACTOR'] = str(actor)
        result, status, stderr = self._json_command(
            [
                *self.config.command(), 'sample', operation,
                sample_name, username,
            ],
            environment=environment,
        )
        if status != 0:
            raise RuntimeError(_command_error(
                'SDB could not update sample collaborators', status, stderr,
            ))
        return result

    def _manifest_packages(self, summary, *, allow_empty=False):
        manifest = Path(summary['manifest']).expanduser().resolve()
        try:
            manifest.relative_to(self.export_root)
        except ValueError as error:
            raise RuntimeError(
                'SDB export manifest lies outside the served root'
            ) from error
        payload = json.loads(manifest.read_text(encoding='utf-8'))
        packages = tuple(dict.fromkeys(
            str(package['directory'])
            for package in payload.get('packages', ())
            if package.get('status') == 'completed' and package.get('directory')
        ))
        if not packages and not allow_empty:
            raise RuntimeError('SDB export produced no complete fitting package')
        return manifest, packages

    def _json_command(self, command, *, environment, timeout_seconds=None):
        timeout_seconds = (
            self.config.timeout_seconds
            if timeout_seconds is None
            else float(timeout_seconds)
        )
        try:
            completed = self.runner(
                command, check=False, capture_output=True, text=True,
                timeout=timeout_seconds, env=environment,
            )
        except (OSError, subprocess.SubprocessError) as error:
            raise RuntimeError('could not run SDB: {}'.format(error)) from error
        stdout = str(completed.stdout or '').strip()
        try:
            payload = json.loads(stdout)
        except json.JSONDecodeError as error:
            raise RuntimeError(_command_error(
                'SDB returned invalid output', completed.returncode,
                completed.stderr or stdout,
            )) from error
        return payload, int(completed.returncode), str(completed.stderr or '')


def run_acquisition_worker(state, importer, *, on_export=None,
                           on_sample_export=None,
                           poll_seconds=5.0, stop_event=None,
                           worker_id=None):
    """Continuously drain the single-writer SDB acquisition queue.

    SDB takes the export-root file lock itself, so callers must not wrap the
    acquisition in the catalog synchronization lock.  The completion callback
    may take that shorter-lived lock while it indexes the finished package.
    """

    worker_id = worker_id or 'sdb-import-{}'.format(os.getpid())
    while stop_event is None or not stop_event.is_set():
        job = state.claim_next_acquisition(worker_id=worker_id)
        sample_job = None
        if job is None:
            sample_job = state.claim_next_sample_acquisition(
                worker_id=worker_id
            )
        if job is None and sample_job is None:
            if stop_event is None:
                import time
                time.sleep(poll_seconds)
            else:
                stop_event.wait(poll_seconds)
            continue
        try:
            if sample_job is not None:
                progress = lambda event: state.update_sample_acquisition_progress(
                    sample_job.id, event,
                )
                result = importer.acquire_sample(
                    sample_job, progress=progress,
                )
                callback = on_sample_export or on_export
                actor = sample_job.actor
            else:
                result = importer.acquire(job.query, actor=job.actor)
                callback = on_export
                actor = job.actor
            if callback is not None:
                if sample_job is not None:
                    progress({
                        'stage': 'publishing',
                        'description': 'Publishing the sample and queueing fits',
                        'current': None,
                        'total': None,
                    })
                enriched = callback(result, actor=actor)
                if enriched is not None:
                    result = enriched
        except Exception as error:
            message = '{}: {}'.format(type(error).__name__, error)
            if sample_job is not None:
                state.finish_sample_acquisition(sample_job.id, error=message)
            else:
                state.finish_acquisition(job.id, error=message)
        else:
            if sample_job is not None:
                state.finish_sample_acquisition(sample_job.id, result=result)
            else:
                state.finish_acquisition(job.id, result=result)


def _watch_progress_file(path, callback, stop, poll_seconds=0.2):
    """Forward the latest complete event from an append-only JSONL file."""

    position = 0
    remainder = ''

    def drain():
        nonlocal position, remainder
        try:
            with open(path, encoding='utf-8') as handle:
                handle.seek(position)
                chunk = handle.read()
                position = handle.tell()
        except FileNotFoundError:
            return
        if not chunk:
            return
        remainder += chunk
        lines = remainder.split('\n')
        remainder = lines.pop()
        latest = None
        for line in lines:
            try:
                latest = json.loads(line)
            except json.JSONDecodeError:
                continue
        if latest is not None:
            try:
                callback(latest)
            except Exception:
                pass

    while not stop.wait(poll_seconds):
        drain()
    drain()


def _optional_path(value):
    return None if value is None else Path(value).expanduser().resolve()


def _command_error(message, status, stderr):
    detail = str(stderr or '').strip().splitlines()
    suffix = '' if not detail else ': ' + detail[-1][:500]
    return '{} (exit {}){}'.format(message, status, suffix)
