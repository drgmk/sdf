"""Durable SDB target acquisition behind the catalog web service."""

from __future__ import annotations

from dataclasses import dataclass
from contextlib import nullcontext
import json
import os
from pathlib import Path
import subprocess


@dataclass(frozen=True)
class SdbImportConfig:
    """Explicit, non-shell command configuration for one SDB database."""

    executable: str
    database: Path
    reference_database: Path | None = None
    cache_database: Path | None = None
    config: Path | None = None
    timeout_seconds: float = 1800.0

    @classmethod
    def create(cls, database, *, executable='sdb', reference_database=None,
               cache_database=None, config=None, timeout_seconds=1800.0):
        database = Path(database).expanduser().resolve()
        return cls(
            executable=str(executable), database=database,
            reference_database=_optional_path(reference_database),
            cache_database=_optional_path(cache_database),
            config=_optional_path(config),
            timeout_seconds=float(timeout_seconds),
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

    def _json_command(self, command, *, environment):
        try:
            completed = self.runner(
                command, check=False, capture_output=True, text=True,
                timeout=self.config.timeout_seconds, env=environment,
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
                           poll_seconds=5.0, stop_event=None,
                           worker_id=None, work_lock=None):
    """Continuously drain the single-writer SDB acquisition queue."""

    worker_id = worker_id or 'sdb-import-{}'.format(os.getpid())
    while stop_event is None or not stop_event.is_set():
        job = state.claim_next_acquisition(worker_id=worker_id)
        if job is None:
            if stop_event is None:
                import time
                time.sleep(poll_seconds)
            else:
                stop_event.wait(poll_seconds)
            continue
        try:
            with work_lock if work_lock is not None else nullcontext():
                result = importer.acquire(job.query, actor=job.actor)
                if on_export is not None:
                    enriched = on_export(result, actor=job.actor)
                    if enriched is not None:
                        result = enriched
        except Exception as error:
            state.finish_acquisition(
                job.id, error='{}: {}'.format(type(error).__name__, error)
            )
        else:
            state.finish_acquisition(job.id, result=result)


def _optional_path(value):
    return None if value is None else Path(value).expanduser().resolve()


def _command_error(message, status, stderr):
    detail = str(stderr or '').strip().splitlines()
    suffix = '' if not detail else ': ' + detail[-1][:500]
    return '{} (exit {}){}'.format(message, status, suffix)
