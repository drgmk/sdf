"""Incrementally synchronize changed export packages with the SDF catalog."""

from __future__ import annotations

from contextlib import contextmanager, nullcontext
import fcntl
import json
from time import monotonic

from .result_index import discover_packages


class CatalogSynchronizer:
    """Poll cheap input metadata and refresh only changed package rows."""

    def __init__(self, service):
        self.service = service
        self.root = service.root
        self._signatures = self._scan_signatures()
        self._marker_signature, self._generation = self._read_generation()

    def sync(self, *, queue_changes=True, autofit_existing=False,
             full_reconcile=False):
        """Refresh changed packages and optionally queue fitting work."""

        # Keep this order consistent with acquisition completion: taking the
        # file lock first can deadlock an SDB export waiting for exclusivity.
        lock = getattr(self.service, 'export_sync_lock', None)
        with lock if lock is not None else nullcontext():
            with _available_export_root(self.root) as available:
                if not available:
                    return {
                        'added': (), 'changed': (), 'removed': (),
                        'queued': (), 'export_busy': True,
                    }
                result = self._sync(
                    queue_changes=queue_changes,
                    autofit_existing=autofit_existing,
                    full_reconcile=full_reconcile,
                )
                result['export_busy'] = False
                return result

    def _sync(self, *, queue_changes=True, autofit_existing=False,
              full_reconcile=False):
        """Run one synchronization while the shared export lock is held."""

        marker_signature, generation = self._read_generation()
        previous = self._signatures
        if (
            not full_reconcile
            and marker_signature is not None
            and marker_signature == self._marker_signature
        ):
            return {
                'added': (), 'changed': (), 'removed': (), 'queued': (),
            }
        changed_paths = None
        if not full_reconcile and generation is not None:
            changed_paths = self._manifest_packages_since(
                self._generation, generation
            )
        if changed_paths is None:
            current = self._scan_signatures()
            added = sorted(set(current) - set(previous))
            changed = sorted(
                package for package in set(current) & set(previous)
                if current[package] != previous[package]
            )
            removed = sorted(set(previous) - set(current))
        else:
            current = dict(previous)
            added = []
            changed = []
            removed = []
            for package_path in sorted(changed_paths):
                package = self.root / package_path
                if not package.is_dir():
                    if package_path in current:
                        current.pop(package_path)
                        removed.append(package_path)
                    continue
                signature = _input_signature(package)
                if package_path not in previous:
                    added.append(package_path)
                elif previous[package_path] != signature:
                    changed.append(package_path)
                current[package_path] = signature
        for package_path in added + changed:
            self.service.index.update_package(self.root / package_path)
        for package_path in removed:
            self.service.index.remove_package(package_path)
        if autofit_existing:
            candidates = self.service.packages()
        elif queue_changes:
            candidates = tuple(
                self.service.package(package_path)
                for package_path in added + changed
            )
        else:
            candidates = ()
        active = {
            job.package_path for job in self.service.state.jobs()
            if job.status in {'queued', 'running'}
        }
        queued = []
        for row in candidates:
            package = self.root / row['package_path']
            if row['status'] not in {'unfit', 'stale'}:
                continue
            if str(package) in active:
                continue
            result = self.service.queue_package(
                row['package_path'], actor='sdf watcher',
                reason='Fit new or changed SDB export {}'.format(
                    ' / '.join(row['target_labels']) or row['package_path']
                ),
                remember_actor=False,
            )
            queued.append(result['job_id'])
            active.add(str(package))
        self._signatures = current
        self._marker_signature = marker_signature
        self._generation = generation
        return {
            'added': tuple(added), 'changed': tuple(changed),
            'removed': tuple(removed), 'queued': tuple(queued),
        }

    def run(self, stop_event, *, poll_seconds=5.0,
            autofit_existing=False, reconcile_seconds=3600.0):
        if autofit_existing:
            self.sync(queue_changes=False, autofit_existing=True)
        next_reconcile = monotonic() + reconcile_seconds
        while not stop_event.wait(poll_seconds):
            now = monotonic()
            self.sync(
                queue_changes=True,
                full_reconcile=now >= next_reconcile,
            )
            if now >= next_reconcile:
                next_reconcile = now + reconcile_seconds

    def _scan_signatures(self):
        return {
            package.relative_to(self.root).as_posix(): _input_signature(package)
            for package in discover_packages(self.root)
        }

    def _read_generation(self):
        path = self.root / '.sdb-export-generation.json'
        try:
            stat = path.stat()
            payload = json.loads(path.read_text(encoding='utf-8'))
            if payload.get('schema') != 'sdb-export-generation':
                return None, None
            return (
                (stat.st_mtime_ns, stat.st_size),
                int(payload['generation']),
            )
        except (FileNotFoundError, KeyError, TypeError, ValueError,
                json.JSONDecodeError, OSError):
            return None, None

    def _manifest_packages_since(self, previous, current):
        if previous is not None and current <= previous:
            return set()
        manifests = []
        for path in self.root.glob('export-*-manifest.json'):
            try:
                run_id = int(path.name[len('export-'):-len('-manifest.json')])
            except ValueError:
                continue
            if (previous is None or run_id > previous) and run_id <= current:
                manifests.append((run_id, path))
        if not manifests or max(run_id for run_id, _path in manifests) != current:
            return None
        packages = set()
        try:
            for _run_id, path in sorted(manifests):
                payload = json.loads(path.read_text(encoding='utf-8'))
                if payload.get('schema') != 'sdb-fit-package-export':
                    return None
                packages.update(
                    str(row['directory']) for row in payload.get('packages', ())
                    if row.get('status') == 'completed'
                )
        except (KeyError, TypeError, ValueError, json.JSONDecodeError, OSError):
            return None
        return packages


def _input_signature(package):
    paths = set(package.glob('*-rawphot.txt'))
    for name in ('identities.yml', 'joint-fit.yml', 'spectra.yml'):
        manifest = package / name
        if manifest.is_file():
            paths.add(manifest)
    spectra = package / 'spectra'
    if spectra.is_dir():
        paths.update(path for path in spectra.rglob('*') if path.is_file())
    return tuple(
        (
            path.relative_to(package).as_posix(),
            path.stat().st_size,
            path.stat().st_mtime_ns,
        )
        for path in sorted(
            paths, key=lambda value: value.relative_to(package).as_posix()
        )
    )


@contextmanager
def _available_export_root(root):
    """Skip a watcher pass while SDB is changing package inputs."""

    path = root / '.sdb-export.lock'
    with path.open('a+', encoding='utf-8') as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_SH | fcntl.LOCK_NB)
        except BlockingIOError:
            yield False
            return
        try:
            yield True
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
