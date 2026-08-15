"""Incrementally synchronize changed export packages with the SDF catalog."""

from __future__ import annotations

from .result_index import discover_packages


class CatalogSynchronizer:
    """Poll cheap input metadata and refresh only changed package rows."""

    def __init__(self, service):
        self.service = service
        self.root = service.root
        self._signatures = self._scan_signatures()

    def sync(self, *, queue_changes=True, autofit_existing=False):
        """Refresh changed packages and optionally queue fitting work."""

        current = self._scan_signatures()
        previous = self._signatures
        added = sorted(set(current) - set(previous))
        changed = sorted(
            package for package in set(current) & set(previous)
            if current[package] != previous[package]
        )
        removed = sorted(set(previous) - set(current))
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
        return {
            'added': tuple(added), 'changed': tuple(changed),
            'removed': tuple(removed), 'queued': tuple(queued),
        }

    def run(self, stop_event, *, poll_seconds=5.0,
            autofit_existing=False):
        if autofit_existing:
            self.sync(queue_changes=False, autofit_existing=True)
        while not stop_event.wait(poll_seconds):
            self.sync(queue_changes=True)

    def _scan_signatures(self):
        return {
            package.relative_to(self.root).as_posix(): _input_signature(package)
            for package in discover_packages(self.root)
        }


def _input_signature(package):
    paths = set(package.glob('*-rawphot.txt'))
    topology = package / 'joint-fit.yml'
    if topology.is_file():
        paths.add(topology)
    return tuple(
        (path.name, path.stat().st_size, path.stat().st_mtime_ns)
        for path in sorted(paths, key=lambda value: value.name)
    )
