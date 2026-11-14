"""Reviewed retention and cleanup for sampler-native work directories."""

from __future__ import annotations

from dataclasses import dataclass
import fcntl
import json
import os
from pathlib import Path
import shutil
import time


@dataclass(frozen=True)
class WorkItem:
    """One sampler work directory and its cleanup decision."""

    path: Path
    package: Path
    category: str
    age_days: float
    size_bytes: int
    eligible: bool
    reason: str
    removed: bool = False


@dataclass(frozen=True)
class CleanupReport:
    """Result of scanning, or scanning and cleaning, an export root."""

    root: Path
    apply: bool
    items: tuple[WorkItem, ...]

    @property
    def total_bytes(self):
        return sum(item.size_bytes for item in self.items)

    @property
    def eligible_bytes(self):
        return sum(item.size_bytes for item in self.items if item.eligible)

    @property
    def removed_bytes(self):
        return sum(item.size_bytes for item in self.items if item.removed)


def cleanup_sampler_work(
        root, *, superseded_days=7, committed_days=30,
        interrupted_days=30, apply=False, active_packages=(), now=None):
    """Inspect sampler work and optionally remove entries past retention.

    Completed compact candidate artifacts do not depend on their native
    MultiNest files.  A grace period is retained for diagnosis, while
    interrupted output is retained separately because it may be resumed.
    """

    root = Path(root).expanduser().resolve()
    if not root.is_dir():
        raise ValueError('export root does not exist: {}'.format(root))
    thresholds = {
        'superseded': _days(superseded_days, 'superseded_days'),
        'committed': _days(committed_days, 'committed_days'),
        'interrupted': _days(interrupted_days, 'interrupted_days'),
    }
    active_packages = {
        Path(package).expanduser().resolve() for package in active_packages
    }
    now = time.time() if now is None else float(now)
    items = []
    references = {}
    for package, directory in _work_directories(root):
        if package not in references:
            references[package] = _referenced_work_directories(package)
        referenced, metadata_valid = references[package]
        category = _category(directory, referenced)
        size, newest = _size_and_newest(directory)
        age_days = max(0.0, (now - newest) / 86400.0)
        threshold = thresholds[category]
        eligible = age_days >= threshold
        reason = (
            'past {}-day {} retention'.format(threshold, category)
            if eligible else
            'within {}-day {} retention'.format(threshold, category)
        )
        if not metadata_valid:
            eligible = False
            reason = 'package has unreadable retained result metadata'
        elif package in active_packages:
            eligible = False
            reason = 'package has a queued or running fit job'
        elif _lock_is_held(directory / '.fit.lock'):
            eligible = False
            reason = 'fit lock is held'

        removed = False
        if apply and eligible:
            removed, reason = _remove_if_still_eligible(
                package, directory, thresholds, now
            )
        items.append(WorkItem(
            path=directory, package=package, category=category,
            age_days=age_days, size_bytes=size, eligible=eligible,
            reason=reason, removed=removed,
        ))
    return CleanupReport(root=root, apply=apply, items=tuple(items))


def _work_directories(root):
    for work in sorted(root.glob('*/fits/.work')):
        package = work.parent.parent.resolve()
        if not work.is_dir() or work.is_symlink():
            continue
        for directory in sorted(work.iterdir()):
            if directory.is_dir() and not directory.is_symlink():
                yield package, directory


def _category(directory, referenced):
    if directory.resolve() in referenced:
        return 'committed'
    if ((directory / 'run_stats.dat').is_file() and
            (directory / 'run_post_equal_weights.dat').is_file()):
        return 'superseded'
    return 'interrupted'


def _referenced_work_directories(package):
    referenced = set()
    valid = True
    for result_path in (package / 'fits').glob('*/result.json'):
        try:
            result = json.loads(result_path.read_text(encoding='utf-8'))
            native_output = result.get('posterior', {}).get('native_output')
            if native_output:
                referenced.add((package / native_output).resolve().parent)
        except (OSError, ValueError, TypeError):
            valid = False
    return referenced, valid


def _size_and_newest(directory):
    size = 0
    newest = directory.stat().st_mtime
    for base, directories, files in os.walk(directory):
        directories[:] = [
            name for name in directories
            if not Path(base, name).is_symlink()
        ]
        for name in files:
            if name == '.fit.lock':
                continue
            path = Path(base, name)
            if path.is_symlink():
                continue
            try:
                stat = path.stat()
            except FileNotFoundError:
                continue
            size += stat.st_size
            newest = max(newest, stat.st_mtime)
    return size, newest


def _lock_is_held(path):
    if not path.is_file():
        return False
    try:
        with path.open('r+') as handle:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                return True
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    except FileNotFoundError:
        return False
    except OSError:
        return True
    return False


def _remove_if_still_eligible(package, directory, thresholds, now):
    lock_path = directory / '.fit.lock'
    try:
        handle = lock_path.open('a+')
    except FileNotFoundError:
        return False, 'work directory disappeared before cleanup'
    try:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return False, 'fit lock became active before cleanup'
        referenced, metadata_valid = _referenced_work_directories(package)
        if not metadata_valid:
            return False, 'retained result metadata became unreadable'
        category = _category(directory, referenced)
        threshold = thresholds[category]
        _, newest = _size_and_newest(directory)
        age_days = max(0.0, (now - newest) / 86400.0)
        if age_days < threshold:
            return False, 'work changed while cleanup was scanning'
        try:
            shutil.rmtree(directory)
        except OSError as error:
            return False, 'could not remove work: {}'.format(error)
        return True, 'removed past {}-day {} retention'.format(
            threshold, category
        )
    finally:
        handle.close()


def _days(value, name):
    value = int(value)
    if value < 1:
        raise ValueError('{} must be at least one day'.format(name))
    return value
