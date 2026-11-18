import fcntl
import json
import os
import time

from .context import sdf
from sdf.work_cleanup import cleanup_sampler_work


DAY = 86400


def _work(package, name, *, age_days, now, complete=True):
    directory = package / 'fits' / '.work' / name
    directory.mkdir(parents=True)
    (directory / 'run_resume.dat').write_text('resume\n', encoding='utf-8')
    if complete:
        (directory / 'run_stats.dat').write_text('stats\n', encoding='utf-8')
        (directory / 'run_post_equal_weights.dat').write_text(
            'posterior\n', encoding='utf-8'
        )
    modified = now - age_days * DAY
    for path in directory.iterdir():
        os.utime(path, (modified, modified))
    os.utime(directory, (modified, modified))
    return directory


def _reference(package, work):
    result = package / 'fits' / 'model--candidate' / 'result.json'
    result.parent.mkdir(parents=True)
    result.write_text(json.dumps({
        'posterior': {
            'native_output': str(
                (work / 'run_').relative_to(package)
            ),
        },
    }), encoding='utf-8')


def test_cleanup_classifies_retains_and_removes_sampler_work(tmp_path):
    now = 2_000_000_000.0
    package = tmp_path / 'target'
    committed = _work(
        package, 'committed-input-fit', age_days=31, now=now
    )
    superseded = _work(
        package, 'superseded-input-fit', age_days=8, now=now
    )
    interrupted = _work(
        package, 'interrupted-input-fit', age_days=31, now=now,
        complete=False,
    )
    recent = _work(package, 'recent-input-fit', age_days=2, now=now)
    _reference(package, committed)

    report = cleanup_sampler_work(tmp_path, now=now)

    decisions = {item.path.name: item for item in report.items}
    assert decisions[committed.name].category == 'committed'
    assert decisions[superseded.name].category == 'superseded'
    assert decisions[interrupted.name].category == 'interrupted'
    assert decisions[recent.name].category == 'superseded'
    assert {item.path.name for item in report.items if item.eligible} == {
        committed.name, superseded.name, interrupted.name,
    }
    assert all(path.is_dir() for path in (
        committed, superseded, interrupted, recent
    ))

    applied = cleanup_sampler_work(tmp_path, now=now, apply=True)

    assert {item.path.name for item in applied.items if item.removed} == {
        committed.name, superseded.name, interrupted.name,
    }
    assert not committed.exists()
    assert not superseded.exists()
    assert not interrupted.exists()
    assert recent.is_dir()


def test_cleanup_protects_held_locks_and_active_packages(tmp_path):
    now = 2_000_000_000.0
    locked_package = tmp_path / 'locked'
    locked = _work(locked_package, 'old', age_days=40, now=now)
    lock_path = locked / '.fit.lock'
    lock_path.touch()
    active_package = tmp_path / 'active'
    active = _work(active_package, 'old', age_days=40, now=now)

    with lock_path.open('r+') as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        report = cleanup_sampler_work(
            tmp_path, now=now, apply=True,
            active_packages=(active_package,),
        )

    decisions = {item.package.name: item for item in report.items}
    assert decisions['locked'].reason == 'fit lock is held'
    assert decisions['active'].reason == (
        'package has a queued or running fit job'
    )
    assert locked.is_dir()
    assert active.is_dir()


def test_cleanup_ignores_symlinked_work_directories(tmp_path):
    outside = tmp_path / 'outside'
    outside.mkdir()
    work = tmp_path / 'target' / 'fits' / '.work'
    work.mkdir(parents=True)
    (work / 'linked').symlink_to(outside, target_is_directory=True)

    report = cleanup_sampler_work(tmp_path, now=2_000_000_000.0, apply=True)

    assert report.items == ()
    assert outside.is_dir()


def test_cleanup_protects_package_with_unreadable_result(tmp_path):
    now = 2_000_000_000.0
    package = tmp_path / 'target'
    work = _work(package, 'old', age_days=40, now=now)
    result = package / 'fits' / 'model--candidate' / 'result.json'
    result.parent.mkdir(parents=True)
    result.write_text('{not json', encoding='utf-8')

    report = cleanup_sampler_work(tmp_path, now=now, apply=True)

    assert report.items[0].reason == (
        'package has unreadable retained result metadata'
    )
    assert work.is_dir()


def test_cleanup_lock_creation_does_not_refresh_work_age(tmp_path):
    now = time.time()
    package = tmp_path / 'target'
    work = _work(package, 'old', age_days=8, now=now)

    report = cleanup_sampler_work(tmp_path, now=now, apply=True)

    assert report.items[0].eligible
    assert report.items[0].removed
    assert report.items[0].reason == (
        'removed past 7-day superseded retention'
    )
    assert not work.exists()


def test_cleanup_recovers_work_refreshed_by_an_earlier_lock(tmp_path):
    now = time.time()
    package = tmp_path / 'target'
    work = _work(package, 'old', age_days=8, now=now)

    # The original cleanup bug left this exact state behind: old sampler
    # content plus a newly-created lock and newly-refreshed directory mtime.
    (work / '.fit.lock').touch()
    os.utime(work, (now, now))

    report = cleanup_sampler_work(tmp_path, now=now, apply=True)

    assert report.items[0].age_days >= 8
    assert report.items[0].eligible
    assert report.items[0].removed
    assert not work.exists()
