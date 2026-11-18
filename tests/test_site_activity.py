from datetime import datetime, timedelta, timezone
import sqlite3
from zoneinfo import ZoneInfo

from .context import sdf
from sdf.site_activity import _access_projection, activity_report, job_health


def test_activity_report_combines_inventory_actions_and_filtered_traffic(
        tmp_path):
    state = sdf.operator_state.OperatorState(
        tmp_path / 'operator.sqlite'
    ).initialize()
    package = tmp_path / 'package'
    package.mkdir()
    state.queue_package(
        package, target_keys=('simbad:HD 123',), actor='grant@example.com',
        reason='Fit imported target',
    )

    sdb_path = tmp_path / 'sdb.sqlite'
    with sqlite3.connect(sdb_path) as connection:
        connection.executescript('''
            CREATE TABLE targets (
                id INTEGER PRIMARY KEY, sdbid TEXT, created_at TEXT
            );
            CREATE TABLE samples (
                id INTEGER PRIMARY KEY, name TEXT, visibility TEXT,
                owner_username TEXT, created_at TEXT, updated_at TEXT
            );
            CREATE TABLE users (
                username TEXT PRIMARY KEY, is_admin INTEGER, disabled INTEGER
            );
            CREATE TABLE normalized_measurements (id INTEGER PRIMARY KEY);
            CREATE TABLE spectral_assets (id INTEGER PRIMARY KEY);
            CREATE TABLE sample_membership_actions (
                id INTEGER PRIMARY KEY, sample_id INTEGER, target_id INTEGER,
                action TEXT, actor TEXT, reason TEXT, created_at TEXT
            );
            CREATE TABLE export_dirty_targets (
                id INTEGER PRIMARY KEY, target_id INTEGER, exported_at TEXT
            );
        ''')
        now = datetime.now(timezone.utc).isoformat()
        connection.execute(
            'INSERT INTO targets VALUES (1, ?, ?)', ('sdb-test', now)
        )
        connection.execute(
            'INSERT INTO samples VALUES (1, ?, ?, ?, ?, ?)',
            ('Nearby', 'private', 'grant@example.com', now, now),
        )
        connection.execute(
            'INSERT INTO users VALUES (?, 1, 0)', ('grant@example.com',)
        )
        connection.execute(
            'INSERT INTO normalized_measurements VALUES (1)'
        )
        connection.execute('INSERT INTO spectral_assets VALUES (1)')
        connection.execute(
            'INSERT INTO sample_membership_actions VALUES '
            '(1, 1, 1, ?, ?, ?, ?)',
            ('add', 'grant@example.com', 'Sample import', now),
        )
        connection.execute(
            'INSERT INTO export_dirty_targets VALUES (1, 1, ?)', (now,)
        )

    log = tmp_path / 'access.log'
    stamp = datetime.now().astimezone().strftime('%d/%b/%Y:%H:%M:%S %z')
    log.write_text('\n'.join((
        '127.0.0.1 - - [{0}] "GET /system/package/ HTTP/1.1" '
        '200 120 "-" "Mozilla/5.0" 12000'.format(stamp),
        '127.0.0.1 - - [{0}] "GET /samples/Nearby/ HTTP/1.1" '
        '200 120 "-" "Mozilla/5.0" 8000'.format(stamp),
        '127.0.0.1 - - [{0}] "GET /system/package/posterior.npz '
        'HTTP/1.1" 200 120 "-" "Mozilla/5.0" 20000'.format(stamp),
        '127.0.0.1 - - [{0}] "GET /system/package/ HTTP/1.1" '
        '200 120 "-" "ExampleBot/1" 4000'.format(stamp),
        '127.0.0.1 - - [{0}] "GET /favicon.ico HTTP/1.1" 204 0 '
        '"-" "Mozilla/5.0" 1000'.format(stamp),
    )), encoding='utf-8')

    report = activity_report(
        index_summary={
            'packages': 1, 'targets': 1, 'observations': 3,
            'candidates': 0, 'components': 0, 'parameters': 0,
        },
        packages=({'package_path': 'package', 'status': 'unfit',
                   'target_labels': ('HD 123',), 'sdbids': ('sdb-test',)},),
        state_path=state.path, sdb_path=sdb_path, access_log=log, days=30,
        display_timezone='Pacific/Auckland',
    )

    assert report['inventory']['database_targets'] == 1
    assert report['inventory']['private_samples'] == 1
    assert report['inventory']['measurements'] == 1
    assert report['changes']['targets_added'] == 1
    assert report['changes']['sample_edits'] == 1
    assert report['changes']['targets_exported'] == 1
    assert report['traffic']['page_views'] == 2
    assert report['traffic']['downloads'] == 1
    assert report['traffic']['top_systems'] == (
        {'label': 'HD 123', 'count': 1},
    )
    assert report['traffic']['top_samples'] == (
        {'label': 'Nearby', 'count': 1},
    )
    assert any(event['kind'] == 'Sample' for event in report['events'])
    assert report['timezone'] == 'Pacific/Auckland'
    assert report['events'][0]['created_at'].endswith(('+12:00', '+13:00'))


def test_common_access_log_without_agent_remains_usable(tmp_path):
    log = tmp_path / 'access.log'
    stamp = datetime.now().astimezone().strftime('%d/%b/%Y:%H:%M:%S %z')
    log.write_text(
        '127.0.0.1 - - [{0}] "GET / HTTP/1.1" 200 12\n'.format(stamp),
        encoding='utf-8',
    )
    report = activity_report(
        index_summary={}, packages=(), state_path=tmp_path / 'missing.sqlite',
        access_log=log, days=7,
    )
    assert report['traffic']['page_views'] == 1
    assert report['traffic']['median_ms'] is None


def test_daily_page_views_are_newest_first(tmp_path):
    log = tmp_path / 'access.log'
    now = datetime.now().astimezone()
    stamps = (
        (now - timedelta(days=2)).strftime('%d/%b/%Y:%H:%M:%S %z'),
        now.strftime('%d/%b/%Y:%H:%M:%S %z'),
    )
    log.write_text('\n'.join(
        '127.0.0.1 - - [{0}] "GET / HTTP/1.1" 200 12'.format(stamp)
        for stamp in stamps
    ), encoding='utf-8')

    report = activity_report(
        index_summary={}, packages=(), state_path=tmp_path / 'missing.sqlite',
        access_log=log, days=7,
    )

    dates = [row['date'] for row in report['traffic']['daily']]
    assert dates == sorted(dates, reverse=True)
    assert report['traffic']['daily'][0]['count'] == 1


def test_daily_page_views_use_configured_timezone_at_utc_boundary(tmp_path):
    log = tmp_path / 'access.log'
    log.write_text(
        '127.0.0.1 - - [11/Nov/2026:07:22:53 +1300] '
        '"GET / HTTP/1.1" 200 12\n',
        encoding='utf-8',
    )

    report = _access_projection(
        log, cutoff=None,
        now=datetime(2026, 11, 10, 18, 30, tzinfo=timezone.utc),
        days=7, package_labels={},
        display_timezone=ZoneInfo('Pacific/Auckland'),
    )

    assert report['daily'][0] == {'date': '2026-11-11', 'count': 1}


def test_job_health_is_aggregate_and_non_identifying(tmp_path):
    state = sdf.operator_state.OperatorState(
        tmp_path / 'operator.sqlite'
    ).initialize()
    package = tmp_path / 'package'
    package.mkdir()
    state.queue_package(
        package, target_keys=('simbad:HD 123',), actor='alice',
        reason='Private reason',
    )
    state.queue_acquisition('Vega', actor='bob')
    state.claim_next_acquisition(worker_id='import-1')
    state.queue_sample_acquisition(
        'Private sample', ('HD 123',), actor='alice'
    )

    health = job_health(
        state.path, display_timezone='Pacific/Auckland'
    )

    assert health['queued'] == 2
    assert health['running'] == 1
    assert health['queues']['fits']['queued'] == 1
    assert health['queues']['target imports']['running'] == 1
    assert health['queues']['sample imports']['queued'] == 1
    assert health['timezone'] == 'Pacific/Auckland'
    assert 'alice' not in str(health)
    assert 'Private sample' not in str(health)
