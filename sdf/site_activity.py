"""Read-only website inventory, activity, and access-log summaries."""

from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
from pathlib import Path
import re
import sqlite3
from statistics import median
from urllib.parse import parse_qs, unquote, urlsplit


_ACCESS_LINE = re.compile(
    r'^(?P<host>\S+) \S+ \S+ \[(?P<when>[^]]+)\] '
    r'"(?P<method>\S+) (?P<target>\S+)(?: [^"]+)?" '
    r'(?P<status>\d{3}) (?P<bytes>\S+)'
    r'(?: "(?P<referer>[^"]*)" "(?P<agent>[^"]*)" '
    r'(?P<duration>\d+))?$'
)
_BOT_AGENT = re.compile(
    r'bot|crawler|spider|slurp|headless|python-requests|curl/|wget/', re.I
)
_DOWNLOAD_SUFFIXES = {
    '.csv', '.fits', '.json', '.npz', '.pdf', '.pickle', '.pkl', '.png',
    '.svg', '.txt', '.xml', '.yaml', '.yml',
}


def activity_report(*, index_summary, packages, state_path, sdb_path=None,
                    access_log=None, days=30, limit=60):
    """Build the complete admin activity projection.

    This deliberately reads existing authoritative/audit databases.  It does
    not copy domain actions into a second event system.  Access-log parsing
    discards client addresses and retains only aggregate route information.
    """

    days = max(0, int(days))
    now = datetime.now(timezone.utc)
    cutoff = None if days == 0 else now - timedelta(days=days)
    packages = tuple(packages)
    package_labels = {
        row['package_path']: (
            row.get('target_labels') or row.get('sdbids')
            or (row['package_path'],)
        )[0]
        for row in packages
    }
    status = Counter(row['status'] for row in packages)

    operator = _operator_projection(
        state_path, cutoff=cutoff, package_labels=package_labels
    )
    sdb = _sdb_projection(sdb_path, cutoff=cutoff)
    traffic = _access_projection(
        access_log, cutoff=cutoff, now=now, days=days,
        package_labels=package_labels,
    )
    events = sorted(
        (*operator['events'], *sdb['events']),
        key=lambda row: _timestamp(row['created_at']) or datetime.min.replace(
            tzinfo=timezone.utc
        ), reverse=True,
    )[:limit]

    inventory = dict(index_summary)
    inventory.update(sdb['inventory'])
    inventory['current_packages'] = status['current']
    inventory['unfit_packages'] = status['unfit']
    inventory['stale_packages'] = status['stale']
    inventory['selection_packages'] = status['needs selection']

    changes = dict(operator['counts'])
    changes.update(sdb['counts'])
    return {
        'days': days,
        'period': 'all recorded time' if days == 0 else 'last {} days'.format(days),
        'inventory': inventory,
        'changes': changes,
        'jobs': operator['jobs'],
        'traffic': traffic,
        'events': events,
    }


def _operator_projection(path, *, cutoff, package_labels):
    empty = {
        'counts': {
            'fits_completed': 0, 'target_imports': 0,
            'sample_imports': 0,
        },
        'jobs': {}, 'events': (),
    }
    path = Path(path)
    if not path.is_file():
        return empty
    with _read_connection(path) as connection:
        tables = _tables(connection)
        events = []
        counts = dict(empty['counts'])
        jobs = {}
        if 'fit_jobs' in tables:
            rows = connection.execute(
                'SELECT * FROM fit_jobs ORDER BY created_at DESC'
            ).fetchall()
            relevant = [row for row in rows if _after(_job_time(row), cutoff)]
            jobs['fits'] = dict(Counter(row['status'] for row in rows))
            counts['fits_completed'] = sum(
                row['status'] == 'completed'
                and _after(row['completed_at'], cutoff)
                for row in rows
            )
            for row in relevant:
                label = package_labels.get(
                    _relative_package(row['package_path'], package_labels),
                    Path(row['package_path']).name,
                )
                events.append(_event(
                    _job_time(row), 'Fit',
                    '{} fit for {}'.format(_job_verb(row['status']), label),
                    row['actor'],
                ))
        if 'acquisition_jobs' in tables:
            rows = connection.execute(
                'SELECT * FROM acquisition_jobs ORDER BY created_at DESC'
            ).fetchall()
            relevant = [row for row in rows if _after(_job_time(row), cutoff)]
            jobs['imports'] = dict(Counter(row['status'] for row in rows))
            counts['target_imports'] = sum(
                row['status'] == 'completed'
                and _after(row['completed_at'], cutoff)
                for row in rows
            )
            events.extend(
                _event(
                    _job_time(row), 'Import',
                    '{} target import for {}'.format(
                        _job_verb(row['status']), row['query']
                    ), row['actor'],
                )
                for row in relevant
            )
        if 'sample_acquisition_jobs' in tables:
            rows = connection.execute(
                'SELECT * FROM sample_acquisition_jobs ORDER BY created_at DESC'
            ).fetchall()
            relevant = [row for row in rows if _after(_job_time(row), cutoff)]
            jobs['samples'] = dict(Counter(row['status'] for row in rows))
            counts['sample_imports'] = sum(
                row['status'] == 'completed'
                and _after(row['completed_at'], cutoff)
                for row in rows
            )
            events.extend(
                _event(
                    _job_time(row), 'Sample',
                    '{} sample {} for {}'.format(
                        _job_verb(row['status']),
                        row['operation'] if 'operation' in row.keys() else 'add',
                        row['sample_name'],
                    ), row['actor'],
                )
                for row in relevant
            )
    return {'counts': counts, 'jobs': jobs, 'events': tuple(events)}


def _sdb_projection(path, *, cutoff):
    empty = {
        'inventory': {
            'database_targets': 0, 'samples': 0, 'public_samples': 0,
            'private_samples': 0, 'measurements': 0, 'spectra': 0,
            'users': 0,
        },
        'counts': {
            'targets_added': 0, 'sample_edits': 0,
            'review_actions': 0, 'targets_exported': 0,
        },
        'events': (),
    }
    if path is None or not Path(path).is_file():
        return empty
    with _read_connection(path) as connection:
        tables = _tables(connection)
        inventory = dict(empty['inventory'])
        inventory['database_targets'] = _count(connection, tables, 'targets')
        inventory['samples'] = _count(connection, tables, 'samples')
        inventory['measurements'] = _count(
            connection, tables, 'normalized_measurements'
        )
        inventory['spectra'] = _count(connection, tables, 'spectral_assets')
        inventory['users'] = _count(connection, tables, 'users')
        if 'samples' in tables:
            visibility = dict(connection.execute(
                'SELECT visibility, count(*) FROM samples GROUP BY visibility'
            ).fetchall())
            inventory['public_samples'] = visibility.get('public', 0)
            inventory['private_samples'] = visibility.get('private', 0)

        counts = dict(empty['counts'])
        counts['targets_added'] = _period_count(
            connection, tables, 'targets', 'created_at', cutoff
        )
        counts['sample_edits'] = _period_count(
            connection, tables, 'sample_membership_actions', 'created_at', cutoff
        )
        review_tables = (
            'measurement_association_actions',
            'measurement_eligibility_actions',
            'catalog_target_association_actions',
            'spectral_product_association_actions',
            'match_decisions', 'catalog_result_decision_history',
        )
        counts['review_actions'] = sum(
            _period_count(connection, tables, table, 'created_at', cutoff)
            for table in review_tables
        )
        if 'export_dirty_targets' in tables:
            sql = (
                'SELECT count(DISTINCT target_id) FROM export_dirty_targets '
                'WHERE exported_at IS NOT NULL'
            )
            parameters = ()
            if cutoff is not None:
                sql += ' AND exported_at >= ?'
                parameters = (_sqlite_time(cutoff),)
            counts['targets_exported'] = connection.execute(
                sql, parameters
            ).fetchone()[0]
        events = _sdb_events(connection, tables, cutoff)
    return {'inventory': inventory, 'counts': counts, 'events': events}


def _sdb_events(connection, tables, cutoff):
    events = []
    if {'sample_membership_actions', 'samples', 'targets'} <= tables:
        rows = connection.execute('''
            SELECT a.created_at, a.actor, a.action, s.name, t.sdbid
            FROM sample_membership_actions AS a
            JOIN samples AS s ON s.id = a.sample_id
            JOIN targets AS t ON t.id = a.target_id
            ORDER BY a.created_at DESC LIMIT 250
        ''').fetchall()
        for row in rows:
            if _after(row['created_at'], cutoff):
                verb = 'Added' if row['action'] == 'add' else 'Removed'
                events.append(_event(
                    row['created_at'], 'Sample',
                    '{} {} {} sample {}'.format(
                        verb, row['sdbid'],
                        'to' if row['action'] == 'add' else 'from', row['name'],
                    ), row['actor'],
                ))
    if 'sample_deletions' in tables:
        rows = connection.execute(
            'SELECT created_at, actor, name FROM sample_deletions '
            'ORDER BY created_at DESC LIMIT 100'
        ).fetchall()
        events.extend(
            _event(
                row['created_at'], 'Sample',
                'Deleted sample {}'.format(row['name']), row['actor'],
            )
            for row in rows if _after(row['created_at'], cutoff)
        )
    action_tables = (
        ('measurement_association_actions', 'Photometry'),
        ('measurement_eligibility_actions', 'Photometry'),
        ('catalog_target_association_actions', 'Catalog source'),
        ('spectral_product_association_actions', 'Spectrum'),
    )
    for table, label in action_tables:
        if table not in tables:
            continue
        columns = _columns(connection, table)
        if not {'created_at', 'actor', 'action'} <= columns:
            continue
        rows = connection.execute(
            'SELECT created_at, actor, action FROM {} '
            'ORDER BY created_at DESC LIMIT 100'.format(table)
        ).fetchall()
        events.extend(
            _event(
                row['created_at'], 'Review',
                '{} {}'.format(label, row['action']), row['actor'],
            )
            for row in rows if _after(row['created_at'], cutoff)
        )
    return tuple(events)


def _access_projection(path, *, cutoff, now, days, package_labels):
    empty = {
        'available': False, 'page_views': 0, 'downloads': 0,
        'server_errors': 0, 'median_ms': None, 'daily': (),
        'daily_max': 0,
        'top_systems': (), 'top_samples': (), 'top_searches': (),
    }
    if path is None or not Path(path).is_file():
        return empty
    daily = Counter()
    systems = Counter()
    samples = Counter()
    searches = Counter()
    durations = []
    page_views = downloads = server_errors = 0
    with Path(path).open(encoding='utf-8', errors='replace') as stream:
        for line in stream:
            record = _parse_access_line(line)
            if record is None or not _after(record['when'], cutoff):
                continue
            if record['status'] >= 500:
                server_errors += 1
            if record['duration_ms'] is not None:
                durations.append(record['duration_ms'])
            if record['method'] != 'GET' or not 200 <= record['status'] < 300:
                continue
            if record['agent'] and _BOT_AGENT.search(record['agent']):
                continue
            split = urlsplit(record['target'])
            path_value = unquote(split.path)
            if _is_download(path_value):
                downloads += 1
                continue
            kind, identifier = _page_identity(path_value)
            if kind is None:
                continue
            page_views += 1
            daily[record['when'].date().isoformat()] += 1
            if kind == 'system' and identifier:
                systems[package_labels.get(identifier, identifier)] += 1
            elif kind == 'sample' and identifier:
                samples[identifier] += 1
            elif kind == 'search':
                query = parse_qs(split.query).get('q', [''])[0].strip()
                if query:
                    searches[query] += 1
    if days and days <= 90:
        start = (now - timedelta(days=days - 1)).date()
        daily_rows = tuple({
            'date': (start + timedelta(days=offset)).isoformat(),
            'count': daily[(start + timedelta(days=offset)).isoformat()],
        } for offset in range(days))
    else:
        daily_rows = tuple(
            {'date': date, 'count': count}
            for date, count in sorted(daily.items())
        )
    return {
        'available': True, 'page_views': page_views,
        'downloads': downloads, 'server_errors': server_errors,
        'median_ms': None if not durations else round(median(durations), 1),
        'daily': daily_rows,
        'daily_max': max((row['count'] for row in daily_rows), default=0),
        'top_systems': _top(systems), 'top_samples': _top(samples),
        'top_searches': _top(searches),
    }


def _parse_access_line(line):
    match = _ACCESS_LINE.match(line.strip())
    if match is None:
        return None
    try:
        when = datetime.strptime(
            match.group('when'), '%d/%b/%Y:%H:%M:%S %z'
        ).astimezone(timezone.utc)
        duration = match.group('duration')
        return {
            'when': when, 'method': match.group('method'),
            'target': match.group('target'),
            'status': int(match.group('status')),
            'agent': match.group('agent'),
            'duration_ms': None if duration is None else int(duration) / 1000,
        }
    except (TypeError, ValueError):
        return None


def _page_identity(path):
    clean = path.rstrip('/') or '/'
    if clean == '/':
        return 'home', None
    if clean in {'/search', '/samples', '/admin', '/admin/search',
                 '/admin/jobs', '/admin/activity', '/imports'}:
        return clean.strip('/').replace('/', '-'), None
    if clean.startswith('/system/'):
        return 'system', clean.removeprefix('/system/').split('/', 1)[0]
    if clean.startswith('/admin/system/'):
        return 'system', clean.removeprefix('/admin/system/').split('/', 1)[0]
    if clean.startswith('/samples/'):
        return 'sample', clean.removeprefix('/samples/').split('/', 1)[0]
    if clean.startswith('/imports/') or clean.startswith('/sample-imports/'):
        return 'import', None
    if clean == '/review' or clean.startswith('/review/'):
        return 'review', None
    if clean == '/docs' or clean.startswith('/docs/'):
        if Path(clean).suffix not in {'', '.html'}:
            return None, None
        return 'docs', None
    if clean == '/herschel' or clean.startswith('/herschel/'):
        if Path(clean).suffix not in {'', '.html'}:
            return None, None
        return 'herschel', None
    return None, None


def _is_download(path):
    lower = path.lower()
    if lower.endswith(('/sample.xml', '/photometry.csv')):
        return True
    if '/system/' not in lower and '/admin/system/' not in lower:
        return False
    return Path(lower).suffix in _DOWNLOAD_SUFFIXES


def _top(counter, limit=10):
    return tuple(
        {'label': label, 'count': count}
        for label, count in counter.most_common(limit)
    )


def _event(created_at, kind, description, actor):
    return {
        'created_at': str(created_at), 'kind': kind,
        'description': description, 'actor': str(actor or 'system'),
    }


def _job_verb(status):
    return {
        'completed': 'Completed', 'failed': 'Failed', 'running': 'Running',
        'queued': 'Queued',
    }.get(status, str(status).title())


def _job_time(row):
    if row['status'] in {'completed', 'failed'} and row['completed_at']:
        return row['completed_at']
    if row['status'] == 'running' and row['started_at']:
        return row['started_at']
    return row['created_at']


def _relative_package(package, labels):
    package = Path(package)
    for relative in labels:
        if package.as_posix().endswith('/' + relative) or package.name == relative:
            return relative
    return package.name


def _after(value, cutoff):
    if cutoff is None:
        return True
    parsed = _timestamp(value)
    return parsed is not None and parsed >= cutoff


def _timestamp(value):
    if isinstance(value, datetime):
        parsed = value
    else:
        text = str(value or '').strip().replace(' ', 'T')
        if not text:
            return None
        try:
            parsed = datetime.fromisoformat(text.replace('Z', '+00:00'))
        except ValueError:
            return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _sqlite_time(value):
    return value.astimezone(timezone.utc).replace(tzinfo=None).isoformat(' ')


def _count(connection, tables, table):
    if table not in tables:
        return 0
    return connection.execute(
        'SELECT count(*) FROM {}'.format(table)
    ).fetchone()[0]


def _period_count(connection, tables, table, column, cutoff):
    if table not in tables or column not in _columns(connection, table):
        return 0
    sql = 'SELECT count(*) FROM {}'.format(table)
    parameters = ()
    if cutoff is not None:
        sql += ' WHERE {} >= ?'.format(column)
        parameters = (_sqlite_time(cutoff),)
    return connection.execute(sql, parameters).fetchone()[0]


def _tables(connection):
    return {
        row[0] for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'"
        )
    }


def _columns(connection, table):
    return {
        row[1] for row in connection.execute(
            'PRAGMA table_info({})'.format(table)
        )
    }


@contextmanager
def _read_connection(path):
    connection = sqlite3.connect(
        'file:{}?mode=ro'.format(Path(path).resolve()), uri=True, timeout=5.0
    )
    connection.row_factory = sqlite3.Row
    try:
        yield connection
    finally:
        connection.close()
