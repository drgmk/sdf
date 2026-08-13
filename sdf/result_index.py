"""Rebuildable SQLite index over retained fit candidate artifacts."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import sqlite3
import tempfile


SCHEMA_VERSION = 1


@dataclass(frozen=True)
class IndexSummary:
    packages: int
    targets: int
    candidates: int
    components: int
    parameters: int


class ResultIndex:
    """A disposable query view of package-local candidate libraries."""

    def __init__(self, path):
        self.path = Path(path).expanduser().resolve()

    def rebuild(self, root):
        """Atomically rebuild the index by scanning an export tree."""

        root = Path(root).expanduser().resolve()
        if not root.is_dir():
            raise ValueError('result index root is not a directory: {}'.format(
                root
            ))
        self.path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary_name = tempfile.mkstemp(
            prefix='.sdf-index-', suffix='.sqlite', dir=self.path.parent
        )
        os.close(descriptor)
        temporary = Path(temporary_name)
        try:
            connection = _connect(temporary)
            try:
                _create_schema(connection)
                connection.execute(
                    'INSERT INTO index_metadata (key, value) VALUES (?, ?)',
                    ('root', str(root)),
                )
                for package in _discover_packages(root):
                    _index_package(connection, root, package)
                connection.commit()
                summary = _summary(connection)
            finally:
                connection.close()
            os.replace(temporary, self.path)
        finally:
            if temporary.exists():
                temporary.unlink()
        return summary

    def summary(self):
        with _connection(self.path) as connection:
            return _summary(connection)

    @classmethod
    def find_for_package(cls, package):
        """Find the nearest existing default index containing ``package``."""

        package = Path(package).expanduser().resolve()
        for parent in (package, *package.parents):
            path = parent / 'sdf-results.sqlite'
            if path.is_file():
                index = cls(path)
                try:
                    root = index.root()
                except (sqlite3.DatabaseError, TypeError, ValueError):
                    continue
                if package == root or root in package.parents:
                    return index
        return None

    def root(self):
        """Return the export root recorded by the index."""

        with _connection(self.path) as connection:
            row = connection.execute(
                "SELECT value FROM index_metadata WHERE key = 'root'"
            ).fetchone()
        if row is None:
            raise ValueError('{} is not an SDF result index'.format(self.path))
        return Path(row[0]).resolve()

    def update_package(self, package):
        """Transactionally refresh one package after a fit or selection."""

        package = Path(package).expanduser().resolve()
        root = self.root()
        if package != root and root not in package.parents:
            raise ValueError('{} is outside index root {}'.format(package, root))
        package_path = package.relative_to(root).as_posix()
        with _connection(self.path) as connection:
            connection.execute(
                'DELETE FROM fit_packages WHERE path = ?', (package_path,)
            )
            if tuple((package / 'fits').glob('*/result.json')):
                _index_package(connection, root, package)
        return self.summary()

    def candidates_for_target(self, target_key):
        """Return indexed candidate rows for a canonical target identity."""

        with _connection(self.path) as connection:
            rows = connection.execute(
                '''
                SELECT DISTINCT c.*, s.position AS display_position,
                       CASE WHEN p.default_candidate_id = c.candidate_id
                            THEN 1 ELSE 0 END AS is_default
                FROM fit_targets AS t
                JOIN fit_candidates AS c ON c.package_path = t.package_path
                JOIN fit_packages AS p ON p.path = c.package_path
                LEFT JOIN fit_selection AS s
                  ON s.package_path = c.package_path
                 AND s.candidate_id = c.candidate_id
                WHERE t.target_key = ?
                ORDER BY is_default DESC, display_position IS NULL,
                         display_position, c.evidence DESC
                ''',
                (target_key,),
            ).fetchall()
        return tuple(dict(row) for row in rows)

    def default_targets(self, identifiers=None):
        """Return physical targets with their package's default candidate."""

        with _connection(self.path) as connection:
            where = _requested_identifiers(
                connection, identifiers, prefix='WHERE'
            )
            if where is None:
                return ()
            rows = connection.execute(
                '''
                SELECT t.package_path, t.target_key, t.main_id, t.sdbid,
                       c.candidate_id, c.result_path, c.display_name,
                       c.evidence, c.evidence_error, c.data_current,
                       c.completed_at
                FROM fit_targets AS t
                JOIN fit_packages AS p ON p.path = t.package_path
                JOIN fit_candidates AS c
                  ON c.package_path = p.path
                 AND c.candidate_id = p.default_candidate_id
                {where}
                ORDER BY t.target_key, t.package_path
                '''.format(where=where),
            ).fetchall()
        return tuple(dict(row) for row in rows)

    def default_parameters(self, identifiers=None):
        """Return indexed scalar results from default candidates."""

        with _connection(self.path) as connection:
            where = _requested_identifiers(
                connection, identifiers, prefix='AND'
            )
            if where is None:
                return ()
            rows = connection.execute(
                '''
                SELECT t.package_path, t.target_key, t.main_id, t.sdbid,
                       fp.candidate_id, fp.component_id, fp.role,
                       fp.quantity_scope, fp.name, fp.median, fp.minus,
                       fp.plus, fp.unit
                FROM fit_targets AS t
                JOIN fit_packages AS p ON p.path = t.package_path
                JOIN fit_parameters AS fp
                  ON fp.package_path = p.path
                 AND fp.candidate_id = p.default_candidate_id
                 AND fp.target_key = t.target_key
                WHERE 1 = 1 {where}
                ORDER BY t.target_key, fp.role, fp.component_id,
                         fp.quantity_scope, fp.name
                '''.format(where=where),
            ).fetchall()
        return tuple(dict(row) for row in rows)

    def verify(self, root=None):
        """Return human-readable consistency problems, or an empty tuple."""

        if not self.path.is_file():
            return ('result index does not exist: {}'.format(self.path),)
        problems = []
        with _connection(self.path) as connection:
            integrity = connection.execute('PRAGMA integrity_check').fetchone()[0]
            if integrity != 'ok':
                problems.append('SQLite integrity check: {}'.format(integrity))
            indexed_root = Path(connection.execute(
                "SELECT value FROM index_metadata WHERE key = 'root'"
            ).fetchone()[0])
            root = indexed_root if root is None else Path(root).expanduser().resolve()
            indexed = {
                row['result_path']
                for row in connection.execute(
                    'SELECT result_path FROM fit_candidates'
                )
            }
            dangling = sorted(
                path for path in indexed if not (root / path).is_file()
            )
            problems.extend(
                'missing candidate artifact: {}'.format(path)
                for path in dangling
            )
            discovered = {
                path.relative_to(root).as_posix()
                for path in _discover_results(root)
            }
            for path in sorted(discovered - indexed):
                problems.append('candidate absent from index: {}'.format(path))
            for path in sorted(indexed - discovered):
                if path not in dangling:
                    problems.append('indexed candidate absent from scan: {}'.format(
                        path
                    ))
        return tuple(problems)


def _connect(path):
    connection = sqlite3.connect(path, timeout=30)
    connection.row_factory = sqlite3.Row
    connection.execute('PRAGMA foreign_keys = ON')
    return connection


@contextmanager
def _connection(path):
    connection = _connect(path)
    try:
        with connection:
            yield connection
    finally:
        connection.close()


def _create_schema(connection):
    connection.executescript('''
        CREATE TABLE index_metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        INSERT INTO index_metadata (key, value)
            VALUES ('schema_version', '{schema_version}');

        CREATE TABLE fit_packages (
            path TEXT PRIMARY KEY,
            input_fingerprint TEXT,
            default_candidate_id TEXT,
            updated_at TEXT NOT NULL
        );

        CREATE TABLE fit_targets (
            package_path TEXT NOT NULL REFERENCES fit_packages(path)
                ON DELETE CASCADE,
            target_key TEXT NOT NULL,
            main_id TEXT,
            sdbid TEXT,
            PRIMARY KEY (package_path, target_key)
        );
        CREATE INDEX fit_targets_identity ON fit_targets(target_key);
        CREATE INDEX fit_targets_sdbid ON fit_targets(sdbid);

        CREATE TABLE fit_candidates (
            package_path TEXT NOT NULL REFERENCES fit_packages(path)
                ON DELETE CASCADE,
            candidate_id TEXT NOT NULL,
            result_path TEXT NOT NULL UNIQUE,
            display_name TEXT NOT NULL,
            input_fingerprint TEXT NOT NULL,
            fit_fingerprint TEXT NOT NULL,
            data_current INTEGER NOT NULL,
            completed_at TEXT NOT NULL,
            evidence REAL,
            evidence_error REAL,
            sampler TEXT NOT NULL,
            PRIMARY KEY (package_path, candidate_id)
        );
        CREATE INDEX fit_candidates_evidence
            ON fit_candidates(package_path, evidence DESC);

        CREATE TABLE fit_components (
            package_path TEXT NOT NULL,
            candidate_id TEXT NOT NULL,
            component_id TEXT NOT NULL,
            owner_key TEXT NOT NULL,
            role TEXT NOT NULL,
            model_name TEXT NOT NULL,
            ordinal INTEGER NOT NULL,
            PRIMARY KEY (package_path, candidate_id, component_id),
            FOREIGN KEY (package_path, candidate_id)
                REFERENCES fit_candidates(package_path, candidate_id)
                ON DELETE CASCADE
        );
        CREATE INDEX fit_components_owner ON fit_components(owner_key);

        CREATE TABLE fit_parameters (
            package_path TEXT NOT NULL,
            candidate_id TEXT NOT NULL,
            target_key TEXT NOT NULL,
            component_id TEXT NOT NULL,
            role TEXT NOT NULL,
            quantity_scope TEXT NOT NULL,
            name TEXT NOT NULL,
            median REAL NOT NULL,
            minus REAL,
            plus REAL,
            unit TEXT,
            PRIMARY KEY (
                package_path, candidate_id, target_key, component_id,
                quantity_scope, name
            ),
            FOREIGN KEY (package_path, candidate_id)
                REFERENCES fit_candidates(package_path, candidate_id)
                ON DELETE CASCADE
        );
        CREATE INDEX fit_parameters_target
            ON fit_parameters(target_key, name);

        CREATE TABLE fit_selection (
            package_path TEXT NOT NULL,
            position INTEGER NOT NULL,
            candidate_id TEXT NOT NULL,
            label TEXT,
            PRIMARY KEY (package_path, position),
            UNIQUE (package_path, candidate_id),
            FOREIGN KEY (package_path, candidate_id)
                REFERENCES fit_candidates(package_path, candidate_id)
                ON DELETE CASCADE
        );
    '''.format(schema_version=SCHEMA_VERSION))


def _requested_identifiers(connection, identifiers, prefix):
    if identifiers is None:
        return ''
    identifiers = tuple(dict.fromkeys(str(value) for value in identifiers))
    if not identifiers:
        return None
    connection.execute(
        'CREATE TEMP TABLE requested_identifiers (value TEXT PRIMARY KEY)'
    )
    connection.executemany(
        'INSERT INTO requested_identifiers VALUES (?)',
        ((value,) for value in identifiers),
    )
    return '''
        {prefix} EXISTS (
            SELECT 1 FROM requested_identifiers AS requested
            WHERE requested.value IN (t.target_key, t.main_id, t.sdbid)
        )
    '''.format(prefix=prefix)


def _discover_results(root):
    return tuple(sorted(
        path for path in root.rglob('result.json')
        if path.parent.parent.name == 'fits'
    ))


def _discover_packages(root):
    return tuple(sorted({path.parents[2] for path in _discover_results(root)}))


def _index_package(connection, root, package):
    package_path = package.relative_to(root).as_posix()
    result_paths = tuple(sorted((package / 'fits').glob('*/result.json')))
    results = [json.loads(path.read_text(encoding='utf-8'))
               for path in result_paths]
    current_input = _package_input_fingerprint(package)
    selection = _read_selection(package)
    default = selection.get('default')
    available = {result['candidate_id'] for result in results}
    if default is not None and default not in available:
        raise ValueError('{} selects unknown default candidate {}'.format(
            package, default
        ))
    displayed = selection.get('displayed', ())
    missing = set(displayed) - available
    if missing:
        raise ValueError('{} selects unknown candidates: {}'.format(
            package, ', '.join(sorted(missing))
        ))

    connection.execute(
        '''
        INSERT INTO fit_packages
            (path, input_fingerprint, default_candidate_id, updated_at)
        VALUES (?, ?, ?, ?)
        ''',
        (
            package_path, current_input, default,
            datetime.now(timezone.utc).isoformat(),
        ),
    )
    targets = {}
    for path, result in zip(result_paths, results):
        posterior = result['posterior']
        connection.execute(
            '''
            INSERT INTO fit_candidates (
                package_path, candidate_id, result_path, display_name,
                input_fingerprint, fit_fingerprint, data_current,
                completed_at, evidence, evidence_error, sampler
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''',
            (
                package_path, result['candidate_id'],
                path.relative_to(root).as_posix(), result['display_name'],
                result['input_fingerprint'], result['fit_fingerprint'],
                int(current_input == result['input_fingerprint']),
                result['completed_at'], posterior.get('evidence'),
                posterior.get('evidence_error'), posterior['sampler'],
            ),
        )
        for component in result.get('components', ()):
            connection.execute(
                '''
                INSERT INTO fit_components VALUES (?, ?, ?, ?, ?, ?, ?)
                ''',
                (
                    package_path, result['candidate_id'],
                    component['component_id'], component['owner_key'],
                    component['role'], component['model_name'],
                    component['ordinal'],
                ),
            )
            targets.setdefault(component['owner_key'], {})
        for target_key, target in result.get('analysis', {}).get(
                'targets', {}).items():
            targets[target_key] = target
        _index_parameters(
            connection, package_path, result['candidate_id'],
            result.get('analysis', {}),
        )

    for target_key, target in sorted(targets.items()):
        connection.execute(
            'INSERT INTO fit_targets VALUES (?, ?, ?, ?)',
            (
                package_path, target_key, target.get('main_id'),
                target.get('sdbid'),
            ),
        )
    labels = selection.get('labels', {})
    for position, candidate_id in enumerate(displayed):
        connection.execute(
            'INSERT INTO fit_selection VALUES (?, ?, ?, ?)',
            (package_path, position, candidate_id, labels.get(candidate_id)),
        )


def _read_selection(package):
    path = package / 'fits' / 'selection.json'
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding='utf-8'))


def _index_parameters(connection, package_path, candidate_id, analysis):
    for component in analysis.get('components', ()):
        for name, summary in component.get('parameters', {}).items():
            _insert_parameter(
                connection, package_path, candidate_id,
                component['owner_key'], component['component_id'],
                component['role'], 'fitted', name, summary,
            )
    for target_key, target in analysis.get('targets', {}).items():
        for name, summary in target.items():
            if _is_summary(summary):
                _insert_parameter(
                    connection, package_path, candidate_id, target_key,
                    '', 'target', 'derived', name, summary,
                )
        for role, field in (
            ('stellar', 'stellar_components'), ('disk', 'disk_components')
        ):
            for component in target.get(field, ()):
                for name, summary in component.items():
                    if _is_summary(summary):
                        _insert_parameter(
                            connection, package_path, candidate_id, target_key,
                            component['component_id'], role, 'derived', name,
                            summary,
                        )


def _insert_parameter(connection, package_path, candidate_id, target_key,
                      component_id, role, quantity_scope, name, summary):
    connection.execute(
        '''
        INSERT INTO fit_parameters VALUES (
            ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
        )
        ''',
        (
            package_path, candidate_id, target_key, component_id, role,
            quantity_scope, name, summary['median'], summary.get('minus'),
            summary.get('plus'), summary.get('unit'),
        ),
    )


def _is_summary(value):
    return isinstance(value, dict) and 'median' in value


def _package_input_fingerprint(package):
    paths = set(package.glob('*-rawphot.txt'))
    topology = package / 'joint-fit.yml'
    if topology.is_file():
        paths.add(topology)
    if not paths:
        return None
    digest = sha256()
    for path in sorted(paths, key=lambda value: value.name):
        digest.update(path.name.encode('utf-8'))
        digest.update(b'\0')
        digest.update(path.read_bytes())
        digest.update(b'\0')
    return digest.hexdigest()


def _summary(connection):
    counts = {}
    for label, table in (
        ('packages', 'fit_packages'), ('targets', 'fit_targets'),
        ('candidates', 'fit_candidates'), ('components', 'fit_components'),
        ('parameters', 'fit_parameters'),
    ):
        counts[label] = connection.execute(
            'SELECT count(*) FROM {}'.format(table)
        ).fetchone()[0]
    return IndexSummary(**counts)
