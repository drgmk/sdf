"""Persistent low-volume operator policy and asynchronous fit-job state."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import sqlite3

from .model_policy import TargetModelPolicy


SCHEMA_VERSION = 4


@dataclass(frozen=True)
class PolicyRecord:
    policy: TargetModelPolicy
    main_id: str | None
    scope: str
    actor: str
    reason: str
    updated_at: str


@dataclass(frozen=True)
class FitJob:
    id: int
    package_path: str
    action: str
    payload: dict
    status: str
    actor: str
    reason: str
    created_at: str
    started_at: str | None
    completed_at: str | None
    error: str | None
    worker_id: str | None


@dataclass(frozen=True)
class AcquisitionJob:
    id: int
    query: str
    status: str
    actor: str
    created_at: str
    started_at: str | None
    completed_at: str | None
    result: dict | None
    error: str | None
    worker_id: str | None


class OperatorState:
    """SQLite-backed operator choices that must survive index rebuilds."""

    def __init__(self, path):
        self.path = Path(path).expanduser().resolve()

    def initialize(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with _connection(self.path) as connection:
            connection.executescript('''
                CREATE TABLE IF NOT EXISTS operator_metadata (
                    key TEXT PRIMARY KEY,
                    value TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS fit_policy (
                    target_key TEXT NOT NULL,
                    scope TEXT NOT NULL,
                    main_id TEXT,
                    stellar_models TEXT NOT NULL,
                    disk_models TEXT NOT NULL,
                    max_disks INTEGER NOT NULL,
                    search_mode TEXT NOT NULL,
                    reddening INTEGER NOT NULL,
                    actor TEXT NOT NULL,
                    reason TEXT NOT NULL,
                    updated_at TEXT NOT NULL,
                    PRIMARY KEY (target_key, scope)
                );
                CREATE TABLE IF NOT EXISTS fit_policy_default (
                    target_key TEXT NOT NULL,
                    scope TEXT NOT NULL,
                    actor TEXT NOT NULL,
                    reason TEXT NOT NULL,
                    updated_at TEXT NOT NULL,
                    PRIMARY KEY (target_key, scope)
                );
                CREATE TABLE IF NOT EXISTS fit_jobs (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    package_path TEXT NOT NULL,
                    action TEXT NOT NULL,
                    payload TEXT NOT NULL,
                    status TEXT NOT NULL,
                    actor TEXT NOT NULL,
                    reason TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    started_at TEXT,
                    completed_at TEXT,
                    error TEXT,
                    worker_id TEXT
                );
                CREATE INDEX IF NOT EXISTS fit_jobs_status
                    ON fit_jobs(status, id);
                CREATE INDEX IF NOT EXISTS fit_jobs_active_package
                    ON fit_jobs(status, package_path);
                CREATE TABLE IF NOT EXISTS acquisition_jobs (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    lookup_key TEXT NOT NULL,
                    query TEXT NOT NULL,
                    status TEXT NOT NULL,
                    actor TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    started_at TEXT,
                    completed_at TEXT,
                    result TEXT,
                    error TEXT,
                    worker_id TEXT
                );
                CREATE INDEX IF NOT EXISTS acquisition_jobs_status
                    ON acquisition_jobs(status, id);
                CREATE INDEX IF NOT EXISTS acquisition_jobs_active_query
                    ON acquisition_jobs(status, lookup_key);
            ''')
            row = connection.execute(
                "SELECT value FROM operator_metadata "
                "WHERE key = 'schema_version'"
            ).fetchone()
            if row is None:
                connection.execute(
                    'INSERT INTO operator_metadata VALUES (?, ?)',
                    ('schema_version', str(SCHEMA_VERSION)),
                )
            elif int(row['value']) in {1, 2, 3}:
                if int(row['value']) == 1:
                    columns = {
                        value['name'] for value in connection.execute(
                            'PRAGMA table_info(fit_jobs)'
                        )
                    }
                    if 'worker_id' not in columns:
                        connection.execute(
                            'ALTER TABLE fit_jobs ADD COLUMN worker_id TEXT'
                        )
                connection.execute(
                    "UPDATE operator_metadata SET value = ? "
                    "WHERE key = 'schema_version'",
                    (str(SCHEMA_VERSION),),
                )
            elif int(row['value']) != SCHEMA_VERSION:
                raise ValueError(
                    'unsupported operator-state schema {}'.format(row['value'])
                )
        return self

    def policy(self, target_key, *, package=None):
        """Return the package-specific policy first, then the global policy."""

        scopes = []
        if package is not None:
            scopes.append(_package_scope(package))
        scopes.append('')
        with _connection(self.path) as connection:
            for scope in scopes:
                row = connection.execute(
                    'SELECT * FROM fit_policy '
                    'WHERE target_key = ? AND scope = ?',
                    (target_key, scope),
                ).fetchone()
                if row is not None:
                    return _policy_record(row)
        return None

    def policies(self, target_keys, *, package=None, scoped_only=False):
        records = {}
        for target_key in target_keys:
            if scoped_only:
                record = self._policy_at_scope(
                    target_key, _package_scope(package)
                )
            else:
                record = self.policy(target_key, package=package)
            if record is not None:
                records[target_key] = record
        return records

    def global_policies(self, target_keys):
        return {
            target_key: record
            for target_key in target_keys
            if (record := self._policy_at_scope(target_key, '')) is not None
        }

    def scoped_policies(self, target_keys, package):
        scope = _package_scope(package)
        return {
            target_key: record
            for target_key in target_keys
            if (record := self._policy_at_scope(target_key, scope)) is not None
        }

    def save_policy_and_queue(self, policy, *, package, main_id=None,
                              scope='global', actor, reason):
        """Save one choice and coalesce a pending package refit."""

        return self.save_system_and_queue(
            ((policy, main_id, scope),), package=package,
            actor=actor, reason=reason,
        )

    def save_system_and_queue(self, policies=(), resets=(), *, package,
                              actor, reason, target_keys=()):
        """Apply system policy changes and enqueue one package fit."""

        self.initialize()
        policies = tuple(policies)
        resets = tuple(resets)
        touched = set(target_keys)
        touched.update(policy.target_key for policy, _, _ in policies)
        touched.update(target_key for target_key, _ in resets)
        if not touched:
            raise ValueError('a system fit needs at least one target')
        package_path = str(Path(package).expanduser().resolve())
        now = _now()
        with _connection(self.path) as connection:
            for policy, main_id, scope in policies:
                stored_scope = (
                    '' if scope == 'global' else _package_scope(package)
                )
                connection.execute(
                    'DELETE FROM fit_policy_default '
                    'WHERE target_key = ? AND scope = ?',
                    (policy.target_key, stored_scope),
                )
                _save_policy(
                    connection, policy, stored_scope, main_id,
                    actor, reason, now,
                )
            for target_key, scope in resets:
                stored_scope = (
                    '' if scope == 'global' else _package_scope(package)
                )
                connection.execute(
                    'DELETE FROM fit_policy WHERE target_key = ? AND scope = ?',
                    (target_key, stored_scope),
                )
                connection.execute(
                    '''
                    INSERT INTO fit_policy_default (
                        target_key, scope, actor, reason, updated_at
                    ) VALUES (?, ?, ?, ?, ?)
                    ON CONFLICT(target_key, scope) DO UPDATE SET
                        actor = excluded.actor,
                        reason = excluded.reason,
                        updated_at = excluded.updated_at
                    ''',
                    (target_key, stored_scope, actor, reason, now),
                )
            self._remember_actor(connection, actor)
            return _queue_job(
                connection, package_path, 'policy-refit',
                {'target_keys': sorted(touched)}, actor, reason, now,
            )

    def seed_policy(self, policy, *, package, main_id=None):
        """Recover a missing package policy from its retained default fit."""

        self.initialize()
        scope = _package_scope(package)
        if (self.policy(policy.target_key, package=package) is not None
                or self._uses_default(policy.target_key, scope)
                or self._uses_default(policy.target_key, '')):
            return False
        with _connection(self.path) as connection:
            connection.execute(
                '''
                INSERT OR IGNORE INTO fit_policy (
                    target_key, scope, main_id, stellar_models, disk_models,
                    max_disks, search_mode, reddening, actor, reason, updated_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''',
                (
                    policy.target_key, scope, main_id,
                    json.dumps(policy.stellar_models),
                    json.dumps(policy.disk_models), policy.max_disks,
                    policy.search_mode, int(policy.reddening),
                    'retained fit',
                    'Recovered from the retained default fit', _now(),
                ),
            )
        return True

    def reset_policy_and_queue(self, target_key, *, package, scope='global',
                               actor, reason):
        return self.save_system_and_queue(
            resets=((target_key, scope),), package=package,
            actor=actor, reason=reason,
        )

    def queue_package(self, package, *, target_keys, actor, reason,
                      remember_actor=False):
        """Queue an automatic package fit without changing model policy."""

        self.initialize()
        target_keys = tuple(dict.fromkeys(str(value) for value in target_keys))
        if not target_keys:
            raise ValueError('a package fit needs at least one target')
        package_path = str(Path(package).expanduser().resolve())
        now = _now()
        with _connection(self.path) as connection:
            if remember_actor:
                self._remember_actor(connection, actor)
            return _queue_job(
                connection, package_path, 'automatic-fit',
                {'target_keys': sorted(target_keys)}, actor, reason, now,
            )

    def jobs(self, *, package=None, status=None):
        clauses = []
        values = []
        if package is not None:
            clauses.append('package_path = ?')
            values.append(str(Path(package).expanduser().resolve()))
        if status is not None:
            clauses.append('status = ?')
            values.append(status)
        where = '' if not clauses else ' WHERE ' + ' AND '.join(clauses)
        with _connection(self.path) as connection:
            rows = connection.execute(
                'SELECT * FROM fit_jobs{} ORDER BY id DESC'.format(where),
                values,
            ).fetchall()
        return tuple(_fit_job(row) for row in rows)

    def claim_next(self, *, worker_id=None):
        """Atomically claim work not already running for the same package."""

        self.initialize()
        connection = _connect(self.path)
        try:
            connection.execute('BEGIN IMMEDIATE')
            row = connection.execute(
                "SELECT queued.* FROM fit_jobs AS queued "
                "WHERE queued.status = 'queued' AND NOT EXISTS ("
                "SELECT 1 FROM fit_jobs AS running "
                "WHERE running.status = 'running' "
                "AND running.package_path = queued.package_path) "
                'ORDER BY queued.id LIMIT 1'
            ).fetchone()
            if row is None:
                connection.commit()
                return None
            started_at = _now()
            connection.execute(
                "UPDATE fit_jobs SET status = 'running', started_at = ?, "
                "worker_id = ? "
                'WHERE id = ?',
                (started_at, worker_id, row['id']),
            )
            connection.commit()
            values = dict(row)
            values['status'] = 'running'
            values['started_at'] = started_at
            values['worker_id'] = worker_id
            return _fit_job(values)
        finally:
            connection.close()

    def finish(self, job_id, *, error=None):
        status = 'failed' if error else 'completed'
        with _connection(self.path) as connection:
            cursor = connection.execute(
                '''
                UPDATE fit_jobs
                SET status = ?, completed_at = ?, error = ?
                WHERE id = ? AND status = 'running'
                ''',
                (status, _now(), error, job_id),
            )
            if cursor.rowcount != 1:
                raise ValueError('fit job {} is not running'.format(job_id))

    def recover_running_jobs(self, *, worker_id=None):
        """Return interrupted jobs, optionally for one worker, to the queue."""

        with _connection(self.path) as connection:
            where = "status = 'running'"
            values = []
            if worker_id is not None:
                where += ' AND worker_id = ?'
                values.append(worker_id)
            cursor = connection.execute(
                '''
                UPDATE fit_jobs
                SET status = 'queued', started_at = NULL,
                    completed_at = NULL,
                    error = 'worker stopped before completion',
                    worker_id = NULL
                WHERE {}
                '''.format(where),
                values,
            )
            return cursor.rowcount

    def queue_acquisition(self, query, *, actor):
        """Queue one exact SIMBAD/SDB acquisition, coalescing active requests."""

        self.initialize()
        query = ' '.join(str(query or '').split())
        actor = str(actor or '').strip()
        if not query:
            raise ValueError('a SIMBAD identifier is required')
        if len(query) > 200:
            raise ValueError('SIMBAD identifier is too long')
        if not actor:
            raise ValueError('actor is required')
        lookup_key = query.casefold()
        now = _now()
        with _connection(self.path) as connection:
            existing = connection.execute(
                '''
                SELECT id FROM acquisition_jobs
                WHERE lookup_key = ? AND actor = ?
                    AND status IN ('queued', 'running')
                ORDER BY id LIMIT 1
                ''',
                (lookup_key, actor),
            ).fetchone()
            if existing is not None:
                return existing['id']
            cursor = connection.execute(
                '''
                INSERT INTO acquisition_jobs (
                    lookup_key, query, status, actor, created_at
                ) VALUES (?, ?, 'queued', ?, ?)
                ''',
                (lookup_key, query, actor, now),
            )
            return cursor.lastrowid

    def acquisitions(self, *, actor=None, status=None):
        clauses = []
        values = []
        if actor is not None:
            clauses.append('actor = ?')
            values.append(str(actor))
        if status is not None:
            clauses.append('status = ?')
            values.append(str(status))
        where = '' if not clauses else ' WHERE ' + ' AND '.join(clauses)
        with _connection(self.path) as connection:
            rows = connection.execute(
                'SELECT * FROM acquisition_jobs{} ORDER BY id DESC'.format(
                    where
                ),
                values,
            ).fetchall()
        return tuple(_acquisition_job(row) for row in rows)

    def acquisition(self, job_id):
        with _connection(self.path) as connection:
            row = connection.execute(
                'SELECT * FROM acquisition_jobs WHERE id = ?',
                (int(job_id),),
            ).fetchone()
        if row is None:
            raise KeyError('acquisition job not found: {}'.format(job_id))
        return _acquisition_job(row)

    def claim_next_acquisition(self, *, worker_id=None):
        self.initialize()
        connection = _connect(self.path)
        try:
            connection.execute('BEGIN IMMEDIATE')
            row = connection.execute(
                "SELECT * FROM acquisition_jobs WHERE status = 'queued' "
                'ORDER BY id LIMIT 1'
            ).fetchone()
            if row is None:
                connection.commit()
                return None
            started_at = _now()
            connection.execute(
                "UPDATE acquisition_jobs SET status = 'running', "
                'started_at = ?, worker_id = ?, error = NULL WHERE id = ?',
                (started_at, worker_id, row['id']),
            )
            connection.commit()
            values = dict(row)
            values.update(
                status='running', started_at=started_at,
                worker_id=worker_id, error=None,
            )
            return _acquisition_job(values)
        finally:
            connection.close()

    def finish_acquisition(self, job_id, *, result=None, error=None):
        status = 'failed' if error else 'completed'
        payload = None if result is None else json.dumps(
            result, sort_keys=True
        )
        with _connection(self.path) as connection:
            cursor = connection.execute(
                '''
                UPDATE acquisition_jobs
                SET status = ?, completed_at = ?, result = ?, error = ?
                WHERE id = ? AND status = 'running'
                ''',
                (status, _now(), payload, error, int(job_id)),
            )
            if cursor.rowcount != 1:
                raise ValueError(
                    'acquisition job {} is not running'.format(job_id)
                )

    def recover_running_acquisitions(self):
        with _connection(self.path) as connection:
            cursor = connection.execute(
                '''
                UPDATE acquisition_jobs
                SET status = 'queued', started_at = NULL,
                    completed_at = NULL,
                    error = 'server stopped before completion',
                    worker_id = NULL
                WHERE status = 'running'
                '''
            )
            return cursor.rowcount

    def last_actor(self):
        if not self.path.is_file():
            return None
        with _connection(self.path) as connection:
            row = connection.execute(
                "SELECT value FROM operator_metadata WHERE key = 'last_actor'"
            ).fetchone()
        return None if row is None else row['value']

    def _policy_at_scope(self, target_key, scope):
        if not self.path.is_file():
            return None
        with _connection(self.path) as connection:
            row = connection.execute(
                'SELECT * FROM fit_policy '
                'WHERE target_key = ? AND scope = ?',
                (target_key, scope),
            ).fetchone()
        return None if row is None else _policy_record(row)

    def _uses_default(self, target_key, scope):
        with _connection(self.path) as connection:
            row = connection.execute(
                'SELECT 1 FROM fit_policy_default '
                'WHERE target_key = ? AND scope = ?',
                (target_key, scope),
            ).fetchone()
        return row is not None

    @staticmethod
    def _remember_actor(connection, actor):
        connection.execute(
            '''
            INSERT INTO operator_metadata (key, value) VALUES ('last_actor', ?)
            ON CONFLICT(key) DO UPDATE SET value = excluded.value
            ''',
            (actor,),
        )


def default_operator_state(package, explicit=None):
    if explicit is not None:
        return OperatorState(explicit)
    from .result_index import ResultIndex

    index = ResultIndex.find_for_package(package)
    if index is not None:
        return OperatorState(index.path.with_name('sdf-operator.sqlite'))
    package = Path(package).expanduser().resolve()
    return OperatorState(package.parent / 'sdf-operator.sqlite')


def _save_policy(connection, policy, stored_scope, main_id, actor, reason,
                 now):
    connection.execute(
        '''
        INSERT INTO fit_policy (
            target_key, scope, main_id, stellar_models, disk_models,
            max_disks, search_mode, reddening, actor, reason, updated_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(target_key, scope) DO UPDATE SET
            main_id = excluded.main_id,
            stellar_models = excluded.stellar_models,
            disk_models = excluded.disk_models,
            max_disks = excluded.max_disks,
            search_mode = excluded.search_mode,
            reddening = excluded.reddening,
            actor = excluded.actor,
            reason = excluded.reason,
            updated_at = excluded.updated_at
        ''',
        (
            policy.target_key, stored_scope, main_id,
            json.dumps(policy.stellar_models),
            json.dumps(policy.disk_models), policy.max_disks,
            policy.search_mode, int(policy.reddening), actor, reason, now,
        ),
    )


def _queue_job(connection, package_path, action, payload, actor, reason, now):
    existing = connection.execute(
        '''
        SELECT id, payload FROM fit_jobs
        WHERE package_path = ? AND action = ? AND status = 'queued'
        ORDER BY id LIMIT 1
        ''',
        (package_path, action),
    ).fetchone()
    if existing is not None:
        combined = json.loads(existing['payload'])
        combined['target_keys'] = sorted(set(
            combined.get('target_keys', ()) + payload.get('target_keys', ())
        ))
        connection.execute(
            '''
            UPDATE fit_jobs
            SET payload = ?, actor = ?, reason = ?, created_at = ?
            WHERE id = ?
            ''',
            (json.dumps(combined, sort_keys=True), actor, reason, now,
             existing['id']),
        )
        return existing['id']
    cursor = connection.execute(
        '''
        INSERT INTO fit_jobs (
            package_path, action, payload, status, actor, reason, created_at
        ) VALUES (?, ?, ?, 'queued', ?, ?, ?)
        ''',
        (
            package_path, action, json.dumps(payload, sort_keys=True), actor,
            reason, now,
        ),
    )
    return cursor.lastrowid


def _policy_record(row):
    row = dict(row)
    return PolicyRecord(
        policy=TargetModelPolicy(
            target_key=row['target_key'],
            stellar_models=tuple(json.loads(row['stellar_models'])),
            disk_models=tuple(json.loads(row['disk_models'])),
            max_disks=row['max_disks'], search_mode=row['search_mode'],
            reddening=bool(row['reddening']), origin='operator',
        ),
        main_id=row['main_id'], scope=row['scope'], actor=row['actor'],
        reason=row['reason'], updated_at=row['updated_at'],
    )


def _fit_job(row):
    row = dict(row)
    return FitJob(
        id=row['id'], package_path=row['package_path'], action=row['action'],
        payload=json.loads(row['payload']), status=row['status'],
        actor=row['actor'], reason=row['reason'],
        created_at=row['created_at'], started_at=row['started_at'],
        completed_at=row['completed_at'], error=row['error'],
        worker_id=row['worker_id'],
    )


def _acquisition_job(row):
    row = dict(row)
    return AcquisitionJob(
        id=row['id'], query=row['query'], status=row['status'],
        actor=row['actor'], created_at=row['created_at'],
        started_at=row['started_at'], completed_at=row['completed_at'],
        result=(None if row['result'] is None else json.loads(row['result'])),
        error=row['error'], worker_id=row['worker_id'],
    )


def _package_scope(package):
    if package is None:
        raise ValueError('package is required for package-scoped policy')
    return str(Path(package).expanduser().resolve())


def _now():
    return datetime.now(timezone.utc).isoformat()


def _connect(path):
    connection = sqlite3.connect(path, timeout=30)
    connection.row_factory = sqlite3.Row
    return connection


@contextmanager
def _connection(path):
    connection = _connect(path)
    try:
        with connection:
            yield connection
    finally:
        connection.close()
