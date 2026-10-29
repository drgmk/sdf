"""Trusted-proxy identities and SDB-backed catalog authorization."""

from __future__ import annotations

import base64
import binascii
from dataclasses import dataclass
from hashlib import sha256
import hmac
import json
from pathlib import Path
import secrets
import sqlite3
import time


class AuthenticationRequired(PermissionError):
    """The requested catalog view requires an authenticated identity."""


class AuthorizationDenied(PermissionError):
    """The authenticated identity cannot access the requested resource."""


@dataclass(frozen=True)
class CatalogPrincipal:
    """One active application user authenticated by the reverse proxy."""

    username: str | None
    is_admin: bool = False

    @property
    def authenticated(self):
        return self.username is not None


@dataclass(frozen=True)
class SamplePolicy:
    name: str
    owner_username: str | None
    visibility: str
    granted: bool


class AccessDatabase:
    """Small read-only authorization projection over the authoritative SDB."""

    def __init__(self, path):
        self.path = Path(path).expanduser().resolve()
        if not self.path.is_file():
            raise ValueError(
                'authorization database does not exist: {}'.format(self.path)
            )

    def user(self, username):
        with self._connect() as connection:
            row = connection.execute(
                'SELECT username, is_admin, disabled FROM users '
                'WHERE username = ?',
                (str(username),),
            ).fetchone()
        if row is None:
            return None
        return CatalogPrincipal(row[0], bool(row[1])) if not row[2] else False

    def sample(self, sample_name, username=None):
        with self._connect() as connection:
            row = connection.execute(
                '''
                SELECT s.name, s.owner_username, s.visibility,
                       CASE WHEN g.username IS NULL THEN 0 ELSE 1 END
                FROM samples AS s
                LEFT JOIN sample_grants AS g
                  ON g.sample_id = s.id AND g.username = ?
                WHERE s.name = ?
                ''',
                (username, str(sample_name)),
            ).fetchone()
        if row is None:
            return None
        return SamplePolicy(row[0], row[1], row[2], bool(row[3]))

    def _connect(self):
        connection = sqlite3.connect(self.path, timeout=5.0)
        connection.execute('PRAGMA query_only = ON')
        return connection


class CatalogAuthorizer:
    """Authorize catalog routes using proxy identity and SDB policy."""

    def __init__(self, *, user_header=None, access_database=None,
                 local_actor=None, session_secret=None,
                 session_seconds=12 * 60 * 60):
        self.user_header = str(user_header or '').strip() or None
        self.local_actor = str(local_actor or '').strip() or 'local'
        if self.user_header and access_database is None:
            raise ValueError(
                'an authorization database is required with a proxy user header'
            )
        if not self.user_header and access_database is not None:
            raise ValueError(
                'a proxy user header is required with an authorization database'
            )
        self.access = (
            AccessDatabase(access_database) if self.user_header else None
        )
        self.session_secret = session_secret or secrets.token_bytes(32)
        if isinstance(self.session_secret, str):
            self.session_secret = self.session_secret.encode('utf-8')
        self.session_seconds = int(session_seconds)
        if self.session_seconds < 1:
            raise ValueError('session lifetime must be positive')

    @property
    def proxy_enabled(self):
        return self.user_header is not None

    def principal(self, headers):
        if not self.proxy_enabled:
            return CatalogPrincipal(self.local_actor, True)
        username = str(headers.get(self.user_header, '')).strip()
        if not username:
            username = self._session_username(headers.get('Cookie')) or ''
        if not username:
            return CatalogPrincipal(None)
        if any(character in username for character in '\r\n\0'):
            raise AuthenticationRequired('invalid authenticated user header')
        principal = self.access.user(username)
        if principal is None:
            raise AuthorizationDenied('application user is not registered')
        if principal is False:
            raise AuthorizationDenied('application user is disabled')
        return principal

    def issue_session(self, username, *, now=None):
        """Return a signed, expiring token for an authenticated active user."""

        username = str(username or '').strip()
        if not username or any(character in username for character in '\r\n\0'):
            raise AuthenticationRequired('invalid authenticated username')
        if self.proxy_enabled:
            principal = self.access.user(username)
            if principal is None:
                raise AuthorizationDenied('application user is not registered')
            if principal is False:
                raise AuthorizationDenied('application user is disabled')
        now = int(time.time() if now is None else now)
        payload = _urlsafe_encode(json.dumps({
            'username': username,
            'expires': now + self.session_seconds,
        }, sort_keys=True, separators=(',', ':')).encode('utf-8'))
        signature = _urlsafe_encode(hmac.new(
            self.session_secret, payload.encode('ascii'), sha256
        ).digest())
        return '{}.{}'.format(payload, signature)

    def is_public_sample(self, sample_name):
        if not self.proxy_enabled:
            return True
        policy = self.access.sample(sample_name)
        return policy is not None and policy.visibility == 'public'

    def require_authenticated(self, principal):
        if not principal.authenticated:
            raise AuthenticationRequired('authentication is required')
        return principal

    def require_admin(self, principal):
        self.require_authenticated(principal)
        if not principal.is_admin:
            raise AuthorizationDenied('administrator access is required')
        return principal

    def can_access_sample(self, principal, sample_name):
        if not self.proxy_enabled:
            return True
        policy = self.access.sample(
            sample_name,
            None if not principal.authenticated else principal.username,
        )
        if policy is None:
            return False
        if policy.visibility == 'public':
            return True
        if not principal.authenticated:
            return False
        return (
            principal.is_admin
            or policy.owner_username == principal.username
            or policy.granted
        )

    def require_sample(self, principal, sample_name):
        if not self.can_access_sample(principal, sample_name):
            if not principal.authenticated:
                raise AuthenticationRequired('authentication is required')
            raise AuthorizationDenied(
                'access to sample {!r} is not permitted'.format(sample_name)
            )
        return principal

    def _session_username(self, cookie_header, *, now=None):
        token = _cookie_value(cookie_header, 'sdf_session')
        if not token:
            return None
        payload, separator, signature = token.partition('.')
        if not separator:
            return None
        try:
            payload_bytes = payload.encode('ascii')
        except UnicodeEncodeError:
            return None
        expected = _urlsafe_encode(hmac.new(
            self.session_secret, payload_bytes, sha256
        ).digest())
        if not hmac.compare_digest(signature, expected):
            return None
        try:
            values = json.loads(_urlsafe_decode(payload))
            expires = int(values['expires'])
            username = str(values['username']).strip()
        except (KeyError, TypeError, ValueError, UnicodeDecodeError,
                binascii.Error, json.JSONDecodeError):
            return None
        now = int(time.time() if now is None else now)
        if expires < now or not username:
            return None
        if any(character in username for character in '\r\n\0'):
            return None
        return username


def read_session_secret(path):
    """Read a stable shared HMAC secret from a non-public file."""

    value = Path(path).expanduser().resolve().read_bytes().strip()
    if len(value) < 32:
        raise ValueError('session secret must contain at least 32 bytes')
    return value


def _cookie_value(header, name):
    for part in str(header or '').split(';'):
        key, separator, value = part.strip().partition('=')
        if separator and key == name:
            return value
    return None


def _urlsafe_encode(value):
    return base64.urlsafe_b64encode(value).decode('ascii').rstrip('=')


def _urlsafe_decode(value):
    padding = '=' * (-len(value) % 4)
    return base64.urlsafe_b64decode(value + padding).decode('utf-8')
