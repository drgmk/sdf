"""Trusted-proxy identities and legacy Apache group-file authorization."""

from __future__ import annotations

import base64
import binascii
from dataclasses import dataclass
from hashlib import sha256
import hmac
import json
from pathlib import Path
import secrets
from threading import RLock
import time


class AuthenticationRequired(PermissionError):
    """The requested catalog view requires an authenticated identity."""


class AuthorizationDenied(PermissionError):
    """The authenticated identity cannot access the requested resource."""


@dataclass(frozen=True)
class CatalogPrincipal:
    """One identity supplied by the trusted reverse proxy."""

    username: str | None
    groups: frozenset[str] = frozenset()

    @property
    def authenticated(self):
        return self.username is not None


class HtGroupFile:
    """Read Apache ``AuthGroupFile`` syntax, reloading after file changes."""

    def __init__(self, path):
        self.path = Path(path).expanduser().resolve()
        self._stamp = None
        self._groups = {}
        self._lock = RLock()
        self._reload()

    def groups_for(self, username):
        self._reload()
        return frozenset(
            group for group, users in self._groups.items()
            if username in users
        )

    def _reload(self):
        with self._lock:
            stat = self.path.stat()
            stamp = (stat.st_mtime_ns, stat.st_size)
            if stamp == self._stamp:
                return
            groups = parse_htgroup(self.path.read_text(encoding='utf-8'))
            self._groups = groups
            self._stamp = stamp


class PublicSampleFile:
    """Read an explicit, fail-closed list of anonymously visible samples."""

    def __init__(self, path):
        self.path = Path(path).expanduser().resolve()
        self._stamp = None
        self._samples = frozenset()
        self._lock = RLock()
        self._reload()

    def contains(self, sample_name):
        self._reload()
        return str(sample_name) in self._samples

    def _reload(self):
        with self._lock:
            stat = self.path.stat()
            stamp = (stat.st_mtime_ns, stat.st_size)
            if stamp == self._stamp:
                return
            self._samples = parse_public_samples(
                self.path.read_text(encoding='utf-8')
            )
            self._stamp = stamp


class CatalogAuthorizer:
    """Authorize catalog routes from either local or proxy identity."""

    def __init__(self, *, user_header=None, group_file=None,
                 admin_group='admin', local_actor=None,
                 public_sample_file=None, session_secret=None,
                 session_seconds=12 * 60 * 60):
        self.user_header = str(user_header or '').strip() or None
        self.admin_group = str(admin_group or '').strip()
        self.local_actor = str(local_actor or '').strip() or 'local'
        if not self.admin_group:
            raise ValueError('admin group is required')
        if self.user_header and group_file is None:
            raise ValueError(
                'a group file is required with a proxy user header'
            )
        if not self.user_header and group_file is not None:
            raise ValueError(
                'a proxy user header is required with a group file'
            )
        self.group_file = (
            HtGroupFile(group_file) if self.user_header else None
        )
        self.public_sample_file = (
            PublicSampleFile(public_sample_file)
            if public_sample_file is not None else None
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
            return CatalogPrincipal(
                self.local_actor, frozenset((self.admin_group,))
            )
        username = str(headers.get(self.user_header, '')).strip()
        if not username:
            username = self._session_username(headers.get('Cookie')) or ''
        if not username:
            return CatalogPrincipal(None)
        if any(character in username for character in '\r\n\0'):
            raise AuthenticationRequired('invalid authenticated user header')
        return CatalogPrincipal(
            username, self.group_file.groups_for(username)
        )

    def issue_session(self, username, *, now=None):
        """Return a signed, expiring token for an Apache-authenticated user."""

        username = str(username or '').strip()
        if not username or any(character in username for character in '\r\n\0'):
            raise AuthenticationRequired('invalid authenticated username')
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
        return (
            self.public_sample_file is not None
            and self.public_sample_file.contains(sample_name)
        )

    def require_authenticated(self, principal):
        if not principal.authenticated:
            raise AuthenticationRequired('authentication is required')
        return principal

    def require_admin(self, principal):
        self.require_authenticated(principal)
        if self.admin_group not in principal.groups:
            raise AuthorizationDenied('administrator access is required')
        return principal

    def can_access_sample(self, principal, sample_name):
        if self.is_public_sample(sample_name):
            return True
        if not principal.authenticated:
            return False
        return (
            self.admin_group in principal.groups
            or str(sample_name) in principal.groups
        )

    def require_sample(self, principal, sample_name):
        if self.is_public_sample(sample_name):
            return principal
        self.require_authenticated(principal)
        if not self.can_access_sample(principal, sample_name):
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


def parse_htgroup(text):
    """Parse ``group: user user`` lines used by Apache group files."""

    groups = {}
    for line_number, raw in enumerate(str(text).splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        group, separator, members = line.partition(':')
        group = group.strip()
        if not separator or not group:
            raise ValueError(
                'invalid group file line {}: expected group: users'.format(
                    line_number
                )
            )
        groups.setdefault(group, set()).update(members.split())
    return {
        group: frozenset(members) for group, members in groups.items()
    }


def parse_public_samples(text):
    """Parse one exact sample name per line, with blank lines and comments."""

    samples = []
    for raw in str(text).splitlines():
        value = raw.strip()
        if not value or value.startswith('#'):
            continue
        samples.append(value)
    return frozenset(samples)


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
