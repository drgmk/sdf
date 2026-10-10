"""Trusted-proxy identities and legacy Apache group-file authorization."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from threading import RLock


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


class CatalogAuthorizer:
    """Authorize catalog routes from either local or proxy identity."""

    def __init__(self, *, user_header=None, group_file=None,
                 admin_group='admin', local_actor=None):
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
            return CatalogPrincipal(None)
        if any(character in username for character in '\r\n\0'):
            raise AuthenticationRequired('invalid authenticated user header')
        return CatalogPrincipal(
            username, self.group_file.groups_for(username)
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
        if not principal.authenticated:
            return False
        return (
            self.admin_group in principal.groups
            or str(sample_name) in principal.groups
        )

    def require_sample(self, principal, sample_name):
        self.require_authenticated(principal)
        if not self.can_access_sample(principal, sample_name):
            raise AuthorizationDenied(
                'access to sample {!r} is not permitted'.format(sample_name)
            )
        return principal


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
