"""Mirror-aware external-service URLs used by the catalog web boundary."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import tomllib
from urllib.parse import quote, urlencode, urlsplit


def _hostname(value, default):
    candidate = str(value or default).strip()
    parsed = urlsplit(
        candidate if '://' in candidate else 'https://' + candidate
    )
    if not parsed.hostname:
        raise ValueError('external service has no hostname: {!r}'.format(
            candidate
        ))
    return '{}:{}'.format(parsed.hostname, parsed.port) \
        if parsed.port is not None else parsed.hostname


@dataclass(frozen=True)
class ExternalServices:
    simbad_server: str = 'simbad.cds.unistra.fr'
    vizier_server: str = 'vizier.cds.unistra.fr'

    @classmethod
    def from_sdb_config(cls, path=None):
        mirrors = {}
        if path is not None:
            config_path = Path(path).expanduser()
            if config_path.is_file():
                with config_path.open('rb') as stream:
                    mirrors = tomllib.load(stream).get('mirrors', {})
        return cls(
            simbad_server=_hostname(
                os.environ.get('SDB_SIMBAD_SERVER') or mirrors.get('simbad'),
                cls.simbad_server,
            ),
            vizier_server=_hostname(
                os.environ.get('SDB_VIZIER_SERVER') or mirrors.get('vizier'),
                cls.vizier_server,
            ),
        )

    def simbad_identifier_url(self, identifier):
        return 'https://{}/simbad/sim-id?{}'.format(
            self.simbad_server,
            urlencode({'submit': 'submit id', 'Ident': str(identifier)}),
        )

    def vizier_catalog_url(self, source):
        return 'https://{}/viz-bin/VizieR?-source={}'.format(
            self.vizier_server, quote(str(source), safe=''),
        )

    def vizier_row_url(self, source, column, value):
        return (
            'https://{}/viz-bin/VizieR-5?'
            '-out.add=.&-source={}&{}==={}'
        ).format(
            self.vizier_server,
            quote(str(source), safe=''),
            quote(str(column), safe=''),
            quote(str(value), safe=''),
        )


def simbad_identifier_route(identifier, root=''):
    return '{}{}?{}'.format(
        str(root or '').rstrip('/'), '/external/simbad/id',
        urlencode({'name': str(identifier)}),
    )
