import pytest

from .context import sdf
import sdf.identity_manifest


def test_read_identity_manifest(tmp_path):
    path = tmp_path / 'identities.yml'
    path.write_text(
        'version: 1\n'
        'targets:\n'
        '  sdb-one:\n'
        '    aliases:\n'
        '      - HD     38\n'
        '      - HIP 473\n',
        encoding='utf-8',
    )

    assert sdf.identity_manifest.read_identity_manifest(path) == {
        'sdb-one': ('HD     38', 'HIP 473'),
    }


def test_identity_manifest_rejects_unknown_fields(tmp_path):
    path = tmp_path / 'identities.yml'
    path.write_text(
        'version: 1\ntargets: {}\nextra: true\n', encoding='utf-8'
    )

    with pytest.raises(ValueError, match='unknown keys: extra'):
        sdf.identity_manifest.read_identity_manifest(path)
