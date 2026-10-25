import pytest

from .context import sdf
import sdf.resource_manifest


def test_resource_manifest_reads_scoped_product_links(tmp_path):
    path = tmp_path / 'resources.yml'
    path.write_text(
        'version: 1\n'
        'resources:\n'
        '  - provider: tldr_pacs\n'
        '    label: Herschel/PACS images\n'
        '    url: https://example.test/herschel/sources/source.html\n'
        '    observation_ids: [1342223848]\n'
        '    observation_sdbids: [AB]\n'
        '    target_sdbids: [A, B]\n',
        encoding='utf-8',
    )

    assert sdf.resource_manifest.read_resource_manifest(path) == ({
        'provider': 'tldr_pacs',
        'label': 'Herschel/PACS images',
        'url': 'https://example.test/herschel/sources/source.html',
        'observation_ids': ('1342223848',),
        'observation_sdbids': ('AB',),
        'target_sdbids': ('A', 'B'),
    },)


def test_resource_manifest_rejects_unknown_version(tmp_path):
    path = tmp_path / 'resources.yml'
    path.write_text('version: 2\nresources: []\n', encoding='utf-8')

    with pytest.raises(ValueError, match='unsupported resource manifest'):
        sdf.resource_manifest.read_resource_manifest(path)
