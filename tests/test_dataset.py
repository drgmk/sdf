import numpy as np
from astropy.table import Table

from .context import sdf
import sdf.dataset


def _write_rawphot(path):
    table = Table(
        rows=[('2MJ', 1.0, 0.1, 0.0, 0, 'Jy', 'test', 'note', 0)],
        names=('Band', 'Phot', 'Err', 'Sys', 'Lim', 'Unit', 'bibcode',
               'Note1', 'exclude'),
    )
    table.meta['keywords'] = {
        'id': {'value': 'sdbid-v3-test'},
        'main_id': {'value': 'HD  123'},
        'sp_type': {'value': 'G2V'},
    }
    table.write(path, format='ascii.ipac')


def test_single_rawphot_becomes_self_mapped_dataset(tmp_path):
    rawphot = tmp_path / 'target-rawphot.txt'
    _write_rawphot(rawphot)

    dataset = sdf.dataset.FitDataset.from_rawphot(rawphot)

    assert len(dataset.targets) == 1
    assert dataset.targets[0].key == 'simbad:HD 123'
    assert dataset.targets[0].main_id == 'HD  123'
    assert dataset.targets[0].sdbid == 'sdbid-v3-test'
    assert len(dataset.observations) == 1
    assert dataset.observations[0].contributor_keys == ('simbad:HD 123',)
    assert dataset.observations[0].source == str(rawphot.resolve())
    assert dataset.data[0].filters.tolist() == ['2MJ']
    assert dataset.primary_keywords['sp_type'] == 'G2V'


def test_target_identity_falls_back_to_sdbid():
    identity = sdf.dataset.TargetIdentity.from_keywords({'id': 'sdb-test'})

    assert identity.key == 'sdb:sdb-test'
    assert identity.main_id is None
