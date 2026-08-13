import numpy as np
import pytest
from astropy.table import Table

from .context import sdf
import sdf.dataset


def _write_rawphot(path, sdbid='sdbid-v3-test', main_id='HD  123',
                   flux=1.0):
    table = Table(
        rows=[('2MJ', flux, 0.1, 0.0, 0, 'Jy', 'test', 'note', 0)],
        names=('Band', 'Phot', 'Err', 'Sys', 'Lim', 'Unit', 'bibcode',
               'Note1', 'exclude'),
    )
    table.meta['keywords'] = {
        'id': {'value': sdbid},
        'main_id': {'value': main_id},
        'sp_type': {'value': 'G2V'},
    }
    table.write(path, format='ascii.ipac')


def _write_empty_rawphot(path, sdbid, main_id):
    table = Table(
        names=('Band', 'Phot', 'Err', 'Sys', 'Lim', 'Unit', 'bibcode',
               'Note1', 'exclude'),
        dtype=('U30', 'f8', 'f8', 'f8', 'i1', 'U20', 'U30', 'U200', 'i1'),
    )
    table.meta['keywords'] = {
        'id': {'value': sdbid},
        'main_id': {'value': main_id},
        'sp_type': {'value': 'DA1'},
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


def test_joint_package_resolves_observations_to_physical_targets(tmp_path):
    _write_rawphot(
        tmp_path / 'AB-rawphot.txt', sdbid='AB', main_id='HD  123'
    )
    _write_rawphot(
        tmp_path / 'B-rawphot.txt', sdbid='B', main_id='HD  123 B'
    )
    _write_rawphot(
        tmp_path / 'A-rawphot.txt', sdbid='A', main_id='HD  123 A'
    )
    (tmp_path / 'joint-fit.yml').write_text(
        'version: 1\nobservations:\n  AB:\n    - A\n    - B\n',
        encoding='utf-8',
    )

    dataset = sdf.dataset.FitDataset.from_package(tmp_path)

    assert [target.key for target in dataset.targets] == [
        'simbad:HD 123 A', 'simbad:HD 123 B'
    ]
    observations = {
        observation.observation_id: observation
        for observation in dataset.observations
    }
    assert observations['A'].contributor_keys == ('simbad:HD 123 A',)
    assert observations['B'].contributor_keys == ('simbad:HD 123 B',)
    assert observations['AB'].contributor_keys == (
        'simbad:HD 123 A', 'simbad:HD 123 B'
    )
    assert observations['AB'].target_key == 'simbad:HD 123'


def test_joint_package_rejects_missing_contributor(tmp_path):
    _write_rawphot(
        tmp_path / 'AB-rawphot.txt', sdbid='AB', main_id='HD  123'
    )
    (tmp_path / 'joint-fit.yml').write_text(
        'version: 1\nobservations:\n  AB:\n    - A\n',
        encoding='utf-8',
    )

    with pytest.raises(sdf.utils.SdfError, match='missing observations: A'):
        sdf.dataset.FitDataset.from_package(tmp_path)


def test_joint_package_allows_contributor_without_resolved_photometry(
        tmp_path):
    _write_rawphot(
        tmp_path / 'AB-rawphot.txt', sdbid='AB', main_id='HD 123'
    )
    _write_rawphot(
        tmp_path / 'A-rawphot.txt', sdbid='A', main_id='HD 123 A'
    )
    _write_empty_rawphot(
        tmp_path / 'B-rawphot.txt', sdbid='B', main_id='HD 123 B'
    )
    (tmp_path / 'joint-fit.yml').write_text(
        'version: 1\nobservations:\n  AB:\n    - A\n    - B\n',
        encoding='utf-8',
    )

    dataset = sdf.dataset.FitDataset.from_package(tmp_path)

    assert [target.key for target in dataset.targets] == [
        'simbad:HD 123 A', 'simbad:HD 123 B'
    ]
    assert [observation.observation_id
            for observation in dataset.observations] == ['A', 'AB']
    assert dataset.observations[1].contributor_keys == (
        'simbad:HD 123 A', 'simbad:HD 123 B'
    )


def test_joint_member_path_discovers_whole_package_once(tmp_path):
    _write_rawphot(
        tmp_path / 'AB-rawphot.txt', sdbid='AB', main_id='HD 123'
    )
    _write_rawphot(
        tmp_path / 'A-rawphot.txt', sdbid='A', main_id='HD 123 A'
    )
    (tmp_path / 'joint-fit.yml').write_text(
        'version: 1\nobservations:\n  AB:\n    - A\n',
        encoding='utf-8',
    )

    sources = sdf.dataset.discover_fit_sources((
        tmp_path / 'AB-rawphot.txt',
        tmp_path / 'A-rawphot.txt',
        tmp_path,
    ))
    dataset = sdf.dataset.FitDataset.from_path(
        tmp_path / 'AB-rawphot.txt'
    )

    assert sources == (tmp_path.resolve(),)
    assert len(dataset.observations) == 2
    assert len(dataset.targets) == 1


def test_ordinary_directory_with_multiple_files_is_not_one_dataset(tmp_path):
    _write_rawphot(tmp_path / 'one-rawphot.txt', sdbid='one')
    _write_rawphot(tmp_path / 'two-rawphot.txt', sdbid='two')

    with pytest.raises(sdf.utils.SdfError, match='contains 2 rawphot files'):
        sdf.dataset.FitDataset.from_path(tmp_path)

    assert sdf.dataset.discover_fit_sources((tmp_path,)) == (
        (tmp_path / 'one-rawphot.txt').resolve(),
        (tmp_path / 'two-rawphot.txt').resolve(),
    )
