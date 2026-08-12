import numpy as np
from astropy.table import Table

from .context import sdf


def test_nused_counts_non_ignored_measurements():
    photometry = sdf.photometry.Photometry(
        filters=np.array(['2MJ', '2MH', '2MKS']),
        ignore=np.array([False, True, False]),
    )

    assert photometry.nused == 2


def test_zero_measurement_and_uncertainties_are_ignored(tmp_path):
    path = tmp_path / 'zero-rawphot.txt'
    table = Table(
        rows=[('2MJ', 0.0, 0.0, 0.0, 0, 'Jy', 'test', '', 0)],
        names=('Band', 'Phot', 'Err', 'Sys', 'Lim', 'Unit', 'bibcode',
               'Note1', 'exclude'),
    )
    table.write(path, format='ascii.ipac')

    photometry = sdf.photometry.Photometry.read_sdb_file(path)

    assert photometry.ignore.tolist() == [True]
