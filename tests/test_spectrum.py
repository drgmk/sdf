import shutil

import numpy as np
from astropy.io import fits

from .context import sdf


def test_read_phoenix(tmp_path):
    model_root = sdf.config.file['model_root']
    for name in ('lte009.5-3.5-0.0.BT-Settl.7.bz2',
                 'lte058-4.5-0.0a+0.0.BT-Settl.7.bz2'):
        model_file = tmp_path / name
        shutil.copyfile(model_root + name, model_file)

        spectrum = sdf.spectrum.ModelSpectrum.read_phoenix(str(model_file))

        assert spectrum.name == 'phoenix'
        assert (tmp_path / (name + '.npy')).is_file()


def test_cassis_low_resolution_reader_splits_and_labels_modules(tmp_path):
    path = tmp_path / 'cassis.fits'
    rows = []
    for module, wavelengths in (
        (0, np.linspace(7.5, 14.0, 8)),
        (2, np.linspace(20.0, 37.0, 9)),
        (3, np.linspace(14.0, 21.0, 7)),
    ):
        for wavelength in wavelengths:
            row = np.zeros(15, dtype=np.float32)
            row[0] = wavelength
            row[1] = wavelength / 10.0
            row[2] = 0.1
            row[8] = module
            rows.append(row)
    hdu = fits.PrimaryHDU(np.asarray(rows))
    names = (
        'wavelength', 'flux', 'error (RMS+SYS)', 'error (RMS)',
        'error (SYS)', 'offset uncertainty (CAL)', 'sky', 'sky error',
        'module', 'order', 'nod', 'status', 'flag1', 'flag2', 'flag3',
    )
    for index, name in enumerate(names, start=1):
        hdu.header['COL{:02d}DEF'.format(index)] = name
    hdu.header['BUNIT'] = 'Jy'
    hdu.header['CASVE'] = '7'
    hdu.writeto(path)

    spectra = sdf.spectrum.ObsSpectrum.read_cassis_lores(
        path, module_split=True
    )

    assert [value.segment_id for value in spectra] == ['SL1', 'LL2', 'LL1']
    assert [len(value.wavelength) for value in spectra] == [8, 7, 9]
    assert all(value.resolving_power == 100.0 for value in spectra)
