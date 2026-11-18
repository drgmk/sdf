from hashlib import sha256

import numpy as np
import pytest

from .context import sdf
from .test_dataset import _write_empty_rawphot, _write_rawphot
import sdf.candidate_library
import sdf.dataset
import sdf.model_selection
import sdf.spectral_manifest
import sdf.result_index


def _write_csv(path, scale=1.0):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        '\n'.join(
            '{},{},{}'.format(wavelength, scale * flux, 0.1)
            for wavelength, flux in ((8.0, 1.0), (9.0, 1.1), (10.0, 1.2))
        ) + '\n',
        encoding='utf-8',
    )


def test_single_target_package_loads_portable_spectrum(tmp_path):
    _write_rawphot(
        tmp_path / 'A-rawphot.txt', sdbid='A', main_id='HD 123 A'
    )
    spectrum_path = tmp_path / 'spectra' / 'irs.csv'
    _write_csv(spectrum_path)
    digest = sha256(spectrum_path.read_bytes()).hexdigest()
    (tmp_path / 'spectra.yml').write_text(
        'version: 1\n'
        'spectra:\n'
        '  - id: custom:irs\n'
        '    file: spectra/irs.csv\n'
        '    format: csv\n'
        '    contributors: [A]\n'
        '    sha256: {}\n'.format(digest),
        encoding='utf-8',
    )

    dataset = sdf.dataset.FitDataset.from_path(tmp_path)

    assert [row.kind for row in dataset.observations] == [
        'photometry', 'spectrum'
    ]
    spectral = dataset.observations[1]
    assert spectral.observation_id == 'custom:irs:spectrum'
    assert spectral.product_id == 'custom:irs'
    assert spectral.segment_id == 'spectrum'
    assert spectral.contributor_keys == ('simbad:HD 123 A',)
    assert spectral.photometry is None
    assert spectral.spectra[0].fnujy.tolist() == [1.0, 1.1, 1.2]
    assert dataset.manifests == (str(tmp_path / 'spectra.yml'),)


def test_manifest_segment_override_changes_only_one_segment(
        tmp_path, monkeypatch):
    _write_rawphot(
        tmp_path / 'A-rawphot.txt', sdbid='A', main_id='HD 123 A'
    )
    _write_rawphot(
        tmp_path / 'B-rawphot.txt', sdbid='B', main_id='HD 123 B'
    )
    spectrum_path = tmp_path / 'irs.fits'
    spectrum_path.write_bytes(b'test spectrum')
    (tmp_path / 'spectra.yml').write_text(
        'version: 1\n'
        'spectra:\n'
        '  - id: cassis:42\n'
        '    file: irs.fits\n'
        '    format: cassis-irs\n'
        '    contributors: [A, B]\n'
        '    segment_overrides:\n'
        '      SL1:\n'
        '        contributors: [A]\n'
        '      LL1:\n'
        '        include: false\n',
        encoding='utf-8',
    )

    def fake_reader(*args, **kwargs):
        return tuple(
            sdf.spectrum.ObsSpectrum(
                wavelength=np.array([8.0, 9.0]),
                fnujy=np.array([1.0, 1.0]),
                e_fnujy=np.array([0.1, 0.1]),
                segment_id=name,
                instrument='Spitzer IRS ' + name,
            )
            for name in ('SL1', 'LL1')
        )

    monkeypatch.setattr(
        sdf.spectrum.ObsSpectrum, 'read_file_of_type', fake_reader
    )

    dataset = sdf.dataset.FitDataset.from_package(tmp_path)
    spectra = [row for row in dataset.observations if row.kind == 'spectrum']

    assert len(spectra) == 1
    assert spectra[0].observation_id == 'cassis:42:SL1'
    assert spectra[0].contributor_keys == ('simbad:HD 123 A',)


def test_manifest_spectral_projection_sums_declared_contributors(
        tmp_path):
    _write_rawphot(
        tmp_path / 'A-rawphot.txt', sdbid='A', main_id='HD 123 A'
    )
    _write_rawphot(
        tmp_path / 'B-rawphot.txt', sdbid='B', main_id='HD 123 B'
    )
    _write_csv(tmp_path / 'irs.csv')
    (tmp_path / 'spectra.yml').write_text(
        'version: 1\n'
        'spectra:\n'
        '  - id: custom:ab\n'
        '    file: irs.csv\n'
        '    format: csv\n'
        '    contributors: [A, B]\n',
        encoding='utf-8',
    )
    dataset = sdf.dataset.FitDataset.from_package(tmp_path)
    candidate = sdf.model_selection.ModelCandidate.from_owned_models((
        ('simbad:HD 123 A', 'kurucz-0.0'),
        ('simbad:HD 123 B', 'kurucz-0.0'),
    ))

    projection = dataset.component_projection(candidate)

    assert projection.tolist() == [
        [True, False],
        [False, True],
        [True, True],
    ]


def test_spectrum_content_and_manifest_change_input_fingerprint(tmp_path):
    _write_rawphot(tmp_path / 'A-rawphot.txt', sdbid='A')
    spectrum_path = tmp_path / 'irs.csv'
    _write_csv(spectrum_path)
    manifest_path = tmp_path / 'spectra.yml'
    manifest_path.write_text(
        'version: 1\n'
        'spectra:\n'
        '  - id: custom:a\n'
        '    file: irs.csv\n'
        '    format: csv\n'
        '    contributors: [A]\n',
        encoding='utf-8',
    )
    first = sdf.dataset.FitDataset.from_package(tmp_path)
    first_fingerprint = sdf.candidate_library.input_fingerprint(first)

    assert sdf.result_index.package_input_fingerprint(
        tmp_path
    ) == first_fingerprint

    _write_csv(spectrum_path, scale=2.0)
    second = sdf.dataset.FitDataset.from_package(tmp_path)
    second_fingerprint = sdf.candidate_library.input_fingerprint(second)
    manifest_path.write_text(
        manifest_path.read_text(encoding='utf-8') + '# changed\n',
        encoding='utf-8',
    )
    third = sdf.dataset.FitDataset.from_package(tmp_path)

    assert second_fingerprint != first_fingerprint
    assert sdf.candidate_library.input_fingerprint(third) != second_fingerprint


def test_manifest_rejects_unsafe_path_and_unknown_segment(tmp_path,
                                                          monkeypatch):
    _write_rawphot(tmp_path / 'A-rawphot.txt', sdbid='A')
    (tmp_path / 'spectra.yml').write_text(
        'version: 1\n'
        'spectra:\n'
        '  - id: bad\n'
        '    file: ../outside.csv\n'
        '    format: csv\n'
        '    contributors: [A]\n',
        encoding='utf-8',
    )

    with pytest.raises(sdf.utils.SdfError, match='package-relative'):
        sdf.spectral_manifest.read_spectra_manifest(
            tmp_path / 'spectra.yml'
        )


def test_nospec_ignores_manifest_data(tmp_path):
    _write_rawphot(tmp_path / 'A-rawphot.txt', sdbid='A')
    _write_csv(tmp_path / 'irs.csv')
    (tmp_path / 'spectra.yml').write_text(
        'version: 1\n'
        'spectra:\n'
        '  - id: custom:a\n'
        '    file: irs.csv\n'
        '    format: csv\n'
        '    contributors: [A]\n',
        encoding='utf-8',
    )

    dataset = sdf.dataset.FitDataset.from_package(tmp_path, nospec=True)

    assert [row.kind for row in dataset.observations] == ['photometry']


def test_legacy_spectrum_only_rawphot_is_not_discarded(
        tmp_path, monkeypatch):
    rawphot = tmp_path / 'A-rawphot.txt'
    _write_empty_rawphot(rawphot, sdbid='A', main_id='HD 123 A')
    observed = sdf.spectrum.ObsSpectrum(
        wavelength=np.array([8.0, 9.0]),
        fnujy=np.array([1.0, 1.1]),
        e_fnujy=np.array([0.1, 0.1]),
        segment_id='SL1',
    )
    monkeypatch.setattr(
        sdf.spectrum.ObsSpectrum, 'read_sdb_file',
        lambda *args, **kwargs: (observed,),
    )

    dataset = sdf.dataset.FitDataset.from_rawphot(rawphot)

    assert len(dataset.observations) == 1
    assert dataset.observations[0].kind == 'spectrum'
    assert dataset.observations[0].segment_id == 'SL1'
