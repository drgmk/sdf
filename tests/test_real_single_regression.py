import importlib.util
from pathlib import Path


PATH = (
    Path(__file__).resolve().parents[1] /
    'examples' / 'fitting' / 'real_single_regression.py'
)
SPEC = importlib.util.spec_from_file_location('real_single_regression', PATH)
real_single_regression = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(real_single_regression)


def test_copy_without_package_markers_preserves_science_header(tmp_path):
    source = tmp_path / 'input.txt'
    source.write_text(
        "\\id='target'\n"
        "\\main_id='HD 123'\n"
        "\\sdf_package='.'\n"
        "\\spectra_manifest='spectra.yml'\n"
        '| Band | Phot |\n',
        encoding='utf-8',
    )
    destination = tmp_path / 'copy.txt'

    real_single_regression._copy_without_package_markers(
        source, destination
    )

    text = destination.read_text()
    assert "\\id='target'" in text
    assert "\\main_id='HD 123'" in text
    assert 'sdf_package' not in text
    assert 'spectra_manifest' not in text
    assert real_single_regression._main_id(destination) == 'HD 123'


def test_summary_row_keeps_only_compact_science_differences(tmp_path):
    comparison = {
        'quantities': [
            {'quantity': 'stellar.Teff', 'z': 0.2},
            {'quantity': 'parameters.norm', 'z': -0.4},
            {'quantity': 'parameters.unused', 'z': 9.0},
        ],
        'evidence': {'difference': 0.3, 'z': 0.8},
        'chi_squared': {'difference': -0.1},
    }

    row = real_single_regression._summary_row(
        'HD 123', tmp_path / 'input.txt', comparison
    )

    assert row['z'] == {
        'stellar.Teff': 0.2, 'parameters.norm': -0.4,
    }
    assert row['max_abs_z'] == 0.4
    assert row['evidence_delta'] == 0.3
    assert 'quantities' not in row


def test_case_list_resolves_sdb_export_paths(tmp_path):
    cases = tmp_path / 'cases.txt'
    cases.write_text('# comment\nsdb-one\n\nsdb-two\n', encoding='utf-8')

    paths = real_single_regression._case_list(cases, tmp_path / 'export')

    assert paths == (
        tmp_path / 'export' / 'sdb-one' / 'sdb-one-rawphot.txt',
        tmp_path / 'export' / 'sdb-two' / 'sdb-two-rawphot.txt',
    )
