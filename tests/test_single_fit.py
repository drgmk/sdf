from pathlib import Path

import pytest

from .context import sdf
import sdf.single_fit


def test_prepare_single_fit_package_normalizes_name_and_copies_input(tmp_path):
    source = tmp_path / 'Synthetic star.tbl'
    source.write_text('first version\n', encoding='utf-8')

    package = sdf.single_fit.prepare_single_fit_package(source)

    assert package.directory == tmp_path / 'Synthetic-star-sdf'
    assert package.rawphot_path == (
        package.directory / 'Synthetic-star-rawphot.txt'
    )
    assert package.rawphot_path.read_text() == 'first version\n'


def test_prepare_single_fit_package_refreshes_only_input(tmp_path):
    source = tmp_path / 'star-rawphot.txt'
    source.write_text('first version\n', encoding='utf-8')
    package = sdf.single_fit.prepare_single_fit_package(source)
    retained = package.directory / 'fits' / 'retained.json'
    retained.parent.mkdir()
    retained.write_text('{}\n', encoding='utf-8')
    source.write_text('second version\n', encoding='utf-8')

    repeated = sdf.single_fit.prepare_single_fit_package(source)

    assert repeated == package
    assert package.rawphot_path.read_text() == 'second version\n'
    assert retained.read_text() == '{}\n'


def test_prepare_single_fit_package_rejects_another_rawphot(tmp_path):
    source = tmp_path / 'star.tbl'
    source.write_text('input\n', encoding='utf-8')
    output = tmp_path / 'output'
    output.mkdir()
    (output / 'other-rawphot.txt').write_text('other\n', encoding='utf-8')

    with pytest.raises(
            sdf.utils.SdfError, match='another rawphot input'):
        sdf.single_fit.prepare_single_fit_package(source, output=output)


def test_prepare_single_fit_package_rejects_detached_package_member(tmp_path):
    source = tmp_path / 'star-rawphot.txt'
    source.write_text(
        "\\sdf_package = '.'\n| Band |\n| char |\n",
        encoding='utf-8',
    )

    with pytest.raises(sdf.utils.SdfError, match='belongs to an SDF package'):
        sdf.single_fit.prepare_single_fit_package(source)


@pytest.mark.parametrize(
    ('name', 'expected'),
    (
        ('target-rawphot.txt', 'target'),
        ('target.tbl', 'target'),
        ('Target name.ipac', 'Target-name'),
    ),
)
def test_input_stem(name, expected):
    assert sdf.single_fit._input_stem(name) == expected
