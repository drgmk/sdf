import sys

import pytest

from .context import sdf
import sdf.scripts


def test_joint_fit_allows_www_and_writes_page(tmp_path, monkeypatch):
    (tmp_path / 'joint-fit.yml').write_text(
        'version: 1\nobservations: {}\n', encoding='utf-8'
    )
    rawphot = tmp_path / 'target-rawphot.txt'
    rawphot.write_text('test input\n', encoding='utf-8')
    calls = []

    monkeypatch.setattr(
        sys, 'argv', ['sdf-fit', '--file', str(rawphot), '--www']
    )
    monkeypatch.setattr(
        sdf.scripts, 'fit_package',
        lambda source, **kwargs: calls.append(('fit', source, kwargs)),
    )
    monkeypatch.setattr(
        sdf.scripts, 'write_joint_results_page',
        lambda source: calls.append(('page', source)),
    )

    sdf.scripts.sdf_fit()

    assert calls == [
        ('fit', str(tmp_path), {'force': False, 'nospec': False}),
        ('page', str(tmp_path)),
    ]


@pytest.mark.parametrize('option', ('--dbwrite', '--update-db'))
def test_joint_fit_rejects_legacy_database_output(
        tmp_path, monkeypatch, option):
    (tmp_path / 'joint-fit.yml').write_text(
        'version: 1\nobservations: {}\n', encoding='utf-8'
    )
    rawphot = tmp_path / 'target-rawphot.txt'
    rawphot.write_text('test input\n', encoding='utf-8')

    monkeypatch.setattr(
        sys, 'argv', ['sdf-fit', '--file', str(rawphot), option]
    )

    with pytest.raises(NotImplementedError, match='result index'):
        sdf.scripts.sdf_fit()
