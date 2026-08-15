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


def test_sdf_fits_can_render_operator_extension(monkeypatch, tmp_path):
    calls = []

    class FakeManager:
        def __init__(self, package):
            calls.append(('init', package))

        def render(self, output=None, fits=None, operator=False):
            calls.append(('render', output, fits, operator))
            return tmp_path / 'operator.html'

    monkeypatch.setattr(
        sys, 'argv', ['sdf-fits', 'render', str(tmp_path), '--operator']
    )
    monkeypatch.setattr(sdf.scripts, 'FitManager', FakeManager)

    sdf.scripts.sdf_fits()

    assert calls == [
        ('init', str(tmp_path)), ('render', None, None, True),
    ]


def test_sdf_fits_serve_starts_catalog_worker_and_watcher(monkeypatch, tmp_path):
    calls = []

    monkeypatch.setattr(
        sys, 'argv', [
            'sdf-fits', 'serve', str(tmp_path), '--state',
            str(tmp_path / 'operator.sqlite'), '--actor', 'Grant',
            '--port', '8123',
        ],
    )
    monkeypatch.setattr(
        sdf.scripts, 'FitManager', lambda package: object(),
    )
    monkeypatch.setattr(
        sdf.scripts, 'serve_catalog',
        lambda package, **kwargs: calls.append((package, kwargs)),
    )

    sdf.scripts.sdf_fits()

    assert calls == [(str(tmp_path), {
        'index_path': None,
        'state_path': str(tmp_path / 'operator.sqlite'),
        'actor': 'Grant', 'rebuild': False,
        'host': '127.0.0.1', 'port': 8123,
        'worker': True, 'watch': True, 'autofit_existing': False,
        'poll_seconds': 5.0,
    })]


def test_sdf_sample_uses_indexed_hierarchy_products(monkeypatch, tmp_path):
    calls = []
    index_path = tmp_path / 'sdf-results.sqlite'
    output_root = tmp_path / 'samples'
    monkeypatch.setattr(
        sys, 'argv', [
            'sdf-sample', '--samples', 'reasons', '--tables', '--plots',
            '--index', str(index_path), '--output-root', str(output_root),
            '--sed-base-url', 'https://example.test/seds',
        ],
    )
    monkeypatch.setattr(
        sdf.scripts.db, 'sample_targets', lambda sample: ['HD 123', 'HD 456']
    )
    monkeypatch.setattr(
        sdf.scripts, 'write_sample_products',
        lambda *args, **kwargs: calls.append((args, kwargs)),
    )

    sdf.scripts.sdf_sample()

    assert len(calls) == 1
    args, kwargs = calls[0]
    assert isinstance(args[0], sdf.result_index.ResultIndex)
    assert args[1:] == (
        'reasons', ['HD 123', 'HD 456'], output_root / 'reasons'
    )
    assert kwargs == {
        'table': True, 'plots': True,
        'sed_base_url': 'https://example.test/seds',
    }
