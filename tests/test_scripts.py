import sys

import pytest

from .context import sdf
import sdf.scripts


@pytest.mark.parametrize(
    ('entrypoint', 'argv', 'expected'),
    (
        (
            sdf.scripts.sdf_index,
            ['sdf-index', '--help'],
            ('rebuild', 'verify', 'disposable SQLite index', 'Examples:'),
        ),
        (
            sdf.scripts.sdf_fits,
            ['sdf-fits', '--help'],
            ('serve', 'select', 'catalog-wide administrative UI', 'Examples:'),
        ),
        (
            sdf.scripts.sdf_jobs,
            ['sdf-jobs', '--help'],
            (
                'list', 'run', 'asynchronous SDF administrative fit queue',
                'Examples:',
            ),
        ),
    ),
)
def test_new_cli_top_level_help_is_operational(
        monkeypatch, capsys, entrypoint, argv, expected):
    monkeypatch.setattr(sys, 'argv', argv)

    with pytest.raises(SystemExit) as exit_info:
        entrypoint()

    assert exit_info.value.code == 0
    output = capsys.readouterr().out
    for text in expected:
        assert text in output


@pytest.mark.parametrize(
    ('argv', 'expected'),
    (
        (
            ['sdf-fits', 'select', '--help'],
            ('does not run a fit', '--default FIT', '--display FIT'),
        ),
        (
            ['sdf-fits', 'fit', '--help'],
            ('exactly one stellar model', 'TARGET=MODEL', 'Examples:'),
        ),
        (
            ['sdf-fits', 'serve', '--help'],
            (
                'Existing unfit packages are not queued',
                '--autofit-existing', '--workers N', '--auth-user-header',
                '--review-base-url URL',
            ),
        ),
        (
            ['sdf-fits', 'reprioritize', '--help'],
            ('never invokes a sampler', '--all', '--threshold LOG_EVIDENCE'),
        ),
        (
            ['sdf-fits', 'refresh-artifacts', '--help'],
            ('self-describing', 'never invokes a sampler', '--all'),
        ),
        (
            ['sdf-fits', 'cleanup', '--help'],
            ('Dry-run is the default', '--apply', '--superseded-days DAYS'),
        ),
    ),
)
def test_sdf_fits_subcommand_help_explains_effects(
        monkeypatch, capsys, argv, expected):
    monkeypatch.setattr(sys, 'argv', argv)

    with pytest.raises(SystemExit) as exit_info:
        sdf.scripts.sdf_fits()

    assert exit_info.value.code == 0
    output = capsys.readouterr().out
    for text in expected:
        assert text in output


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
        'worker': True, 'workers': 1,
        'watch': True, 'autofit_existing': False,
        'poll_seconds': 5.0,
        'auth_user_header': None, 'auth_database': None,
        'auth_session_secret': None, 'auth_cookie_secure': False,
        'sdb_import_config': None,
        'review_base_url': None,
        'access_log': None,
        'activity_timezone': 'UTC',
    })]


def test_sdf_fits_serve_passes_parallel_worker_count(monkeypatch, tmp_path):
    calls = []

    monkeypatch.setattr(
        sys, 'argv', [
            'sdf-fits', 'serve', str(tmp_path), '--workers', '4',
            '--no-watch',
        ],
    )
    monkeypatch.setattr(sdf.scripts, 'FitManager', lambda package: object())
    monkeypatch.setattr(
        sdf.scripts, 'serve_catalog',
        lambda package, **kwargs: calls.append(kwargs),
    )

    sdf.scripts.sdf_fits()

    assert calls[0]['workers'] == 4
    assert calls[0]['worker'] is True


def test_sdf_fits_serve_builds_sdb_import_configuration(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(
        sys, 'argv', [
            'sdf-fits', 'serve', str(tmp_path),
            '--sdb-database', str(tmp_path / 'sdb.sqlite'),
            '--sdb-reference-database', str(tmp_path / 'reference.sqlite'),
            '--sdb-cache-database', str(tmp_path / 'cache.sqlite'),
            '--sdb-config', str(tmp_path / 'sdb.toml'),
            '--sdb-executable', '/opt/sdb',
        ],
    )
    monkeypatch.setattr(
        sdf.scripts, 'serve_catalog',
        lambda package, **kwargs: calls.append(kwargs),
    )

    sdf.scripts.sdf_fits()

    config = calls[0]['sdb_import_config']
    assert config.executable == '/opt/sdb'
    assert config.database == (tmp_path / 'sdb.sqlite').resolve()
    assert config.reference_database == (tmp_path / 'reference.sqlite').resolve()


def test_sdf_fits_reprioritizes_every_package_and_rebuilds_once(
        monkeypatch, tmp_path, capsys):
    packages = (tmp_path / 'one', tmp_path / 'two')
    calls = []

    monkeypatch.setattr(
        sys, 'argv', ['sdf-fits', 'reprioritize', str(tmp_path), '--all'],
    )
    monkeypatch.setattr(
        sdf.scripts, 'discover_packages', lambda root: packages,
    )
    monkeypatch.setattr(
        sdf.scripts, 'reprioritize_package',
        lambda package, **kwargs: calls.append((package, kwargs)) or object(),
    )

    class FakeIndex:
        def __init__(self, path):
            calls.append(('index', path))

        def rebuild(self, root):
            calls.append(('rebuild', root))

    monkeypatch.setattr(sdf.scripts, 'ResultIndex', FakeIndex)

    sdf.scripts.sdf_fits()

    assert [call[0] for call in calls[:2]] == list(packages)
    assert all(call[1]['result_index'] is False for call in calls[:2])
    assert calls[-1] == ('rebuild', tmp_path)
    assert '"reprioritized": 2' in capsys.readouterr().out


def test_sdf_fits_refreshes_numerical_artifacts_for_every_package(
        monkeypatch, tmp_path, capsys):
    packages = (tmp_path / 'one', tmp_path / 'two')
    calls = []

    monkeypatch.setattr(
        sys, 'argv',
        ['sdf-fits', 'refresh-artifacts', str(tmp_path), '--all'],
    )
    monkeypatch.setattr(
        sdf.scripts, 'discover_packages', lambda root: packages,
    )

    class FakeLibrary:
        def __init__(self, package):
            calls.append(package)

        def refresh_artifacts(self):
            return 2

    monkeypatch.setattr(sdf.scripts, 'CandidateLibrary', FakeLibrary)

    sdf.scripts.sdf_fits()

    assert calls == list(packages)
    assert capsys.readouterr().out == (
        '{"packages": 2, "refreshed": 4}\n'
    )


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
