from pathlib import Path

import pytest

from .context import sdf


def test_bundled_toml_contains_complete_defaults(tmp_path, monkeypatch):
    monkeypatch.delenv('SDF_CONFIG', raising=False)

    parser, sources = sdf.config.load_config(home=tmp_path)

    assert sources == (sdf.config._DEFAULT_CONFIG,)
    assert set(parser.sections()) == {
        'calc', 'file', 'models', 'fitting', 'db', 'sqlite', 'mysql',
        'www', 'filter_plot',
    }
    assert parser['calc'].getint('cpu') == 8
    assert parser['models']['star'] == 'phoenix,kurucz,koester,star'
    assert parser['fitting']['models'] == '()'
    assert parser['fitting'].getfloat('ev_threshold') == 5.0
    assert parser['fitting'].getint('n_live') == 200
    assert parser['www'].getint('tablemax') == 5000
    assert parser['filter_plot']['g17'].endswith('SPIRE500')


def test_user_dotfile_overrides_only_selected_defaults(tmp_path, monkeypatch):
    home = tmp_path / 'home'
    home.mkdir()
    user_config = home / '.sdf.conf'
    user_config.write_text(
        '[file]\nmodel_root = models\n'
        '[fitting]\nn_live = 321\n',
        encoding='utf-8',
    )
    monkeypatch.delenv('SDF_CONFIG', raising=False)

    parser, sources = sdf.config.load_config(home=home)

    assert sources[-1] == user_config
    assert parser['file']['model_root'] == str(home / 'models') + '/'
    assert parser['fitting'].getint('n_live') == 321
    assert parser['fitting'].getfloat('ev_threshold') == 5.0


def test_explicit_toml_applies_partial_override(tmp_path):
    config_path = tmp_path / 'override.toml'
    config_path.write_text(
        """
[file]
model_root = "models"
spectra = "spectra"

[models]
star = ["photosphere", "wd"]
default_resolution = 150

[fitting]
models = [["photosphere"], ["photosphere", "disk"]]
ev_threshold = 6.0
n_live = 321
n_update = 12345
n_samples_max = 456
exclude_filters = ["BAD"]
upperlim_filters = ["LIMIT1", "LIMIT2"]
verb = true
""".strip()
        + '\n',
        encoding='utf-8',
    )

    parser, sources = sdf.config.load_config(config_path)

    assert sources[-1] == config_path.resolve()
    assert parser['file']['model_root'] == str(tmp_path / 'models') + '/'
    assert parser['file']['spectra'] == str(tmp_path / 'spectra') + '/'
    assert parser['models']['star'] == 'photosphere,wd'
    assert parser['models'].getfloat('default_resolution') == 150
    assert parser['fitting']['models'] == (
        "(('photosphere',), ('photosphere', 'disk'))"
    )
    assert parser['fitting'].getfloat('ev_threshold') == 6.0
    assert parser['fitting'].getint('n_live') == 321
    assert parser['fitting'].getint('n_update') == 12345
    assert parser['fitting'].getint('n_samples_max') == 456
    assert parser['fitting']['exclude_filters'] == 'BAD'
    assert parser['fitting']['upperlim_filters'] == 'LIMIT1,LIMIT2'
    assert parser['fitting'].getboolean('verb') is True


def test_unknown_toml_setting_is_rejected(tmp_path):
    config_path = tmp_path / 'override.toml'
    config_path.write_text(
        '[fitting]\nlive_point = 200\n', encoding='utf-8'
    )

    with pytest.raises(ValueError, match='live_point'):
        sdf.config.load_config(config_path)


def test_bundled_model_root_is_absolute_and_present(tmp_path, monkeypatch):
    monkeypatch.delenv('SDF_CONFIG', raising=False)

    parser, _ = sdf.config.load_config(home=tmp_path)
    model_root = Path(parser['file']['model_root'])

    assert model_root.is_absolute()
    assert model_root.is_dir()
