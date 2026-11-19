from pathlib import Path

import pytest

from .context import sdf


def test_explicit_toml_maps_supported_runtime_settings(tmp_path):
    config_path = tmp_path / 'sdf.toml'
    config_path.write_text(
        """
[models]
root = "models"
stellar_patterns = ["photosphere", "wd"]
resolution = 150

[fitting]
evidence_threshold = 6.0
live_points = 321
update_interval = 12345
posterior_samples = 456
model_normalization = [-8.0, 7.0]
spectrum_normalization = [0.5, 2.0]
excluded_filters = ["BAD"]
upper_limit_filters = ["LIMIT1", "LIMIT2"]
verbose = true

[paths]
spectra = "spectra"
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
    assert parser['fitting'].getfloat('ev_threshold') == 6.0
    assert parser['fitting'].getint('n_live') == 321
    assert parser['fitting'].getint('n_update') == 12345
    assert parser['fitting'].getint('n_samples_max') == 456
    assert parser['fitting'].getfloat('model_om_lo') == -8.0
    assert parser['fitting'].getfloat('model_om_hi') == 7.0
    assert parser['fitting'].getfloat('spectra_norm_lo') == 0.5
    assert parser['fitting'].getfloat('spectra_norm_hi') == 2.0
    assert parser['fitting']['exclude_filters'] == 'BAD'
    assert parser['fitting']['upperlim_filters'] == 'LIMIT1,LIMIT2'
    assert parser['fitting'].getboolean('verb') is True


def test_project_toml_overrides_user_toml_and_legacy_ini(tmp_path, monkeypatch):
    home = tmp_path / 'home'
    project = tmp_path / 'project'
    user_config = home / '.config' / 'sdf' / 'config.toml'
    user_config.parent.mkdir(parents=True)
    project.mkdir()
    (home / '.sdf.conf').write_text(
        '[file]\nmodel_root = legacy-models\n', encoding='utf-8'
    )
    user_config.write_text(
        '[models]\nroot = "user-models"\n', encoding='utf-8'
    )
    project_config = project / 'sdf.toml'
    project_config.write_text(
        '[models]\nroot = "project-models"\n', encoding='utf-8'
    )
    monkeypatch.delenv('SDF_CONFIG', raising=False)

    parser, sources = sdf.config.load_config(cwd=project, home=home)

    assert parser['file']['model_root'] == (
        str(project / 'project-models') + '/'
    )
    assert sources[-3:] == (
        home / '.sdf.conf', user_config.resolve(), project_config.resolve()
    )


def test_unknown_toml_setting_is_rejected(tmp_path):
    config_path = tmp_path / 'sdf.toml'
    config_path.write_text(
        '[fitting]\nlive_point = 200\n', encoding='utf-8'
    )

    with pytest.raises(ValueError, match='live_point'):
        sdf.config.load_config(config_path)


def test_bundled_model_root_is_absolute_and_present(tmp_path):
    config_path = tmp_path / 'sdf.toml'
    config_path.write_text('', encoding='utf-8')

    parser, _ = sdf.config.load_config(config_path)
    model_root = Path(parser['file']['model_root'])

    assert model_root.is_absolute()
    assert model_root.is_dir()
