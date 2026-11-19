"""Runtime defaults and optional local configuration.

The bundled ``sdf.conf`` remains the complete internal default set. Normal
installations only need an ``sdf.toml`` when they use an external model
library or want to change fitting controls. The configuration can be selected
with ``SDF_CONFIG``; otherwise ``~/.config/sdf/config.toml`` and
``./sdf.toml`` are applied in that order. The old ``~/.sdf.conf`` is still
read as a legacy fallback when ``SDF_CONFIG`` is not set.

The processed values can be accessed by importing this module:

>>> from sdf import config as cfg
>>> wave = cfg.models['default_wave']

Plotting things that are harder (or more laborious) to configure so are
done fully here for now. Most of these are dictionaries with a list of
keywords to be given to the plotting routines. An (incomplete) example
would look something like: 

>>> from bokeh.plotting import figure
>>> import config as cfg
>>> fig = figure()
>>> fig.circle('x', 'y', source=data, **cfg.pl['obs_ph'])
"""

import ast
import configparser
import os
from pathlib import Path
import tomllib

import numpy as np

_PACKAGE_ROOT = Path(__file__).resolve().parent
_DEFAULT_CONFIG = _PACKAGE_ROOT / 'sdf.conf'


def _directory(value, relative_to):
    """Return one expanded absolute directory string with a trailing slash."""

    path = Path(os.path.expandvars(str(value))).expanduser()
    if not path.is_absolute():
        path = Path(relative_to) / path
    return str(path.resolve()) + os.sep


def _string_list(value, name):
    if not isinstance(value, list) or not all(
            isinstance(item, str) for item in value):
        raise ValueError('{} must be an array of strings'.format(name))
    return ','.join(value)


def _number_pair(value, name):
    if (not isinstance(value, list) or len(value) != 2
            or not all(isinstance(item, (int, float)) for item in value)):
        raise ValueError('{} must be an array of two numbers'.format(name))
    return value


def _table(payload, name):
    value = payload.get(name, {})
    if not isinstance(value, dict):
        raise ValueError('{} must be a TOML table'.format(name))
    return value


def _check_keys(table, allowed, name):
    unknown = sorted(set(table) - set(allowed))
    if unknown:
        raise ValueError(
            'unknown {} setting(s): {}'.format(name, ', '.join(unknown))
        )


def _apply_toml(parser, path):
    """Apply the small supported TOML surface to legacy internal defaults."""

    with Path(path).open('rb') as handle:
        payload = tomllib.load(handle)
    _check_keys(payload, {'models', 'fitting', 'paths'}, 'top-level')

    models_table = _table(payload, 'models')
    _check_keys(
        models_table,
        {
            'root', 'stellar_patterns', 'disk_patterns',
            'radial_disk_patterns', 'min_wavelength_micron',
            'max_wavelength_micron', 'resolution',
        },
        'models',
    )
    if 'root' in models_table:
        parser['file']['model_root'] = _directory(
            models_table['root'], Path(path).parent
        )
    for source, target in (
            ('stellar_patterns', 'star'),
            ('disk_patterns', 'disk'),
            ('radial_disk_patterns', 'disk_r')):
        if source in models_table:
            parser['models'][target] = _string_list(
                models_table[source], 'models.{}'.format(source)
            )
    for source, target in (
            ('min_wavelength_micron', 'min_wave_micron'),
            ('max_wavelength_micron', 'max_wave_micron'),
            ('resolution', 'default_resolution')):
        if source in models_table:
            parser['models'][target] = str(models_table[source])

    fitting_table = _table(payload, 'fitting')
    _check_keys(
        fitting_table,
        {
            'evidence_threshold', 'live_points', 'update_interval',
            'posterior_samples', 'model_normalization',
            'spectrum_normalization', 'excluded_filters',
            'upper_limit_filters', 'verbose',
        },
        'fitting',
    )
    for source, target in (
            ('evidence_threshold', 'ev_threshold'),
            ('live_points', 'n_live'),
            ('update_interval', 'n_update'),
            ('posterior_samples', 'n_samples_max'),
            ('verbose', 'verb')):
        if source in fitting_table:
            parser['fitting'][target] = str(fitting_table[source])
    for source, low, high in (
            ('model_normalization', 'model_om_lo', 'model_om_hi'),
            ('spectrum_normalization', 'spectra_norm_lo', 'spectra_norm_hi')):
        if source in fitting_table:
            values = _number_pair(
                fitting_table[source], 'fitting.{}'.format(source)
            )
            parser['fitting'][low] = str(values[0])
            parser['fitting'][high] = str(values[1])
    for source, target in (
            ('excluded_filters', 'exclude_filters'),
            ('upper_limit_filters', 'upperlim_filters')):
        if source in fitting_table:
            parser['fitting'][target] = _string_list(
                fitting_table[source], 'fitting.{}'.format(source)
            )

    paths_table = _table(payload, 'paths')
    _check_keys(paths_table, {'spectra'}, 'paths')
    if 'spectra' in paths_table:
        parser['file']['spectra'] = _directory(
            paths_table['spectra'], Path(path).parent
        )


def _read_ini_override(parser, path):
    """Read an old INI override and make its model path location-independent."""

    previous_root = parser['file']['model_root']
    parser.read(path)
    if parser['file']['model_root'] != previous_root:
        parser['file']['model_root'] = _directory(
            parser['file']['model_root'], Path(path).parent
        )


def load_config(path=None, *, cwd=None, home=None):
    """Load defaults and return ``(ConfigParser, source_paths)``.

    An explicit *path*, or ``SDF_CONFIG`` when *path* is omitted, replaces
    normal user/project discovery. INI paths are accepted to ease migration;
    new configuration should use TOML.
    """

    parser = configparser.ConfigParser()
    parser.read(_DEFAULT_CONFIG)
    parser['file']['model_root'] = _directory(
        'data/models', _PACKAGE_ROOT
    )
    sources = [_DEFAULT_CONFIG]
    explicit = path if path is not None else os.environ.get('SDF_CONFIG')
    if explicit:
        explicit = Path(explicit).expanduser().resolve()
        if not explicit.is_file():
            raise FileNotFoundError(
                'SDF configuration not found: {}'.format(explicit)
            )
        if explicit.suffix.lower() == '.toml':
            _apply_toml(parser, explicit)
        else:
            _read_ini_override(parser, explicit)
        sources.append(explicit)
        return parser, tuple(sources)

    cwd = Path.cwd() if cwd is None else Path(cwd)
    home = Path.home() if home is None else Path(home)
    legacy = home / '.sdf.conf'
    if legacy.is_file():
        _read_ini_override(parser, legacy)
        sources.append(legacy)
    for candidate in (
            home / '.config' / 'sdf' / 'config.toml',
            cwd / 'sdf.toml'):
        if candidate.is_file():
            _apply_toml(parser, candidate)
            sources.append(candidate.resolve())
    return parser, tuple(sources)


cfg, config_sources = load_config()

# simple dictionaries with only strings
file = cfg['file']

# databases
db = cfg['db']
if db['type'] == 'mysql': 
    db.update(cfg['mysql'])
elif db['type'] == 'sqlite': 
    db.update(cfg['sqlite'])

# calculation stuff
calc = {
    'cpu': cfg['calc'].getint('cpu')
}

# find all the model names and derive the locations
model_root = Path(cfg['file']['model_root'])
if model_root.is_dir():
    model_names = sorted(path.name for path in model_root.iterdir()
                         if path.is_dir())
else:
    model_names = []

model_loc = {}
for name in model_names:
    model_loc[name] = _directory(name, model_root)

# classify each model
star_names = cfg['models']['star'].split(',')
disk_names = cfg['models']['disk'].split(',')
disk_r_names = cfg['models']['disk_r'].split(',')
star = []
disk = []
disk_r = []
for n in model_names: 
    for m in star_names: 
        if m in n: 
            star.append(n)
    for m in disk_names: 
        if m in n: 
            disk.append(n)
    for m in disk_r_names: 
        if m in n: 
            disk_r.append(n)

# default wavelengths
default_wave = 10**np.arange(np.log10(cfg['models'].getfloat('min_wave_micron')),
                             np.log10(cfg['models'].getfloat('max_wave_micron')),
                np.log10(1+1/float(cfg['models'].getfloat('default_resolution'))))
if np.max(default_wave) != cfg['models'].getfloat('max_wave_micron'): 
    default_wave = np.append(default_wave, cfg['models'].getfloat('max_wave_micron'))

models = {
    'names': model_names,
    'loc': model_loc,
    'default_wave': default_wave,
    'default_resolution': cfg['models'].getfloat('default_resolution'),
    'min_wav_micron': cfg['models'].getfloat('min_wave_micron'),
    'max_wav_micron': cfg['models'].getfloat('max_wave_micron'),
    'star': star,
    'disk': disk,
    'disk_r': disk_r
}

# Solar radius at 1pc as base unit for emitting area in sr
ssr = 1.596074069110538e-15  # ( np.pi*u.solRad**2 / (u.pc)**2 ).si.value

# details for fitting
fitting = {
    'models': ast.literal_eval(cfg['fitting']['models']),
    'extra_models': ast.literal_eval(cfg['fitting']['extra_models']),
    'exclude_filters': cfg['fitting']['exclude_filters'].split(','),
    'upperlim_filters': cfg['fitting']['upperlim_filters'].split(','),
    'model_join': cfg['fitting']['model_join'],
    'pmn_dir_suffix': cfg['fitting']['pmn_dir_suffix'],
    'pmn_model_suffix': cfg['fitting']['pmn_model_suffix'],
    'ev_threshold': cfg['fitting'].getfloat('ev_threshold'),
    'n_live': cfg['fitting'].getint('n_live'),
    'n_update': cfg['fitting'].getint('n_update'),
    'verb': cfg['fitting'].getboolean('verb'),
    'n_samples_max': cfg['fitting'].getint('n_samples_max'),
    'model_om_range': (cfg['fitting'].getfloat('model_om_lo'),
                       cfg['fitting'].getfloat('model_om_hi')),
    'spectra_norm_range': (cfg['fitting'].getfloat('spectra_norm_lo'),
                           cfg['fitting'].getfloat('spectra_norm_hi')),
    'mn_oldest': cfg['fitting'].getfloat('mn_oldest'),
    'an_oldest': cfg['fitting'].getfloat('an_oldest')
}

# www stuff
www = {
    'base_url': cfg['www']['base_url'],
    'sdb_path': cfg['www']['sdb_path'],
    'sdb_url': cfg['www']['base_url'] + cfg['www']['sdb_path'],
    'tablemax': cfg['www'].getint('tablemax'),
    'votmax': cfg['www'].getint('votmax'),
    'plotmax': cfg['www'].getint('plotmax'),
    'cal_samples': cfg['www']['cal_samples'].split(','),
    'min_cal_wav': cfg['www'].getfloat('min_cal_wav'),
    'max_cal_wav': cfg['www'].getfloat('max_cal_wav')
}

# some numbers to be used in the plotting dict below
# TODO: tidy and move to sdf.conf

# colours for each model, first is total model
model_colours = ['navy', 'firebrick', 'green', 'red', 'darkslategray']
phot_alpha = [0.5, 0.3, 0.3, 0.3, 0.3]
line_alpha = [0.3, 0.3, 0.3, 0.3, 0.5]

line_thin = 2
line_thick = 4
fill_alpha = 0.1
ob_sz = 10
ph_sz = ob_sz+line_thin+2

# plotting dict
pl = {

    # SEDs #

    # standard plot size
    'x_size': 850,
    'y_top_size': 550,
    'y_bot_size': 130,

    # photometry
    'obs_ph':     {'fill_color': '#444444', 'fill_alpha': 0.8, 'size': ob_sz},
    'obs_ig_ph':  {'fill_color': '#444444', 'fill_alpha': 0.4, 'line_alpha': 0.4, 'size': ob_sz},
    'obs_e_ph':   {'line_color': '#444444', 'line_alpha': 0.8, 'line_width': line_thin},
    'obs_e_ig_ph': {'line_color': '#444444', 'line_alpha': 0.4, 'line_width': line_thin},
    'obs_sp':     {'line_color': '#444444', 'line_alpha': 0.5, 'line_width': line_thin},
    'obs_e_sp':   {'line_color': '#444444', 'line_alpha': 0.2, 'line_width': line_thin},

    # lines in residuals plot
    'guide_dash': {'line_width': line_thin, 'line_alpha': 0.5, 'line_dash': 'dashed'},

    # models
    'mod_ph': [{'size': ph_sz, 'line_color': model_colours[0],
               'line_alpha': phot_alpha[0], 'fill_alpha': fill_alpha,
               'line_width': line_thin},
              {'size': ph_sz, 'line_color': model_colours[1],
               'line_alpha': phot_alpha[1], 'fill_alpha': fill_alpha,
               'line_width': line_thin},
              {'size': ph_sz, 'line_color': model_colours[2],
               'line_alpha': phot_alpha[2], 'fill_alpha': fill_alpha,
               'line_width': line_thin},
              {'size': ph_sz, 'line_color': model_colours[3],
               'line_alpha': phot_alpha[3], 'fill_alpha': fill_alpha,
               'line_width': line_thin},
              {'size': ph_sz, 'line_color': model_colours[4],
               'line_alpha': phot_alpha[4], 'fill_alpha': fill_alpha,
               'line_width': line_thin}],
    
    'mod_sp': [{'line_color': model_colours[0], 'line_alpha': line_alpha[0],
                    'line_width': line_thick},
              {'line_color': model_colours[1], 'line_alpha': line_alpha[1],
                    'line_width': line_thin},
              {'line_color': model_colours[2], 'line_alpha': line_alpha[2],
                    'line_width': line_thin},
              {'line_color': model_colours[3], 'line_alpha': line_alpha[3],
                    'line_width': line_thin},
              {'line_color': model_colours[4], 'line_alpha': line_alpha[4],
                    'line_width': line_thin}],

    # HR diagram #
    'hr_dot':   {'size': 10, 'fill_alpha': 0.6, 'line_alpha': 1},
    'hr_e_dot': {'line_alpha': 0.4, 'line_width': line_thin},
    'hr_track': {'line_color': model_colours[0], 'line_alpha': 0.1,
                 'line_width': line_thin},
    'fvsr_join': {'line_color': 'darkslategray', 'line_dash': 'dashed',
                 'line_alpha': 0.2, 'line_width': line_thin}
}

# filters for filter plot
filter_plot = {
                'groups': [cfg['filter_plot'].get(k).replace('\n', '').split(',')
                           for k in cfg['filter_plot'].keys()]
}


# for convenience
tiny = 1e-99
huge = 1e99
