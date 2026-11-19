"""Runtime defaults and optional local configuration.

The bundled ``sdf.toml`` is the complete, authoritative default configuration
for both command-line entry points and the Python package. Machine-specific
changes belong in ``~/.sdf.conf``. An explicit TOML or INI override can be
selected with ``SDF_CONFIG``.

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
_DEFAULT_CONFIG = _PACKAGE_ROOT / 'sdf.toml'


def _directory(value, relative_to):
    """Return one expanded absolute directory string with a trailing slash."""

    path = Path(os.path.expandvars(str(value))).expanduser()
    if not path.is_absolute():
        path = Path(relative_to) / path
    return str(path.resolve()) + os.sep


def _parser_value(section, key, value):
    """Convert one TOML value to the historical ConfigParser interface."""

    name = '{}.{}'.format(section, key)
    if isinstance(value, list):
        if section == 'fitting' and key in {'models', 'extra_models'}:
            if not all(
                    isinstance(group, list)
                    and all(isinstance(item, str) for item in group)
                    for group in value):
                raise ValueError(
                    '{} must be an array of arrays of strings'.format(name)
                )
            return repr(tuple(tuple(group) for group in value))
        if not all(isinstance(item, str) for item in value):
            raise ValueError('{} must be an array of strings'.format(name))
        return ','.join(value)
    if isinstance(value, (str, int, float, bool)):
        return str(value)
    raise ValueError('{} has an unsupported TOML value'.format(name))


def _apply_toml(parser, path, *, defaults=False):
    """Apply a complete default TOML file or a partial TOML override."""

    path = Path(path)
    with path.open('rb') as handle:
        payload = tomllib.load(handle)
    if not all(isinstance(table, dict) for table in payload.values()):
        raise ValueError('SDF TOML top-level entries must be tables')

    for section, table in payload.items():
        if section not in parser:
            if not defaults:
                raise ValueError('unknown SDF configuration table: {}'.format(
                    section
                ))
            parser.add_section(section)
        unknown = sorted(set(table) - set(parser[section]))
        if unknown and not defaults:
            raise ValueError(
                'unknown {} setting(s): {}'.format(
                    section, ', '.join(unknown)
                )
            )
        for key, value in table.items():
            parsed = _parser_value(section, key, value)
            if section == 'file':
                parsed = _directory(parsed, path.parent)
            parser[section][key] = parsed


def _read_ini_override(parser, path):
    """Read a concise INI override and resolve any changed file paths."""

    previous_paths = dict(parser['file'])
    parser.read(path)
    for key, value in parser['file'].items():
        if value != previous_paths.get(key):
            parser['file'][key] = _directory(value, Path(path).parent)


def load_config(path=None, *, cwd=None, home=None):
    """Load defaults and return ``(ConfigParser, source_paths)``.

    The complete bundled TOML file is always loaded first. With no explicit
    *path* (or ``SDF_CONFIG``), ``~/.sdf.conf`` is then applied as the concise
    per-user override. Explicit TOML and INI override files are also accepted.
    *cwd* remains accepted for callers of the previous discovery API but is no
    longer used.
    """

    del cwd
    parser = configparser.ConfigParser()
    _apply_toml(parser, _DEFAULT_CONFIG, defaults=True)
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

    home = Path.home() if home is None else Path(home)
    user_config = home / '.sdf.conf'
    if user_config.is_file():
        _read_ini_override(parser, user_config)
        sources.append(user_config)
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
# TODO: move these plotting constants into sdf.toml.

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
