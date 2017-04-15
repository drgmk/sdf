"""File locations and other configuration.
    
This is mostly configured with a .conf file, but some post-processing is
done here. The imformation contained herein can be accessed from 
elsewhere by importing this module:

>>> import config as cfg
>>> wave = cfg.models['default_wave']

Plotting things that are harder (or more laborious) to configure so are
done fully here for now. Most of these are dictionaries with a list of
keywords to be given to the plotting routines. An (incomplete) example
would look something like:

>>> form bokeh.plotting import figure
>>> import config as cfg
>>> fig = figure()
>>> fig.circle('x','y',source=data,**cfg.pl['obs_ph'])
"""

import os
import glob
import configparser

import numpy as np
import bokeh.palettes
import binarytree as bt

from . import utils

# info from the config files, defaults are in same dir as this file
cfg = configparser.ConfigParser()
cfg.read([
          os.path.dirname(os.path.realpath(__file__))+'/sdf.conf',
          os.path.expanduser('~')+'/.sdf.conf'
         ])

# simple dictionaries with only strings
file = cfg['file']
mysql = cfg['mysql']

# calculation stuff
calc = {
    'cpu':cfg['calc'].getint('cpu')
}

# find all the model names and derive the locations
if os.path.exists(cfg['file']['model_root']):
    model_names = [i for i in os.walk(cfg['file']['model_root'])][0][1]
else:
    raise utils.SdfError("model directory {} doesn't exist".\
                         format(cfg['file']['model_root']))
model_loc = {}
for name in model_names:
    model_loc[name] = file['model_root']+name+'/'

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

models = {
    'names': model_names,
    'loc': model_loc,
    'default_wave': 10**np.arange(np.log10(cfg['models'].getfloat('min_wave_micron')),
                                  np.log10(cfg['models'].getfloat('max_wave_micron')),
                    np.log10(1+1/float(cfg['models'].getfloat('default_resolution')))),
    'min_wav_micron': cfg['models'].getfloat('min_wave_micron'),
    'max_wav_micron': cfg['models'].getfloat('max_wave_micron'),
    'star': star,
    'disk': disk,
    'disk_r': disk_r
}

# Solar radius at 1pc as base unit for emitting area in sr
ssr = 1.596074069110538e-15 #( np.pi*u.solRad**2 / (u.pc)**2 ).si.value

# binary tree for model execution and selection
star_type = ('phoenix_m',)
t = bt.Node('start')
t.left = bt.Node( star_type )
t.right = bt.Node( (star_type+('amsil_r',) ) )
t.right.left = bt.Node( star_type+('amsil_r',) )
t.right.right = bt.Node( star_type+('amsil_r','amsil_r') )
t.right.left.left = bt.Node( star_type+('amsil_r',) )
t.right.left.right = bt.Node( star_type+('modbb_disk_r',) )
t.right.right.left = bt.Node( star_type+('amsil_r','amsil_r') )
t.right.right.right = bt.Node( star_type+('modbb_disk_r','modbb_disk_r') )

# details for fitting
fitting = {
    'models':(
#              ('kurucz_m',),
#        ('phoenix_m','amsil_r'),
#        ('phoenix_m','bb_disk_r'),
    ),
    'tree': t,
    'exclude_filters': cfg['fitting']['exclude_filters'].split(','),
    'upperlim_filters': cfg['fitting']['upperlim_filters'].split(','),
    'model_join': cfg['fitting']['model_join'],
    'pmn_dir_suffix': cfg['fitting']['pmn_dir_suffix'],
    'pmn_model_suffix': cfg['fitting']['pmn_model_suffix'],
    'ev_threshold': cfg['fitting'].getfloat('ev_threshold'),
    'n_live': cfg['fitting'].getint('n_live'),
    'n_update': cfg['fitting'].getint('n_update'),
    'verb': cfg['fitting'].getboolean('verb'),
    'n_samples': cfg['fitting'].getint('n_samples'),
    'model_om_range': (cfg['fitting'].getfloat('model_om_lo'),
                       cfg['fitting'].getfloat('model_om_hi')),
    'spectra_norm_range': (cfg['fitting'].getfloat('spectra_norm_lo'),
                           cfg['fitting'].getfloat('spectra_norm_hi')),
    'mn_oldest': cfg['fitting'].getfloat('mn_oldest'),
    'an_oldest': cfg['fitting'].getfloat('an_oldest')
}

# www stuff
www = {
    'site_root': cfg['www']['site_root'],
    'tablemax': cfg['www'].getint('tablemax'),
    'votmax': cfg['www'].getint('votmax'),
    'plotmax': cfg['www'].getint('plotmax'),
    'cal_samples': cfg['www']['cal_samples'].split(',')
}

# some numbers to be used in the plotting dict below
# TODO: tidy and move to sdf.conf

# colours for each model, first is total model
model_colours = ['navy','firebrick','green','red','darkslategray']
phot_alpha = [0.5,0.3,0.3,0.3,0.3]
line_alpha = [0.3,0.3,0.3,0.3,0.5]

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
    'obs_ph':     {'fill_color':'#444444','fill_alpha':0.8,'size':ob_sz},
    'obs_ig_ph':  {'fill_color':'#444444','fill_alpha':0.4,'line_alpha':0.4,'size':ob_sz},
    'obs_e_ph':   {'line_color':'#444444','line_alpha':0.8,'line_width':line_thin},
    'obs_e_ig_ph':{'line_color':'#444444','line_alpha':0.4,'line_width':line_thin},
    'obs_sp':     {'line_color':'#444444','line_alpha':0.5,'line_width':line_thin},
    'obs_e_sp':   {'line_color':'#444444','line_alpha':0.2,'line_width':line_thin},

    # lines in residuals plot
    'guide_dash': {'line_width':line_thin,'line_alpha':0.5,'line_dash':'dashed'},

    # models
    'mod_ph':[{'size':ph_sz,'line_color':model_colours[0],
               'line_alpha':phot_alpha[0],'fill_alpha':fill_alpha,
               'line_width':line_thin},
              {'size':ph_sz,'line_color':model_colours[1],
               'line_alpha':phot_alpha[1],'fill_alpha':fill_alpha,
               'line_width':line_thin},
              {'size':ph_sz,'line_color':model_colours[2],
               'line_alpha':phot_alpha[2],'fill_alpha':fill_alpha,
               'line_width':line_thin},
              {'size':ph_sz,'line_color':model_colours[3],
               'line_alpha':phot_alpha[3],'fill_alpha':fill_alpha,
               'line_width':line_thin},
              {'size':ph_sz,'line_color':model_colours[4],
               'line_alpha':phot_alpha[4],'fill_alpha':fill_alpha,
               'line_width':line_thin}],
    
    'mod_sp':[{'line_color':model_colours[0],'line_alpha':line_alpha[0],
                    'line_width':line_thick},
              {'line_color':model_colours[1],'line_alpha':line_alpha[1],
                    'line_width':line_thin},
              {'line_color':model_colours[2],'line_alpha':line_alpha[2],
                    'line_width':line_thin},
              {'line_color':model_colours[3],'line_alpha':line_alpha[3],
                    'line_width':line_thin},
              {'line_color':model_colours[4],'line_alpha':line_alpha[4],
                    'line_width':line_thin}],

    # HR diagram #

    'hr_dot':   {'size':10,'fill_alpha':0.6,'line_alpha':1},
    'hr_e_dot': {'line_alpha':0.4,'line_width':line_thin}
}

# filters for filter plot
filter_plot = {
                'groups': [cfg['filter_plot'].get(k).split(',')
                           for k in cfg['filter_plot'].keys()]
}


# for convenience
tiny = 1e-99
huge = 1e99
