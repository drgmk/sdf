"""Generate www material for results from a single target.""" 

import io
import os
import glob
from datetime import datetime

import numpy as np

import jinja2
from bokeh.resources import CDN,INLINE

from . import plotting
from . import db
from . import utils
from . import config as cfg


def www_all(results, write_path=None, sed_file='index.html',
            f_limits_file='f_limits.html', update=False):
    """Generate www material.
    
    Parameters
    ----------
    results : list of sdf.result.Result
        List of Result objects.
    write_path : str, optional
        Path to write files to, location of photometry file by default.
    sed_file : str, optional
        Name of sed file.
    f_limits_file : str, optional
        Name of f_limits file.
    update : bool, optional
        Force update of html output.
    """

    print(" Web")

    if write_path is None:
        write_path = results[0].path

    file = write_path+'/'+sed_file
    f_file = write_path+'/'+f_limits_file

    # see whether index.html needs updating (unless update enforced)
    if os.path.exists(file):
        mtime = []
        for r in results:
            mtime.append( r.pickle_time )

        if os.path.getmtime(file) > np.max(mtime):
            if not update:
                print("   no update needed")
                return
        print("   updating")
    else:
        print("   generating")

    sed_page(results, file=file)
    f_limits_page(results, file=f_file)


def home_page(file=cfg.file['www_root']+'index.html'):
    """Make the home page.
    
    Parameters
    ----------
    file : str, optional
        File name.
    """

    env = jinja2.Environment(autoescape=False,
         loader=jinja2.PackageLoader('sdf',package_path='www/templates'))
    template = env.get_template("home.html")

    html = template.render()

    with io.open(file, mode='w', encoding='utf-8') as f:
        f.write(html)


def sed_page(results,file='index.html',cdn=True):
    """Make html page with an SED and other info.
        
    Parameters
    ----------
    results : list of sdf.result.Result
        List of Result objects.
    file : str, optional
        File to write html to.
    cdn : bool, optional
        Use CDN for bokeh javascript, otherwise inline.
    """

    script,div = plotting.sed_components(results)

    env = jinja2.Environment(autoescape=False,
         loader=jinja2.PackageLoader('sdf',package_path='www/templates'))
    template = env.get_template("sed_page.html")

    if cdn:
        bokeh_js = CDN.render_js()
        bokeh_css = CDN.render_css()
    else:
        bokeh_js = INLINE.render_js()
        bokeh_css = INLINE.render_css()

    html = template.render(
               js=[bokeh_js],
               css=[bokeh_css],
               plot_script=script,
               plot_div=div,
               phot_file=os.path.basename(results[0].rawphot),
               json_file=os.path.basename(results[0].pmn_dir)+'/'+\
                         os.path.basename(results[0].json),
               main_id=results[0].obs_keywords['main_id'],
               spty=results[0].obs_keywords['sp_type'],
               ra=results[0].obs_keywords['raj2000'],
               dec=results[0].obs_keywords['dej2000'],
               iau_coord=utils.iau_coord(results[0].obs_keywords['raj2000'],
                                         results[0].obs_keywords['dej2000']),
               plx=results[0].obs_keywords['plx_value'],
               xids=db.sdb_xids(results[0].obs_keywords['id']),
               best_fit=results[0].main_results_text(),
               par_dist=os.path.basename(results[0].pmn_dir)+'/'+\
                        os.path.basename(results[0].corner_plot),
               derived_dist=os.path.basename(results[0].pmn_dir)+'/'+\
                        os.path.basename(results[0].distributions_plot),
               h_obsid=','.join(utils.get_herschel_obsid(results[0].obs)),
               creation_time=datetime.utcnow().strftime("%d/%m/%y %X")
                           )

    with io.open(file, mode='w', encoding='utf-8') as f:
        f.write(html)


def f_limits_page(results,file='f_limits.html',cdn=True):

    script,div = plotting.f_limits(results[0])

    env = jinja2.Environment(autoescape=False,
         loader=jinja2.PackageLoader('sdf',package_path='www/templates'))
    template = env.get_template("sed_page.html")

    if cdn:
        bokeh_js = CDN.render_js()
        bokeh_css = CDN.render_css()
    else:
        bokeh_js = INLINE.render_js()
        bokeh_css = INLINE.render_css()

    html = template.render(
               js=[bokeh_js],
               css=[bokeh_css],
               plot_script=script,
               plot_div=div,
               phot_file=os.path.basename(results[0].rawphot),
               json_file=os.path.basename(results[0].pmn_dir)+'/'+\
                         os.path.basename(results[0].json),
               main_id=results[0].obs_keywords['main_id'],
               spty=results[0].obs_keywords['sp_type'],
               ra=results[0].obs_keywords['raj2000'],
               dec=results[0].obs_keywords['dej2000'],
               iau_coord=utils.iau_coord(results[0].obs_keywords['raj2000'],
                                         results[0].obs_keywords['dej2000']),
               plx=results[0].obs_keywords['plx_value'],
               xids=db.sdb_xids(results[0].obs_keywords['id']),
               best_fit=results[0].main_results_text(),
               par_dist=os.path.basename(results[0].pmn_dir)+'/'+\
                        os.path.basename(results[0].corner_plot),
               derived_dist=os.path.basename(results[0].pmn_dir)+'/'+\
                        os.path.basename(results[0].distributions_plot),
               creation_time=datetime.utcnow().strftime("%d/%m/%y %X")
                           )

    with io.open(file, mode='w', encoding='utf-8') as f:
        f.write(html)



def cleanup_sample_dirs():
    """Remove dirs for samples that no longer exist."""

    dirs = glob.glob(cfg.file['www_root']+'samples/*/')
    samp = db.get_samples()
    for d in dirs:
        dname = os.path.basename(d.rstrip('/'))
        if dname not in samp:
            fs = glob.glob(d+'/*')
            fs += glob.glob(d+'/.*')
            print("  {} removed (and files {})".format(d,fs))
            [os.remove(f) for f in fs]
            os.rmdir(d)
        else:
            print("  {} ok".format(d))


def cleanup_calibration_dirs():
    """Remove plots for calibrations that are no longer required."""
    
    fs = glob.glob(cfg.file['www_root']+'calibration/*')
    for f in fs:
        fname = os.path.basename(f.rstrip('.html'))
        if fname not in cfg.www['cal_samples']:
            if not os.path.isdir(f):
                print("  {} removed ".format(f))
                os.remove(f)
        else:
            print("  {} ok".format(f))


def create_dir(wwwroot,sample):
    """Create sample directories and .htaccess if necessary."""

    # make dir and .htaccess if dir doesn't exist
    if not os.path.isdir(wwwroot+sample):
        os.mkdir(wwwroot+sample)

    # make .htaccess if needed, don't put one in those ending with "_"
    # so those directories remain visible to those not logged in
    if  sample[-1] != '_':
        fd = open(wwwroot+sample+'/.htaccess','w')
        fd.write('AuthName "Must login"\n')
        fd.write('AuthType Basic\n')
        fd.write('AuthUserFile '+cfg.file['www_root']+'.htpasswd\n')
        fd.write('AuthGroupFile '+cfg.file['www_root']+'.htgroup\n')
        fd.write('require group admin '+sample+'\n')
        fd.close()
