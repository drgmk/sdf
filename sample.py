"""Generate HTML pages to browse database.

Uses non-standard version of astropy to ensure html anchors retained in
jsviewer tables. This is a virtualenv created into which the modified
version of astropy is installed.

"""

import sys
if sys.version_info[0] < 3:
    raise Exception("Please use python3, not {}!".format(sys.version))

import io
from datetime import datetime
import glob
from os.path import isdir,isfile,basename
from os import mkdir,remove,write,rmdir
import argparse

import numpy as np
from astropy.table import Table
from astropy.utils import xml
import mysql.connector

from jinja2 import Template
import bokeh.resources

from sdf import db
from sdf import plotting
from sdf import tables
from sdf import templates
from sdf import www
from sdf import utils
from sdf import config as cfg


def cleanup_sample_dirs():
    """Remove dirs for samples that no longer exist."""

    dirs = glob.glob(cfg.file['www_root']+'samples/*/')
    samp = db.get_samples()
    for d in dirs:
        dname = basename(d.rstrip('/'))
        if dname not in samp:
            fs = glob.glob(d+'/*')
            fs += glob.glob(d+'/.*')
            print("  {} removed (and files {})".format(d,fs))
            [remove(f) for f in fs]
            rmdir(d)
        else:
            print("  {} ok".format(d))


def cleanup_calibration_dirs():
    """Remove dirs for calibrations that are no longer required."""
    
    fs = glob.glob(cfg.file['www_root']+'calibration/*')
    for f in fs:
        fname = basename(f.rstrip('.html'))
        if fname not in cfg.www['cal_samples']:
            print("  {} removed ".format(f))
            remove(f)
        else:
            print("  {} ok".format(f))


def sample_tables():
    """Generate tables for all samples."""

    # set up connection
    cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                  password=cfg.mysql['passwd'],
                                  host=cfg.mysql['host'],
                                  database=cfg.mysql['db_sdb'])
    cursor = cnx.cursor(buffered=True)

    # get a list of samples and generate their pages
    samples = db.get_samples()
    for sample in samples:
        print("  sample:",sample)
        tables.sample_table_www(cursor,sample)
        tables.sample_table_votable(cursor,sample)

    cursor.close()
    cnx.close()


def sample_plots():

    # set up connection
    cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                  password=cfg.mysql['passwd'],
                                  host=cfg.mysql['host'],
                                  database=cfg.mysql['db_sdb'])
    cursor = cnx.cursor(buffered=True)

    template = Template(templates.sample_plot_wide)
    bokeh_js = bokeh.resources.CDN.render_js()
    bokeh_css = bokeh.resources.CDN.render_css()

    # get a list of samples and generate their pages
    samples = db.get_samples()
    for sample in samples:
        
        print("  sample:",sample)
        
        wwwroot = cfg.file['www_root']+'samples/'

        # create dir and .htaccess if neeeded
        www.create_dir(wwwroot,sample)
        file = wwwroot+sample+"/hr.html"
        script,div = plotting.sample_plot(cursor,sample)

        html = template.render(bokeh_js=bokeh_js,
                               bokeh_css=bokeh_css,
                               css=templates.css,
                               plot_script=script,
                               plot_div=div,
                               title=sample,
                               creation_time=datetime.utcnow().strftime("%d/%m/%y %X"))

        with io.open(file, mode='w', encoding='utf-8') as f:
            f.write(html)

    cursor.close()
    cnx.close()


def flux_size_plots():

    # set up connection
    cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                  password=cfg.mysql['passwd'],
                                  host=cfg.mysql['host'],
                                  database=cfg.mysql['db_sdb'])
    cursor = cnx.cursor(buffered=True)


    template = Template(templates.sample_plot)
    bokeh_js = bokeh.resources.CDN.render_js()
    bokeh_css = bokeh.resources.CDN.render_css()

    # get a list of samples and generate their pages
    samples = db.get_samples()
    for sample in samples:
        
        print("  sample:",sample)
        
        wwwroot = cfg.file['www_root']+'samples/'

        # create dir and .htaccess if neeeded
        www.create_dir(wwwroot,sample)
        file = wwwroot+sample+"/fnuvsr.html"
        out = plotting.flux_size_plot(cursor,sample)
        if out is not None:
            script,div = out
        else:
            continue        

        html = template.render(bokeh_js=bokeh_js,
                               bokeh_css=bokeh_css,
                               css=templates.css,
                               plot_script=script,
                               plot_div=div,
                               title=sample,
                               creation_time=datetime.utcnow().strftime("%d/%m/%y %X"))

        with io.open(file, mode='w', encoding='utf-8') as f:
            f.write(html)
            
    cursor.close()
    cnx.close()


# run from the command line
if __name__ == "__main__":

    # inputs
    parser = argparse.ArgumentParser(description='Update web pages')
    parser.add_argument('--tables','-t',action='store_true',
                        help='Update sample tables')
    parser.add_argument('--plots','-p',action='store_true',
                        help='Update sample plots')
    parser.add_argument('--calibration','-l',action='store_true',
                        help='Update calibration plots')
    parser.add_argument('--cleanup','-c',action='store_true',
                        help='Remove unneccessary sample dirs')
    args = parser.parse_args()

    if args.tables:
        print("Updating sample tables")
        sample_tables()

    if args.plots:
        print("Updating sample plots")
        flux_size_plots()
        sample_plots()

    if args.calibration:
        print("Updating calibration plots")
        for sample in cfg.www['cal_samples']:
            plotting.calibration(sample=sample)

    if args.cleanup:
        print("Cleaning up")
        cleanup_sample_dirs()
        cleanup_calibration_dirs()
