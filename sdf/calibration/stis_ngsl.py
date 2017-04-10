"""Generate material to compare STIS NGSL spectra to models."""

import os
import io
from datetime import datetime

from jinja2 import Template
import bokeh.resources
import mysql.connector
import astropy.io.fits
import astropy.units as u

from .. import utils
from .. import db
from .. import fitting
from .. import plotting
from .. import tables
from .. import templates
from .. import config as cfg


def add_obs_spec_fits(fig,r,fits=None):
    """Add a spectrum contained in a fits file to a figure.
        
    Parameters
    ----------
    fig : bokeh.plotting.figure
        Figure in which to plot.
    r : sdf.result.Result
        Result to plot data from.
    fits_spec : str
        Name of fits file to plot.
        
    See Also
    --------
    sdf.plotting.sed_components
    sdf.calibration.stis_ngsl.generate_cal_seds
    """
    # TODO: this is STIS NGSL specific for now, could make more general

    hdu = astropy.io.fits.open(fits)
    wave = hdu[1].data['WAVELENGTH'] / 1e4
    flux = hdu[1].data['FLUX'] * u.Unit('erg/(cm2*s*Angstrom)').to('Jy',
                           equivalencies=u.spectral_density(wave*u.micron))

    fig.line(wave,flux,legend='ngsl',**cfg.pl['mod_sp'][-1])


def generate_cal_seds(out_dir=cfg.file['www_root']+'calibration/stis_ngsl/',
                      cdn=True):
    """Generate a set of SEDs with NGSL STIS spectra added.
        
    Parameters
    ----------
    out_dir : string, optional
        Where to put results.
    """

    print("STIS NGSL calibration plots")

    # set up connection
    try:
        cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                      password=cfg.mysql['passwd'],
                                      host=cfg.mysql['host'],
                                      database=cfg.mysql['db_samples'])
        cursor = cnx.cursor(buffered=True)

    except mysql.connector.InterfaceError:
        print("   Can't connect to {} at {}".format(cfg.mysql['db_samples'],
                                                    cfg.mysql['host']) )
        return

    # get the ids and spectra names
    cursor.execute("SELECT sdbid,spec_file FROM stis_ngsl_ "
                   "WHERE sdbid IS NOT NULL")

    for (sdbid,fits_name) in cursor:

        print(" {}".format(sdbid))

        results = fitting.fit_results(utils.rawphot_path(sdbid))

        fits = cfg.file['spectra']+'stis_ngsl/'+fits_name
        script, div = plotting.sed_components(results,extra_func=add_obs_spec_fits,
                                              extra_kwargs={'fits':fits},
                                              model_spec_kwargs={'plot_wave':None})

        template = Template(templates.index)

        if cdn:
            bokeh_js = bokeh.resources.CDN.render_js()
            bokeh_css = bokeh.resources.CDN.render_css()
        else:
            bokeh_js = bokeh.resources.INLINE.render_js()
            bokeh_css = bokeh.resources.INLINE.render_css()

        html = template.render(bokeh_js=bokeh_js,
                               bokeh_css=bokeh_css,
                               css=templates.css,
                               plot_script=script,
                               plot_div=div,
                               main_id=results[0].obs_keywords['main_id'],
                               spty=results[0].obs_keywords['sp_type'],
                               ra=results[0].obs_keywords['raj2000'],
                               dec=results[0].obs_keywords['dej2000'],
                               plx=results[0].obs_keywords['plx_value'],
                               xids=db.sdb_xids(results[0].obs_keywords['id']),
                               best_fit=results[0].main_results_text(),
                               creation_time=datetime.utcnow().strftime("%d/%m/%y %X")
                               )

        with io.open(out_dir+sdbid+'.html', mode='w', encoding='utf-8') as f:
            f.write(html)


def generate_cal_table():
    """Generate a table that helps browse calibration SEDs."""

    print("STIS NGSL calibration table")

    # set up connection
    try:
        cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                      password=cfg.mysql['passwd'],
                                      host=cfg.mysql['host'],
                                      database=cfg.mysql['db_sdb'])
        cursor = cnx.cursor(buffered=True)

    except mysql.connector.InterfaceError:
        print("   Can't connect to {} at {}".format(cfg.mysql['db_sdb'],
                                                    cfg.mysql['host']) )
        return

    tables.sample_table_www(cursor,'stis_ngsl_',absolute_paths=False,
                            file=cfg.file['www_root']+'calibration/stis_ngsl/index.html')
