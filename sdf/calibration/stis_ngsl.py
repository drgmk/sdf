"""Generate material to compare STIS NGSL spectra to models."""

import os
import io
from datetime import datetime

import numpy as np
import jinja2
from bokeh.plotting import ColumnDataSource
import bokeh.resources
import bokeh.palettes
import mysql.connector
import astropy.io.fits
import astropy.units as u

from .. import utils
from .. import db
from .. import photometry
from .. import filter
from .. import fitting
from .. import plotting
from .. import tables
from .. import result
from .. import config as cfg


def add_obs_spec_fits(fig,r,fits=None,parameters=None):
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

    if parameters is not None:
        label = 'ngsl:'+str(parameters)
    else:
        label = 'ngsl'

    fig.line(wave,flux,legend_label=label,**cfg.pl['mod_sp'][-1])


def add_fixed_model(fig,r,model,parameters=None):
    """Add a model with specific parameters.

    Parameters
    ----------
    fig : bokeh.plotting.figure
        Figure in which to plot.
    r : sdf.result.Result
        A result object, used for info.
    model : tuple
        Model to plot.
    parameters : tuple of lists
        Model paramters, excluding normalisations.
    """

    if parameters is not None:
        rf = result.FixedResult(r.rawphot,model,parameters)
        plotting.add_model_spec(fig,rf)


def add_filters(fig,r):
    """Add filters to a figure.
        
    Parameters
    ----------
    fig : bokeh.plotting.figure
        Figure in which to plot.
    r : sdf.result.Result
        Result to plot data from.

    See Also
    --------
    sdf.plotting.sed_components
    sdf.calibration.stis_ngsl.generate_cal_seds
    """

    c_micron = u.micron.to(u.Hz,equivalencies=u.spectral())
    cols = bokeh.palettes.Category10[10]

    # get filters for this plot, converting colours to filters
    fs = []
    for obs in r.obs:
        if isinstance(obs,photometry.Photometry):
            for fname in obs.filters:
                if filter.iscolour(fname):
                    col = filter.Colour.get(fname)
                    fs.extend(col.filters)
                else:
                    fs.append(fname)

    for i,fname in enumerate(fs):

        f = filter.Filter.get(fname)

        if f.mean_wavelength > 1:
            continue

        data = {'wave': c_micron/f.nu_hz,
                'response':-3 + 6 * f.response/np.max(f.response)}
        pldata = ColumnDataSource(data=data)

        fig.line('wave','response',source=pldata,
                 line_color=cols[i % 10],line_width=2)


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
    cursor.execute("SELECT sdbid,spec_file,IFNULL(stis_Teff,0), "
                   "stis_logg,stis_logz FROM stis_ngsl_ WHERE sdbid IS NOT NULL")

    for (sdbid,fits_name,teff,logg,logz) in cursor:

        print("{}".format(sdbid))

        results = fitting.fit_results(utils.rawphot_path(sdbid))

        print(" Plotting")
        
        # parameters for fixed model
        if teff > 0:
            mod = ('phoenix_m',)
            par = ([teff,logg,logz],)
            rf = result.FixedResult(results[0].rawphot,mod,par)
        else:
            par = None

        print(par)

        fits = cfg.file['spectra']+'stis_ngsl/'+fits_name

        # plot both best-fit and 'true' parameter models
        script, div = plotting.sed_components(
                      results,
                      main_extra_func=(add_obs_spec_fits,
                                       add_fixed_model),
                      main_extra_kwargs=({'fits':fits,
                                         'parameters':par},
                                         {'model':mod,
                                         'parameters':par}),
                      res_extra_func=add_filters,
                      model_spec_kwargs={'plot_wave':None}
                                              )

        # plot only the model with the 'true' parameters
#        script, div = plotting.sed_components(
#                        [rf],
#                        main_extra_func=(add_obs_spec_fits,),
#                        main_extra_kwargs=({'fits':fits},),
#                        res_extra_func=add_filters,
#                        model_spec_kwargs={'plot_wave':None}
#                                              )

        env = jinja2.Environment(autoescape=False,
             loader=jinja2.PackageLoader('sdf',package_path='www/templates'))
        template = env.get_template("sed_page.html")

        if cdn:
            bokeh_js = bokeh.resources.CDN.render_js()
            bokeh_css = bokeh.resources.CDN.render_css()
        else:
            bokeh_js = bokeh.resources.INLINE.render_js()
            bokeh_css = bokeh.resources.INLINE.render_css()

        html = template.render(js=[bokeh_js],
                               css=[bokeh_css],
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
    """Generate a table that helps browse calibration SEDs.
    
    .. todo:: fix broken www links.
    """

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


def generate_cal_hr_diag(cdn=True):
    """Generate HR diagram and f-r plot to help browse calibration SEDs.
    
    .. todo:: fix www so link goes back to calibration table.
    .. todo:: fix thumbs
    """

    print("STIS NGSL HR diagram")

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

    file = cfg.file['www_root']+'calibration/stis_ngsl/hr.html'
    script,div = plotting.sample_plot(cursor,'stis_ngsl_',
                                      absolute_paths=False)

    env = jinja2.Environment(autoescape=False,
         loader=jinja2.PackageLoader('sdf',package_path='www/templates'))
    template = env.get_template("sample_plot.html")

    if cdn:
        bokeh_js = bokeh.resources.CDN.render_js()
        bokeh_css = bokeh.resources.CDN.render_css()
    else:
        bokeh_js = bokeh.resources.INLINE.render_js()
        bokeh_css = bokeh.resources.INLINE.render_css()

    html = template.render(js=[bokeh_js],
                           css=[bokeh_css],
                           body_class='wide',
                           plot_script=script,
                           plot_div=div,
                           title='stis_ngsl_',
                           creation_time=datetime.utcnow().strftime("%d/%m/%y %X"))

    with io.open(file, mode='w', encoding='utf-8') as f:
        f.write(html)

