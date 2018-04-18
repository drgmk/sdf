import io
from datetime import datetime
from os import path,remove

import ast
import numpy as np
import matplotlib.pyplot as plt
import mysql.connector
import astropy.units as u

import jinja2
import bokeh.resources
from bokeh.plotting import figure,ColumnDataSource
from bokeh.models.widgets import Panel,Tabs
from bokeh.models import HoverTool,TapTool,OpenURL
from bokeh.layouts import gridplot,layout
import bokeh.palettes
import bokeh.embed

from . import analytics
from . import model
from . import spectrum
from . import photometry
from . import filter
from . import fitting
from . import db
from . import www
from . import utils
from . import config as cfg


def sed_components(results,tab_order=None,
                   main_extra_func=None,main_extra_kwargs={},
                   res_extra_func=None,res_extra_kwargs={},
                   model_spec_kwargs={}):
    """Return bokeh script/div components for an sed.
        
    Parameters
    ----------
    results : list of sdf.result.Result
        List of Results.
    tab_order : list, optional
        Order in which the list of Results should appear in tabs.
    main_extra_func : tuple of functions, optional
        Function that takes a figure instance, a list of results, and
        optionally arbitrary extra keywords. Adds to the main plot.
    main_extra_kwargs : tuple of dicts, optional
        Dict of keywords to be passed to main_extra_func.
    res_extra_func : function, optional
        Function that takes a figure instance, a list of results, and
        optionally arbitrary extra keywords. Adds to the residuals plot.
    res_extra_kwargs : dict, optional
        Dict of keywords to be passed to res_extra_func.
    model_spec_kwargs : dict, optional
        Dict of keywords to be passed to add_model_spec.
    """

    # results should be list so we can loop
    if not isinstance(results,list):
        raise utils.SdfError("expected list of results, not {}".format(type(results)))
    
    # figure plot limits
    xlim,ylim = sed_limits(results)

    # loop to create tabs in plot
    sed = []
    res = []
    tabs = []
    for i,r in enumerate(results):

        if hasattr(r,'evidence'):
            title = "Evidence: {:.2f} | ".format(r.evidence)
        else:
            title = ''

        for j,par in enumerate(r.parameters):
            if 'norm' not in par:
                title += "{}:{:.1f}, ".format(par,r.best_params[j])

        # sed plot
        hover = HoverTool(names=['phot'],tooltips=[('band',"@filter"),
                                                   ('meas',"@flux"),
                                                   ('unc',"@err")])
        tools = ['wheel_zoom,box_zoom,save,reset',hover]

        sed.append( figure(title=title,
                           x_axis_label='Wavelength / micron',x_axis_type='log',
                           x_range=xlim,y_range=ylim,
                           y_axis_label='Flux density / Jy',y_axis_type='log',
                           tools=tools,active_scroll='wheel_zoom',
                           width=cfg.pl['x_size'],height=cfg.pl['y_top_size']) )

        # add models
        add_model_phot(sed[i],r)
        add_model_spec(sed[i],r,**model_spec_kwargs)

        # add obs
        add_obs_phot(sed[i],r)
        add_obs_spec(sed[i],r)

        # optional extra function to add other stuff
        if main_extra_func:
            for func,kw in zip(main_extra_func,main_extra_kwargs):
                func(sed[i],r,**kw)

        # residuals
        hover = HoverTool(names=['resid'],tooltips=[('band',"@filter"),
                                                    ('sigma',"@res")])
        tools = ['wheel_zoom,box_zoom,save,reset',hover]

        res.append( figure(x_axis_type='log',x_axis_location='above',
                           x_range=sed[i].x_range,
                           y_axis_label='Residuals',
                           tools=tools,active_scroll='wheel_zoom',
                           toolbar_location='right',
                           width=cfg.pl['x_size'],height=cfg.pl['y_bot_size']) )

        # add residuals
        res[i].line(x=xlim,y=[0 , 0],**cfg.pl['guide_dash'])
        res[i].line(x=xlim,y=[-3,-3],**cfg.pl['guide_dash'])
        res[i].line(x=xlim,y=[3 , 3],**cfg.pl['guide_dash'])
        add_res(res[i],r)

        if res_extra_func:
            res_extra_func(res[i],r,**res_extra_kwargs)

        grid = gridplot([[sed[i]],[res[i]]],toolbar_location='below')

        tabs.append( Panel(child=grid, title=r.model_info['name']) )
    
    if tab_order is not None:
        if len(tab_order) != len(results):
         raise SdfError("tab_order {} not the same length as results ({})".
                        format(tab_order,len(results)))
        tabs = [tabs[i] for i in tab_order]
    tab = Tabs(tabs=tabs)

    return bokeh.embed.components(tab)


def add_obs_phot(fig,r):
    """Add observed photometry to an SED plot.
        
    Parameters
    ----------
    fig : bokeh.plotting.figure
        Figure in which to plot.
    r : sdf.result.Result
        Result to plot data from.
        
    Notes
    -----
    Only data in photometry.Photometry objects is plotted.
    
    See Also
    --------
    sdf.plotting.sed_components
    """
    
    for p in r.obs:
        if not isinstance(p,photometry.Photometry):
            continue
        
        # photometry, not upper limits or ignored
        data = {}
        ok = np.invert(np.logical_or(p.upperlim,p.ignore))
        data['filter'] = p.filters[ok]
        data['wave'] = p.mean_wavelength()[ok]
        data['flux'] = p.fnujy[ok]
        data['err'] = p.e_fnujy[ok]
        pldata = ColumnDataSource(data=data)
        xs,_,_,err_ys = utils.plot_err(data['wave'],data['wave'],
                                       data['flux'],data['err'])
        fig.multi_line(xs,err_ys,**cfg.pl['obs_e_ph'])
        fig.circle('wave','flux',source=pldata,name='phot',
                   **cfg.pl['obs_ph'])

        # ignored photometry
        ok = np.logical_and(p.ignore,np.invert(p.upperlim))
        data['filter'] = p.filters[ok]
        data['wave'] = p.mean_wavelength()[ok]
        data['flux'] = p.fnujy[ok]
        data['err'] = p.e_fnujy[ok]
        pldata = ColumnDataSource(data=data)
        xs,_,_,err_ys = utils.plot_err(data['wave'],data['wave'],
                                       data['flux'],data['err'])
        fig.multi_line(xs,err_ys,**cfg.pl['obs_e_ig_ph'])
        fig.circle('wave','flux',source=pldata,name='phot',
                   **cfg.pl['obs_ig_ph'])

        # upper limits
        ok = np.logical_and(p.upperlim,np.invert(p.ignore))
        data['filter'] = p.filters[ok]
        data['wave'] = p.mean_wavelength()[ok]
        data['flux'] = p.fnujy[ok]
        data['err'] = p.e_fnujy[ok]
        pldata = ColumnDataSource(data=data)
        fig.inverted_triangle('wave','flux',source=pldata,name='phot',
                              **cfg.pl['obs_ph'])

        # ignored upper limits
        ok = np.logical_and(p.upperlim,p.ignore)
        data['filter'] = p.filters[ok]
        data['wave'] = p.mean_wavelength()[ok]
        data['flux'] = p.fnujy[ok]
        data['err'] = p.e_fnujy[ok]
        pldata = ColumnDataSource(data=data)
        fig.inverted_triangle('wave','flux',source=pldata,name='phot',
                              **cfg.pl['obs_ig_ph'])


def add_obs_spec(fig,r):
    """Add observed spectra to an SED plot.
        
    Parameters
    ----------
    fig : bokeh.plotting.figure
        Figure in which to plot.
    r : sdf.result.Result
        Result to plot data from.
        
    Notes
    -----
    Only data in spectrum.ObsSpectrum objects is plotted.
    
    See Also
    --------
    sdf.plotting.sed_components
    """
    ispec = -1
    for s in r.obs:
        if not isinstance(s,spectrum.ObsSpectrum):
            continue
        data = {}
        data['wave'] = s.wavelength
        data['flux'] = s.fnujy * r.best_params[ispec]
        data['loerr'] = (s.fnujy - s.e_fnujy) * r.best_params[ispec]
        data['hierr'] = (s.fnujy + s.e_fnujy) * r.best_params[ispec]
        ispec -= 1
        pldata = ColumnDataSource(data=data)

        fig.line('wave','flux',source=pldata,**cfg.pl['obs_sp'])
        fig.line('wave','loerr',source=pldata,**cfg.pl['obs_e_sp'])
        fig.line('wave','hierr',source=pldata,**cfg.pl['obs_e_sp'])


def add_model_phot(fig,r):
    """Add model photometry to an SED plot.

    Parameters
    ----------
    fig : bokeh.plotting.figure
        Figure in which to plot.
    r : sdf.result.Result
        Result to plot data from.
        
    See Also
    --------
    sdf.plotting.sed_components
    """

    # filters and colours (spectra have None for filter)
    filt = np.array([isinstance(f,(str,np.str_)) for f in r.filters])
    
    # individual model components, ignoring colours (if more than one)
    if r.n_comps > 1:
        
        # indices that are filters but not colours/indices
        notcol = np.array(filt)
        for i,f in enumerate(r.filters[filt]):
            if filter.iscolour(f):
                notcol[i] = False

        for i in range(r.n_comps):
            data = {}
            data['filter'] = r.filters[notcol]
            data['wave'] = filter.mean_wavelength(data['filter'])
            data['flux'] = r.model_comp_fnujy[i,notcol]

            pldata = ColumnDataSource(data=data)

            fig.circle('wave','flux',source=pldata, **cfg.pl['mod_ph'][i+1])

    # full models
    data = {}
    data['filter'] = r.filters[filt]
    data['wave'] = filter.mean_wavelength(data['filter'])
    data['flux'] = r.model_fnujy[filt]
    pldata = ColumnDataSource(data=data)

    fig.circle('wave','flux',source=pldata, **cfg.pl['mod_ph'][0])


def add_model_spec(fig,r,plot_wave=cfg.models['default_wave']):
    """Add model spectra to an SED plot.
        
    Parameters
    ----------
    fig : bokeh.plotting.figure
        Figure in which to plot.
    r : sdf.result.Result
        Result to plot data from.
    plot_wave : numpy.ndarray, optional
        Wavelengths to use for plot, spectra are not resampled if None.
    
    See Also
    --------
    sdf.plotting.sed_components
    """

    if len(r.comp_spectra) == 0:
        return
    
    # one thick line if one model
    elif len(r.comp_spectra) == 1:
        
        s = r.comp_spectra[0].copy()
        if plot_wave is not None:
            s.resample(plot_wave)
        
        data = {'wave':s.wavelength,
                'flux':s.fnujy}
        pldata = ColumnDataSource(data=data)
        fig.line('wave','flux',source=pldata,**cfg.pl['mod_sp'][0])

    # thin lines plus thick total for >1 models
    else:
        for i,s in enumerate(r.comp_spectra):

            data = {}
            if plot_wave is not None:
                s.resample(plot_wave)
            data['wave'] = s.wavelength
            data['flux'] = s.fnujy
            pldata = ColumnDataSource(data=data)
            
            fig.line('wave','flux',source=pldata,**cfg.pl['mod_sp'][i+1])

        # sum of all models
        s = r.total_spec.copy()
        if plot_wave is not None:
            s.resample(plot_wave)
        data = {'wave':s.wavelength,'flux':s.fnujy}
        pldata = ColumnDataSource(data=data)
        fig.line('wave','flux',source=pldata,**cfg.pl['mod_sp'][0])


def add_res(fig,r):
    """Add residuals to a plot.
        
    Parameters
    ----------
    fig : bokeh.plotting.figure
        Figure in which to plot.
    r : sdf.result.Result
        Result to plot data from.
        
    See Also
    --------
    sdf.plotting.sed_components
    """

    # used photometry
    data = {}
    ok = np.invert(r.filters_ignore)
    data['filter'] = r.filters[ok]
    data['wave'] = r.wavelengths[ok]
    data['res'] = r.residuals[ok]
    pldata = ColumnDataSource(data=data)
    fig.circle('wave','res',source=pldata,name='resid',**cfg.pl['obs_ph'])

    # ignored photometry
    ok = r.filters_ignore
    data['filter'] = r.filters[ok]
    data['wave'] = r.wavelengths[ok]
    data['res'] = r.residuals[ok]
    pldata = ColumnDataSource(data=data)
    fig.circle('wave','res',source=pldata,name='resid',**cfg.pl['obs_ig_ph'])


def hardcopy_sed(r,file=None,fig=None,xsize=8,ysize=6,dpi=100,
                 axis_labels=True):
    """Make a hardcopy SED for a specific result object.

    Parameters
    ----------
    r : sdf.result.Result object
        Result object with details to be plotted.
    file : str, optional
        Name of file to write plot to.
    fig : matplotlib.pyplot.Figure object
        Overplot on this figure.
    xsize : float, optional
        X size of plot.
    ysize : float, optional
        Y size of plot.
    dpi : int, optional
        DPI of output figure.
    axis_labels : bool, optional
        Show axis labels or not.

    .. todo: this is very simple, needs work.
    """

    if fig is None:
        fig,ax = plt.subplots(figsize=(xsize,ysize))
    else:
        ax = fig.axes[0]

    # model spectra
    for s in r.comp_spectra:
        ax.loglog(s.wavelength,s.fnujy)

    ax.loglog(r.total_spec.wavelength,r.total_spec.fnujy)


    # photometry
    for p in r.obs:
        if not isinstance(p,photometry.Photometry):
            continue

        ok = np.invert(np.logical_or(p.upperlim,p.ignore))
        ax.errorbar(p.mean_wavelength()[ok],p.fnujy[ok],yerr=p.e_fnujy[ok],
                    fmt='o')

    # cosmetics
    xl,yl = sed_limits((r,))
    ax.set_xlim(0.3,2e3)
    ax.set_ylim(yl)
    if axis_labels:
        ax.set_xlabel('Wavelength / $\mu$m')
        ax.set_ylabel('Flux density / Jy')
    else:
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])

    if file is not None:
        fig.tight_layout()
        fig.savefig(file,dpi=dpi)
        plt.close(fig)

    return fig


def sed_limits(results):
    """Figure out plotting limits given observations."""

    pad = 1.5
    padlg = 2.
    xlims = [-1,-1]
    ylims = [-1,-1]
    
    for r in results:
        obs = r.obs
        for o in obs:

            if isinstance(o,photometry.Photometry):
                xlim = [np.min(o.mean_wavelength()),np.max(o.mean_wavelength())]
            elif isinstance(o,spectrum.ObsSpectrum):
                xlim = [np.min(o.wavelength),np.max(o.wavelength)]

            if xlim[0] < xlims[0] or xlims[0] == -1:
                xlims[0] = xlim[0]
            if xlim[1] > xlims[1] or xlims[1] == -1:
                xlims[1] = xlim[1]

            ylim = [np.min(o.fnujy),np.max(o.fnujy)]

            # set lower y limit to uncertainty if negative or zero fluxes
            if ylim[0] <= 0:
                ylim[0] = np.min(o.e_fnujy)

            if ylim[0] < ylims[0] or ylims[0] == -1:
                ylims[0] = ylim[0]
            if ylim[1] > ylims[1] or ylims[1] == -1:
                ylims[1] = ylim[1]

        for comp in r.comp_spectra:
            if ylims[1] < np.max(comp.fnujy) or ylims[1] == -1:
                ylims[1] = np.max(comp.fnujy)

    xlims = [xlims[0]/pad,  xlims[1]*pad]
    ylims = [ylims[0]/padlg,ylims[1]*pad]
    return xlims,ylims


def f_limits(r):
    '''Make a plot of what could have been detected.'''

    cols = bokeh.palettes.Category20_20
    ncols = len(cols)

    # check we have a star
    if hasattr(r,'star'):

        # and a distance
        if 'plx_arcsec' in r.star[0].keys():
            dist = 1/r.star[0]['plx_arcsec']
            lstar = np.sum( [s['lstar'] for s in r.star] )
            bb = analytics.BB_Disk(lstar=lstar, distance=dist)
        else:
            bb = analytics.BB_Disk()

        # observed limits
        lim,fname = bb.f_limits_from_result(r)

        # "theoretical" photometric limit
        th_lim = 0.05 * \
                 bb.temperatures**3 / r.star[0]['Teff']**3
        
        # approximate ALMA limit
        lstar_1pc = np.sum( [s['lstar_1pc'] for s in r.star] )
        alma_lim = bb.f_limits(np.array([880]),
                               flux_limits=np.array([20e-6]),
                               fwhm=np.array([1]),
                               lstar_1pc=lstar_1pc)

        yrange = [2e-7,0.1]

        hover = HoverTool(names=['lim'],tooltips=[('band',"@filter")])
        tools = ['wheel_zoom,box_zoom,save,reset',hover]

        # vs disk temperature
        fig = figure(title='detection limits',
                     x_axis_label='Temperature / Kelvin',x_axis_type='log',
                     x_range=[0.99*np.min(bb.temperatures),
                              1.01*np.max(bb.temperatures)],
                     y_range=yrange,
                     y_axis_label='Fractional luminosity',y_axis_type='log',
                     tools=tools,active_scroll='wheel_zoom',
                     width=cfg.pl['x_size'],height=cfg.pl['y_top_size'] )

        for i in range(lim.shape[1]):
            data = {'temp':bb.temperatures,'f':lim[:,i],
                    'filter':np.repeat(fname[i],len(lim[:,i]))}
            pldata = ColumnDataSource(data=data)
            fig.line('temp','f',source=pldata,color=cols[i % ncols],
                     line_width=3,muted_alpha=0.2,
                     name='lim',legend=fname[i])

        # extra limits
        for l,n,s in zip([th_lim,alma_lim],
                         ['Phot lim 5%','1mm/1" 20uJy'],
                         ['solid', 'dashed']):
            data = {'temp':bb.temperatures,'f':l,
                    'filter':np.repeat(n, len(l))}
            pldata = ColumnDataSource(data=data)
            fig.line('temp','f',source=pldata,color='lightgrey',
                     line_width=3,muted_alpha=0.2, line_dash=s,
                     name='lim',legend=n)

        fig.legend.click_policy="mute"
        fig.legend.location = 'top_left'
        fig.legend.label_text_font_size = '10pt'

        tab1 = Panel(child=fig, title='vs. temperature')

        # vs disk radius, need a stellar luminosity for this
        if bb.lstar is not None:
            fig = figure(title='detection limits',
                         x_axis_label='Blackbody radius / au',x_axis_type='log',
                         x_range=[np.min(bb.blackbody_radii()),
                                  np.max(bb.blackbody_radii())],
                         y_range=yrange,
                         y_axis_label='Fractional luminosity',y_axis_type='log',
                         tools=tools,active_scroll='wheel_zoom',
                         width=cfg.pl['x_size'],height=cfg.pl['y_top_size'] )

            for i in range(lim.shape[1]):
                data = {'rad':bb.blackbody_radii(),'f':lim[:,i],
                        'filter':np.repeat(fname[i],len(lim[:,i]))}
                pldata = ColumnDataSource(data=data)
                fig.line('rad','f',source=pldata,color=cols[i % ncols],
                         line_width=3,muted_alpha=0.2,
                         name='lim',legend=fname[i])

            # extra limits
            for l,n,s in zip([th_lim,alma_lim],
                             ['Phot lim 5%','1mm/1" 20uJy'],
                             ['solid', 'dashed']):
                data = {'rad':bb.blackbody_radii(),'f':l,
                        'filter':np.repeat('limit',len(l))}
                pldata = ColumnDataSource(data=data)
                fig.line('rad','f',source=pldata,color='lightgrey',
                         line_width=3,muted_alpha=0.2, line_dash=s,
                         name='lim',legend=n)

            fig.legend.click_policy="mute"
            fig.legend.location = 'top_left'
            fig.legend.label_text_font_size = '10pt'

            tab2 = Panel(child=fig, title='vs. radius')
        else:
            fig = figure(title='L_star unknown, radius limits not available',
                         width=cfg.pl['x_size'],height=cfg.pl['y_top_size'])
            tab2 = Panel(child=fig, title='vs. blackbody radius')

        if bb.lstar:
            tab = Tabs(tabs=[tab2,tab1])
        else:
            tab = Tabs(tabs=[tab1,tab2])

        return bokeh.embed.components(tab)

    else:
        fig = figure(title='limits could not be derived')
        return bokeh.embed.components(fig)


def sample_plots():

    # set up connection
    cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                  password=cfg.mysql['passwd'],
                                  host=cfg.mysql['host'],
                                  database=cfg.mysql['db_sdb'])
    cursor = cnx.cursor(buffered=True)

    env = jinja2.Environment(autoescape=False,
         loader=jinja2.PackageLoader('sdf',package_path='www/templates'))
    template = env.get_template("sample_plot.html")

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
        script,div = sample_plot(cursor,sample)

        html = template.render(
                   js=[bokeh_js],
                   css=[bokeh_css],
                   body_class='wide',
                   plot_script=script,
                   plot_div=div,
                   title=sample,
                   creation_time=datetime.utcnow().strftime("%d/%m/%y %X")
                               )

        with io.open(file, mode='w', encoding='utf-8') as f:
            f.write(html)

    cursor.close()
    cnx.close()


def sample_plot(cursor,sample,absolute_paths=True):
    """Return bokeh componenets for HR diagram + f vs. r sample plots.

    Parameters
    ----------
    cursor : mysql.connector.connect.cursor
        Connection to database.
    sample : string
        Name of sample in config.mysql['db_samples'] to plot.
    absolute_paths : bool, optional
        Use absolute urls, otherwise relative to this plot.
    """

    # get data, ensure primary axes are not nans else bokeh will
    # complain about the data
    if sample == 'everything' or sample == 'public':
        sel1 = " FROM "+cfg.mysql['db_sdb']+".sdb_pm"
    else:
        sel1 = (" FROM "+cfg.mysql['db_samples']+"."+sample+
                " LEFT JOIN "+cfg.mysql['db_sdb']+".sdb_pm USING (sdbid)")

    sel = sel1 + (" LEFT JOIN "+cfg.mysql['db_sdb']+".simbad USING (sdbid)"
                  " LEFT JOIN "+cfg.mysql['db_results']+".star ON sdbid=star.id"
                  " LEFT JOIN "+cfg.mysql['db_results']+".disk_r ON sdbid=disk_r.id")

    # statement for selecting stuff to plot
    sel += " WHERE teff IS NOT NULL AND lstar IS NOT NULL"

    selall = ("SELECT sdbid,main_id,teff,e_teff,lstar,e_lstar,"
              " IFNULL(ldisk_lstar,-1) as ldisklstar,"
              " IFNULL(e_ldisk_lstar,-1) as e_ldisklstar,"
              " IFNULL(temp,-1) as tdisk,"
              " IFNULL(e_temp,-1) as e_tdisk") + sel

    # limit table sizes
    if sample != 'everything':
        selall += " LIMIT "+str(cfg.www['plotmax'])+";"

    # number we could have plotted if we knew their distances
    selnum = "SELECT COUNT(*)" + sel1
    cursor.execute(selnum)
    ntot = cursor.fetchall()[0][0]

    cursor.execute(selall)
    t = {}
    allsql = cursor.fetchall()
    ngot = len(allsql)

    if ngot == 0:
        print("    HR + f vs. r: nothing to plot")
        p = gridplot([[figure(),figure()]],sizing_mode='scale_width',toolbar_location='above')
        return bokeh.embed.components(p)
    else:
        print("    got ",ngot," rows for plot")

    # organise these into a dict of ndarrays
    l = list(zip(*allsql))
    keys = cursor.column_names
    dtypes = [None,None,float,float,float,float,
              float,float,float,float]
    for i in range(len(keys)):
        col = np.array(l[i],dtype=dtypes[i])
        t[keys[i]] = col

    # colour scale
    col,cr = colours_for_list(t['ldisklstar'],bokeh.palettes.plasma,log=True)
    _,tr = colours_for_list(t['tdisk'],bokeh.palettes.plasma,log=True)

    t['col'] = col
    n_unique = len(np.unique(t['sdbid']))
    data = ColumnDataSource(data=t)

    # set up hover/tap tools
    # TODO: hover in one highlights in the other
    hover1 = HoverTool(names=['dot'],tooltips=[("name","@main_id")])
    hover2 = HoverTool(names=['dot'],tooltips=[("name","@main_id")])
    tap1 = TapTool(names=['dot'])
    tap2 = TapTool(names=['dot'])
    tools1 = ['wheel_zoom,box_zoom,box_select,save,reset',hover1,tap1]
    tools2 = ['wheel_zoom,box_zoom,box_select,save,reset',hover2,tap2]

    # hr diagram
    hr = figure(title="HR diagram ("+str(n_unique)+" of "+str(ntot)+")",
                tools=tools1,active_scroll='wheel_zoom',
                x_axis_label='Effective temperature / K',y_axis_label='Stellar luminosity / Solar',
                y_axis_type="log",y_range=(0.5*min(t['lstar']),max(t['lstar'])*2),
                x_axis_type="log",x_range=(1.1*max(t['teff']),min(t['teff'])/1.1),
                width=750,height=800)
    xs,err_xs,ys,err_ys = utils.plot_err(t['teff'],t['e_teff'],t['lstar'],t['e_lstar'])
    hr.multi_line(xs,err_ys,line_color=t['col'],**cfg.pl['hr_e_dot'])
    hr.multi_line(err_xs,ys,line_color=t['col'],**cfg.pl['hr_e_dot'])
    hr.circle('teff','lstar',source=data,name='dot',line_color='col',
              fill_color='col',**cfg.pl['hr_dot'])

    # f vs temp (if we have any)
    if np.max(t['ldisklstar']) > 0:
        ft = figure(title="fractional luminosity vs disk temperature",
                    tools=tools2,active_scroll='wheel_zoom',
                    x_axis_label='Disk temperature / K',
                    y_axis_label='Disk fractional luminosity',
                    y_axis_type="log",y_range=(0.5*cr[0],2*cr[1]),
                    x_axis_type="log",x_range=(0.5*tr[0],2*tr[1]),
                    width=750,height=800)
                    
        xs,ys = utils.plot_join_line(t,'sdbid','tdisk','ldisklstar')
        if xs is not None:
            ft.multi_line(xs,ys,**cfg.pl['fvsr_join'])
        
        xs,err_xs,ys,err_ys = utils.plot_err(t['tdisk'],t['e_tdisk'],
                                             t['ldisklstar'], t['e_ldisklstar'])
        ft.multi_line(xs,err_ys,line_color=t['col'],**cfg.pl['hr_e_dot'])
        ft.multi_line(err_xs,ys,line_color=t['col'],**cfg.pl['hr_e_dot'])
        ft.circle('tdisk','ldisklstar',source=data,name='dot',fill_color='col',
                  line_color='col',**cfg.pl['hr_dot'])
    else:
        ft = figure(title='no IR excesses')
            
    p = gridplot([[hr,ft]],sizing_mode='scale_width',toolbar_location='above')

    # taptool callback
    if absolute_paths:
        url = "/sdb/seds/masters/@sdbid/public"
    else:
        url = "@sdbid.html"
    taptool = hr.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    taptool = ft.select(type=TapTool)
    taptool.callback = OpenURL(url=url)

    return bokeh.embed.components(p)


def flux_size_plots():

    # set up connection
    cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                  password=cfg.mysql['passwd'],
                                  host=cfg.mysql['host'],
                                  database=cfg.mysql['db_sdb'])
    cursor = cnx.cursor(buffered=True)

    env = jinja2.Environment(autoescape=False,
         loader=jinja2.PackageLoader('sdf',package_path='www/templates'))
    template = env.get_template("sample_plot.html")

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
        out = flux_size_plot(cursor,sample)
        if out is not None:
            script,div = out
        else:
            continue        

        html = template.render(
                   js=[bokeh_js],
                   css=[bokeh_css],
                   plot_script=script,
                   plot_div=div,
                   title=sample,
                   creation_time=datetime.utcnow().strftime("%d/%m/%y %X")
                               )

        with io.open(file, mode='w', encoding='utf-8') as f:
            f.write(html)
            
    cursor.close()
    cnx.close()


def flux_size_plot(cursor,sample):
    """Return bokeh components for flux vs. size sample plots.

    Parameters
    ----------
    cursor : mysql.connector.connect.cursor
        Connection to database
    sample : string
        Name of sample in config.mysql['db_samples'] to plot
    """
    # TODO: multiple compoments are plotted at the total disk flux =bad

    # bands to show disk fluxes at
    filters = ['MIRI.F1000W','MIRI.F1800W','WISE22','MIRI.F2550W',
               'PACS70','WAV880']

    # sensitivities from Smith & Wyatt 2010
    miri_10_sens = np.array([[1000.,0.15],[100,0.15],[1,0.2],[0.01,0.3],
                             [0.001,0.7],[0.001,1.5],[0.01,10],[0.1,50]])
    miri_18_sens = np.array([[1000,0.25],[100.,0.3],[10,0.4],[1,0.7],
                             [0.5,1],[0.1,2],[0.05,4],[0.05,10],[0.05,40]])
    miri_25_sens = np.array([[1000.,0.3],[100,0.3],[10,0.4],[3,0.5],[2,1],
                             [1,1.5],[0.5,2],[0.2,4],[0.2,10],[0.2,40]])
    sens = [miri_10_sens,miri_18_sens,miri_25_sens,miri_25_sens,
            None,None]

    tabs = []
    for i,f in enumerate(filters):

        # get the results
        if sample == 'everything' or sample == 'public':
            sel1 = "FROM "+cfg.mysql['db_sdb']+".sdb_pm "
        else:
            sel1 = "FROM "+cfg.mysql['db_samples']+"."+sample+" "
         
        sel = sel1 + ("LEFT JOIN "+cfg.mysql['db_results']+".model ON sdbid=id "
                      "LEFT JOIN "+cfg.mysql['db_results']+".disk_r USING (id) "
                      "LEFT JOIN "+cfg.mysql['db_results']+".star USING (id) "
                      "LEFT JOIN "+cfg.mysql['db_results']+".phot "
                      "ON (model.id=phot.id AND disk_r.disk_r_comp_no=phot.comp_no) "
                      "LEFT JOIN "+cfg.mysql['db_sdb']+".simbad USING (sdbid) "
                      "WHERE filter = %s AND phot.id IS NOT NULL "
                      "AND model_jy IS NOT NULL AND rdisk_bb IS NOT NULL "
                      "AND plx_arcsec IS NOT NULL ")
        # limit table sizes
        if sample != 'everything':
            sel += " LIMIT "+str(cfg.www['tablemax'])+";"

        selall = ("SELECT coalesce(main_id,sdbid) as id,sdbid,chisq, "
                "1e3*model_jy as flux,rdisk_bb*plx_arcsec as rdisk ") + sel

        # number in sample
        selnum = "SELECT COUNT(*)" + sel1
        cursor.execute(selnum)
        ntot = cursor.fetchall()[0][0]

        # fill a table with fluxes
        cursor.execute(selall,(str(f),))
        allsql = cursor.fetchall()
        l = list(zip(*allsql))

        if len(l) == 0:
            print("    flux vs r: nothing to plot")
            pl = figure()
            tabs.append( Panel(child=pl, title=f) )
            continue

        keys = cursor.column_names
        dtypes = [None,None,float,float,float]
        t = {}
        for j in range(len(keys)):
            col = np.array(l[j],dtype=dtypes[j])
            t[keys[j]] = col

        ngot = len(t['id'])
        print("    got ",ngot," rows for filter ",f)

        # colour scale
        col,cr = colours_for_list(t['chisq'],bokeh.palettes.plasma,log=True)
        t['col'] = col

        data = ColumnDataSource(data=t)

        hover = HoverTool(names=['dot'],tooltips=[("name","@id")])
        tools = ['wheel_zoom,box_zoom,box_select,tap,save,reset',hover]

        pl = figure(title="disk flux vs radius ("+str(ngot)+" of "+str(ntot)+")",
                    tools=tools,active_scroll='wheel_zoom',toolbar_location='above',
                    x_axis_label='Disk black body radius / arcsec',
                    y_axis_label='Disk flux / mJy',
                    y_axis_type="log",y_range=(0.5*min(t['flux']),max(t['flux'])*2),
                    x_axis_type="log",x_range=(0.5*min(t['rdisk']),max(t['rdisk'])*2),
                    width=850,height=600)

        if sens[i] is not None:
            pl.line(sens[i][:,1],sens[i][:,0],
                    line_color='red',line_alpha=0.3,line_width=4)

        xs,ys = utils.plot_join_line(t,'sdbid','rdisk','flux')
        if xs is not None:
            pl.multi_line(xs,ys,**cfg.pl['fvsr_join'])
        
        pl.circle('rdisk','flux',source=data,size=10,fill_color='col',
                  name='dot',fill_alpha=0.6,line_color='col',line_alpha=1)

        url = "/sdb/seds/masters/@sdbid/public"
        taptool = pl.select(type=TapTool)
        taptool.callback = OpenURL(url=url)

        tabs.append( Panel(child=pl, title=f) )
    
    tab = Tabs(tabs=tabs)
    return bokeh.embed.components(tab)


def colours_for_list(values_in,palette,log=False):
    """Return plotting colours for an array of values."""
    
    values = np.array(values_in)
    nval = len(values)
    col = np.repeat('#969696',nval)

    # if there's only one value
    if len(np.unique(values)) == 1:
        return col,np.array([0.5*values[0],2*values[0]])

    # figure what's OK
    ok = np.isfinite(values)
    if log:
        ok1 = values[ok] <= 0
        ok[ok1] = False

    # set ranges
    if np.sum(ok) > 0:
        if log:
            range = np.array([np.min(np.log(values[ok])),
                              np.max(np.log(values[ok]))])
            if range[0] == range[1]:
                range = np.array([range[0]-0.3,range[0]+0.3])
            ci = 0.999*(np.log(values[ok])-range[0])/(range[1]-range[0])
            range = np.exp(range)
        else:
            range = np.array([np.min(values[ok]),np.max(values[ok])])
            if range[0] == range[1]:
                range = np.array([range[0]/2.,range[0]*2])
            ci = 0.999*(values[ok]-range[0])/(range[1]-range[0])
    
        # assign colours, don't use the whole range since the top end
        # results in very whiteish (invisible) symbols
        col[ok] = np.array(palette(100))[np.floor(80*ci).astype(int)]

    return col,range


def calibration(sample='zpo_cal_',
                fileroot=cfg.file['www_root']+'calibration/'):
    """Diagnostic plot showing quality of photometric calibration.

    Parameters
    ----------
    sample : string
        Name of sample in config.mysql['db_samples'] to plot.
    fileroot : string, optional
        Location of directory in which to place plots.
    """

    try:
        cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                      password=cfg.mysql['passwd'],
                                      host=cfg.mysql['host'],
                                      database=cfg.mysql['db_results'])
        cursor = cnx.cursor(buffered=True)
        print("  "+sample)
            
    except mysql.connector.InterfaceError:
        print("Can't connect to {} at {}".format(cfg.mysql['db_results'],
                                                 cfg.mysql['host']) )
        return

    # get a wavelength-sorted list of filters.
    cursor.execute("SELECT DISTINCT filter FROM "+cfg.mysql['phot_table']+" "
                   "WHERE obs_jy IS NOT NULL LIMIT 100")
    filters = cursor.fetchall()
    filters = [f for (f,) in filters]
    wav = filter.mean_wavelength(filters)

    # use fields to allow sorting by wave and then filter name
    f_w = np.array([i for i in zip(filters,wav)],
                   dtype=[('filter','S20'),('wav',float)])
    srt = np.argsort(f_w,order=('wav','filter'))

    filters = np.array(filters)[srt]
    wav = np.array(wav)[srt]

    # now loop to make the plots
    flux = []
    rhist = []
    chist = []
    for i,f in enumerate(filters):

        if wav[i] > 26 or wav[i] < 0.017:
            continue

        # grab the data for this filter
        data = {'name':[],'sdbid':[],'chi':[],'R':[],'Teff':[]}
        col = np.array([])
        stmt = ("SELECT name, sdbid, chi_star,"
                "IF(R_star BETWEEN -100 and 100,R_star,100),"
                "parameters, chisq/IF(dof<1,1,dof) as cdof FROM "
                +cfg.mysql['model_table']+" ""LEFT JOIN "
                +cfg.mysql['phot_table']+" USING (id) "
                "LEFT JOIN "+cfg.mysql['db_samples']+'.'+sample+" ON id=sdbid "
                "WHERE sdbid IS NOT NULL AND filter='"+f+"' "
                "AND obs_upperlim=0 and chi_star != 0")
        cursor.execute(stmt)
        for (name,sdbid,chi,R,par,cdof) in cursor.fetchall():
            data['name'].append(name)
            data['sdbid'].append(sdbid)
            data['chi'].append(chi)
            data['R'].append(R)
            pars = ast.literal_eval(par)
            data['Teff'].append(pars[0])
            col = np.append(col, cdof )
        
#        print("  ",f,":",len(col))

        if len(col) == 0:
            continue

        # set colour range, clipped at chisq/dof=10
        col = np.clip(255 * col / 10.,0,255)
        data['col'] = np.array(bokeh.palettes.Viridis256)[col.astype(int)]

        pldata = ColumnDataSource(data=data)

        # set range for flux ratios
        std = np.percentile(data['R'],[5,95])

        # flux ratio plot
        hover = HoverTool(names=['pl'],tooltips=[('id',"@name")])
        tap = TapTool(names=['pl'])
        tools = ['wheel_zoom,pan,save,reset']
        flux.append( figure(x_axis_label='Teff / K',y_axis_label=f,
                            y_range=std.tolist(),
                            tools=tools+[hover,tap],active_scroll='wheel_zoom',
                            width=1100,height=200) )

        center = (0 if filter.iscolour(f) else 1)
        flux[-1].line(x=[ np.min(data['Teff']) , np.max(data['Teff']) ],
                    y=[center,center],**cfg.pl['guide_dash'])
        flux[-1].circle('Teff','R',source=pldata,name='pl',
                       fill_color='col')

        # ratio histograms (charts.Histogram uselessly slow)
        hhist,hedges = np.histogram(data['R'],bins='auto')
        hcen = (hedges[0:-1]+hedges[1:])/2.
        rhist.append( figure(x_axis_label='R',x_range=flux[-1].y_range,
                             tools=tools,
                             active_scroll='wheel_zoom',width=200,height=200) )

        rhist[-1].vbar(x=hcen,top=hhist,width=hedges[1]-hedges[0],bottom=0,
                     line_color='#333333')

        # chi histogram
        if np.max(np.array(data['chi'])>0):
            hhist,hedges = np.histogram(data['chi'],bins='auto')
            hcen = (hedges[0:-1]+hedges[1:])/2.
            chist.append( figure(x_axis_label='chi',x_range=[-5,5],
                                 tools=tools,
                                 active_scroll='wheel_zoom',width=200,height=200) )

            chist[-1].vbar(x=hcen,top=hhist,width=hedges[1]-hedges[0],bottom=0,
                         line_color='#333333')
        else:
            chist.append( figure(width=200,height=200) )

        # taptool callback
        url = "/sdb/seds/masters/@sdbid/public"
        taptool = flux[-1].select(type=TapTool)
        taptool.callback = OpenURL(url=url)

    # link all the x ranges
    for i in range(len(flux)-1):
        flux[i].x_range = flux[-1].x_range

    pl = []
    for i in range(len(flux)):
        pl.append([flux[i],rhist[i],chist[i]])

    # TODO: future bokeh release might allow sizing_mode='scale_both'
    # with different element widths
    grid = gridplot(pl,toolbar_location='above')

    script,div = bokeh.embed.components(grid)

    # now write the html
    env = jinja2.Environment(autoescape=False,
         loader=jinja2.PackageLoader('sdf',package_path='www/templates'))
    template = env.get_template("sample_plot.html")

    bokeh_js = bokeh.resources.CDN.render_js()
    bokeh_css = bokeh.resources.CDN.render_css()

    html = template.render(
                   js=[bokeh_js],
                   css=[bokeh_css],
                   body_class='wide',
                   plot_script=script,
                   plot_div=div,
                   title=sample,
                   creation_time=datetime.utcnow().strftime("%d/%m/%y %X")
                           )

    file = fileroot + sample + '.html'
    with io.open(file, mode='w', encoding='utf-8') as f:
        f.write(html)


def filter_plot(file=cfg.file['www_root']+'filters.html'):
    """Plot all filters.
        
    Parameters
    ----------
    file : string, optional
        File to write plot to.
    """

    c_micron = u.micron.to(u.Hz,equivalencies=u.spectral())
    cols = bokeh.palettes.Category20_12

    # groups in which to split the filters
    groups = cfg.filter_plot['groups']

    hover = HoverTool(tooltips=[('',"@filter")])
    tools = ['wheel_zoom,pan,save,reset']

    pl = []
    for g in groups:

        pl.append( figure(x_axis_label='wavelength / micron',
                         y_axis_label='response',
                         tools=tools,active_scroll='wheel_zoom',
                         width=1100,height=200) )

        for i,fname in enumerate(g):
            
            f = filter.Filter.get(fname)
            
            data = {'wave': c_micron/f.nu_hz,
                    'response':f.response/np.max(f.response)}
            pldata = ColumnDataSource(data=data)

            pl[-1].line('wave','response',
                        line_color=cols[i],line_width=2,
                        source=pldata,
                        legend=fname,muted_alpha=0.2)

            pl[-1].legend.click_policy="mute"
            pl[-1].legend.label_text_font_size = '8pt'

    grid = gridplot(pl,ncols=1,sizing_mode='scale_width',
                    toolbar_location='above')

    script,div = bokeh.embed.components(grid)

    # now write the html
    env = jinja2.Environment(autoescape=False,
         loader=jinja2.PackageLoader('sdf',package_path='www/templates'))
    template = env.get_template("sample_plot.html")

    bokeh_js = bokeh.resources.CDN.render_js()
    bokeh_css = bokeh.resources.CDN.render_css()

    html = template.render(
                   js=[bokeh_js],
                   css=[bokeh_css],
                   body_class='wide',
                   plot_script=script,
                   plot_div=div,
                   title='filters',
                   creation_time=datetime.utcnow().strftime("%d/%m/%y %X")
                           )

    with io.open(file, mode='w', encoding='utf-8') as f:
        f.write(html)
