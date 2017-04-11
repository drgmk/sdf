import io
from datetime import datetime
from os import path,remove

import ast
import numpy as np
import matplotlib.pyplot as plt
import mysql.connector
import astropy.units as u

from jinja2 import Template
import bokeh.resources
from bokeh.plotting import figure,ColumnDataSource
from bokeh.models.widgets import Panel,Tabs
from bokeh.models import HoverTool,TapTool,OpenURL
from bokeh.layouts import gridplot,layout
import bokeh.palettes
import bokeh.embed

from . import model
from . import spectrum
from . import photometry
from . import filter
from . import fitting
from . import db
from . import templates
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
    main_extra_func : function, optional
        Function that takes a figure instance, a list of results, and
        optionally arbitrary extra keywords. Adds to the main plot.
    main_extra_kwargs : dict, optional
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

        title = "Evidence: {:.2f} | ".format(r.evidence)
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
            main_extra_func(sed[i],r,**main_extra_kwargs)

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


def hardcopy_sed(r,file='sed.pdf'):
    """Make a hardcopy SED for a specific result object."""
    # TODO: this is very simple, needs work

    fig,ax = plt.subplots(figsize=(8,6))

    # model spectra, using star and disk means multiple components will
    # be added together
    if r.star_spec is not None:
        ax.loglog(r.star_spec.wavelength,r.star_spec.fnujy)
    if r.disk_spec is not None:
        ax.loglog(r.disk_spec.wavelength,r.disk_spec.fnujy)

    # photometry
    for p in r.obs:
        if not isinstance(p,photometry.Photometry):
            continue

        ok = np.invert(p.upperlim)
        ax.errorbar(p.mean_wavelength()[ok],p.fnujy[ok],yerr=p.e_fnujy[ok],
                    fmt='o')
        ok = p.upperlim
        ax.plot(p.mean_wavelength()[ok],p.fnujy[ok],'v')

    # cosmetics
    xl,yl = sed_limits((r,))
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_xlabel('Wavelength / $\mu$m')
    ax.set_ylabel('Flux density / Jy')

    fig.savefig(file)
    plt.close(fig)


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
                if xlim[0] < xlims[0] or xlims[0] == -1:
                    xlims[0] = xlim[0]
                if xlim[1] > xlims[1] or xlims[1] == -1:
                    xlims[1] = xlim[1]
                
                ok = np.invert(o.upperlim)
                ylim = [np.min(o.fnujy[ok]),np.max(o.fnujy[ok])]
                # set lower y limit to uncertainty if negative or zero fluxes
                if ylim[0] <= 0:
                    ylim[0] = np.min(o.e_fnujy[ok])
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
        sel = " FROM sdb_pm"
    else:
        sel = (" FROM "+cfg.mysql['db_samples']+"."+sample+
               " LEFT JOIN sdb_pm USING (sdbid)")

    sel += (" LEFT JOIN simbad USING (sdbid)"
            " LEFT JOIN sdb_results.star ON sdbid=star.id"
            " LEFT JOIN sdb_results.disk_r ON sdbid=disk_r.id")

    # statement for selecting stuff to plot
    sel += " WHERE teff IS NOT NULL AND lstar IS NOT NULL"

    selall = ("SELECT sdbid,main_id,teff,e_teff,lstar,e_lstar,"
              " IFNULL(ldisk_lstar,-1) as ldisklstar,"
              " e_ldisk_lstar as e_ldisklstar,"
              " IFNULL(temp,-1) as tdisk,"
              " e_temp as e_tdisk") + sel

    # limit table sizes
    if sample != 'everything':
        selall += " LIMIT "+str(cfg.www['plotmax'])+";"

    # number we could have plotted (more if some were nan)
    selnum = "SELECT COUNT(*)" + sel
    cursor.execute(selnum)
    ntot = cursor.fetchall()[0][0]

    cursor.execute(selall)
    t = {}
    allsql = cursor.fetchall()
    ngot = len(allsql)
    print("    got ",ngot," rows for plot")
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
    hr = figure(title="HR diagram ("+str(ngot)+" of "+str(ntot)+")",
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
        url = "/~grant/sdb/seds/masters/@sdbid/public"
    else:
        url = "@sdbid.html"
    taptool = hr.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    taptool = ft.select(type=TapTool)
    taptool.callback = OpenURL(url=url)

    return bokeh.embed.components(p)


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
    filters = ['WISE3P4','AKARI9','WISE12','AKARI18','WISE22',
               'PACS70','WAV880']

    # sensitivities from Smith & Wyatt 2010
    miri_10_sens = np.array([[1000.,0.15],[100,0.15],[1,0.2],[0.01,0.3],
                             [0.001,0.7],[0.001,1.5],[0.01,10],[0.1,50]])
    miri_18_sens = np.array([[1000,0.25],[100.,0.3],[10,0.4],[1,0.7],
                             [0.5,1],[0.1,2],[0.05,4],[0.05,10],[0.05,40]])
    miri_25_sens = np.array([[1000.,0.3],[100,0.3],[10,0.4],[3,0.5],[2,1],
                             [1,1.5],[0.5,2],[0.2,4],[0.2,10],[0.2,40]])
    sens = [None,miri_10_sens,miri_10_sens,miri_18_sens,miri_25_sens,
            None,None]

    tabs = []
    for i,f in enumerate(filters):

        # get the results
        stmt = "SELECT coalesce(main_id,id),sdbid,chisq,1e3*disk_jy,rdisk_bb*plx_arcsec "
        if sample == 'everything' or sample == 'public':
            stmt += "FROM sdb_pm "
        else:
            stmt += "FROM "+cfg.mysql['db_samples']+"."+sample+" "
         
        stmt += ("LEFT JOIN "+cfg.mysql['db_results']+".model ON sdbid=id "
                 "LEFT JOIN "+cfg.mysql['db_results']+".phot USING (id) "
                 "LEFT JOIN "+cfg.mysql['db_results']+".star USING (id) "
                 "LEFT JOIN "+cfg.mysql['db_results']+".disk_r USING (id) "
                 "LEFT JOIN "+cfg.mysql['db_sdb']+".simbad USING (sdbid) "
                 "WHERE filter = %s AND id IS NOT NULL")
        # limit table sizes
        if sample != 'everything':
            stmt += " LIMIT "+str(cfg.www['tablemax'])+";"

        cursor.execute(stmt,(str(f),))

        # fill a table with fluxes
        t = {'id':[],'sdbid':[],'chisq':[],'flux':[],'rdisk':[]}
        ntot = 0
        for (id,sdbid,chisq,flux,rdisk) in cursor:
            ntot += 1
            if flux is not None and rdisk is not None:
                if flux > 0:
                    t['id'].append(id)
                    t['sdbid'].append(sdbid)
                    t['chisq'].append(chisq)
                    t['flux'].append(flux)
                    t['rdisk'].append(rdisk)

        ngot = len(t['id'])
        if ntot == 0 or ngot == 0:
            print("    flux vs r: nothing to plot")
            return
        else:
            print("    got ",ngot," rows for filter ",f)

        # colour scale
        col,cr = colours_for_list(t['chisq'],bokeh.palettes.plasma,log=True)
        t['col'] = col

        data = ColumnDataSource(data=t)

        hover = HoverTool(tooltips=[("name","@id")])
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

        pl.circle('rdisk','flux',source=data,size=10,fill_color='col',
                  fill_alpha=0.6,line_color='col',line_alpha=1)

        url = "/~grant/sdb/seds/masters/@sdbid/public"
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
            ci = 0.999*(np.log(values[ok])-range[0])/(range[1]-range[0])
            range = np.exp(range)
        else:
            range = np.array([np.min(values[ok]),np.max(values[ok])])
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
        stmt = ("SELECT name,sdbid,chi,IF(R BETWEEN -100 and 100,R,100),parameters, "
                "chisq/IF(dof<1,1,dof) as cdof FROM "
                +cfg.mysql['model_table']+" ""LEFT JOIN "
                +cfg.mysql['phot_table']+" USING (id) "
                "LEFT JOIN "+cfg.mysql['db_samples']+'.'+sample+" ON id=sdbid "
                "WHERE sdbid IS NOT NULL AND filter='"+f+"' "
                "AND obs_upperlim=0 and chi != 0")
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

        center = (1 if std[0] > 0.0 else 0)
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
        url = "/~grant/sdb/seds/masters/@sdbid/public"
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
    template = Template(templates.sample_plot_wide)
    bokeh_js = bokeh.resources.CDN.render_js()
    bokeh_css = bokeh.resources.CDN.render_css()

    html = template.render(bokeh_js=bokeh_js,
                           bokeh_css=bokeh_css,
                           css=templates.css,
                           plot_script=script,
                           plot_div=div,
                           title=sample,
                           creation_time=datetime.utcnow().strftime("%d/%m/%y %X"))

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
    cols = bokeh.palettes.Category10[10]

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
                        source=pldata,legend=fname)

            pl[-1].legend.label_text_font_size = '8pt'

    grid = gridplot(pl,ncols=1,sizing_mode='scale_width',
                    toolbar_location='above')

    script,div = bokeh.embed.components(grid)

    # now write the html
    template = Template(templates.sample_plot_wide)
    bokeh_js = bokeh.resources.CDN.render_js()
    bokeh_css = bokeh.resources.CDN.render_css()

    html = template.render(bokeh_js=bokeh_js,
                           bokeh_css=bokeh_css,
                           css=templates.css,
                           plot_script=script,
                           plot_div=div,
                           title='filters',
                           creation_time=datetime.utcnow().strftime("%d/%m/%y %X"))

    with io.open(file, mode='w', encoding='utf-8') as f:
        f.write(html)
