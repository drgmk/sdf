from os import path,remove

import ast
import numpy as np
import matplotlib.pyplot as plt
import mysql.connector

from bokeh.plotting import figure,output_file,save,ColumnDataSource
from bokeh.models.widgets import Panel,Tabs
from bokeh.models import HoverTool
from bokeh.layouts import gridplot
from bokeh import palettes

from . import model
from . import spectrum
from . import photometry
from . import filter
from . import fitting
from .utils import SdfError
from . import config as cfg


def sed(results,tab_order=None,file='sed.html'):
    """Plot an SED with observations and models.
        
    If multiple models and parameters are passed (i.e. the length of the
    results list is more than 1 then these are placed in tabs.
    """

    # results should be list so we can loop
    if not isinstance(results,list):
        raise SdfError("expected list of results, not {}".format(type(results)))
    
    # remove file to avoid errors
    if path.exists(file):
        remove(file)
    output_file(file,mode='inline')

    # figure plot limits
    xlim,ylim = sed_limits(results)

    # loop to create tabs in plot
    sed = []
    res = []
    tabs = []
    for i,r in enumerate(results):

        title = "Evidence: {:.2f} | ".format(r.evidence)
        for j,par in enumerate(r.parameters):
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
                           width=800,height=500) )
        sed[i].xaxis.visible=False

        # add models
        add_model_phot(sed[i],r)
        add_model_spec(sed[i],r)

        # add obs
        add_obs_phot(sed[i],r)
        add_obs_spec(sed[i],r)

        # residuals
        hover = HoverTool(names=['resid'],tooltips=[('band',"@filter"),
                                                    ('sigma',"@res")])
        tools = ['wheel_zoom,box_zoom,save,reset',hover]

        res.append( figure(x_axis_label='Wavelength / micron',x_axis_type='log',
                           x_range=sed[i].x_range,
                           y_axis_label='Significance',
                           tools=tools,active_scroll='wheel_zoom',
                           toolbar_location='right',width=800,height=130) )

        # add residuals
        res[i].line(x=xlim,y=[0 , 0],**cfg.pl['guide_dash'])
        res[i].line(x=xlim,y=[-3,-3],**cfg.pl['guide_dash'])
        res[i].line(x=xlim,y=[3 , 3],**cfg.pl['guide_dash'])
        add_res(res[i],r)

        grid = gridplot([[sed[i]],[res[i]]],toolbar_location='right')

        tabs.append( Panel(child=grid, title=r.model_info['name']) )
    
    if tab_order is not None:
        tabs = [tabs[i] for i in tab_order]
    tab = Tabs(tabs=tabs)
    save(tab)


def add_obs_phot(fig,r):
    """Add observed photometry to an SED plot
    
    If the object passed is not photometry.Photometry,
    don't do anything
    """
    
    for p in r.obs:
        if not isinstance(p,photometry.Photometry):
            continue
        data = {}
        ok = np.invert(p.upperlim)
        data['filter'] = p.filters[ok]
        data['wave'] = p.mean_wavelength()[ok]
        data['flux'] = p.fnujy[ok]
        data['err'] = p.e_fnujy[ok]
        pldata = ColumnDataSource(data=data)

        # plot detections
        fig.circle('wave','flux',source=pldata,name='phot',**cfg.pl['obs_ph'])

        # create and plot errorbars
        err_xs = []
        err_ys = []
        for x, y, yerr in zip(data['wave'], data['flux'], data['err']):
            err_xs.append((x, x))
            err_ys.append((y - yerr, y + yerr))
        fig.multi_line(err_xs,err_ys,**cfg.pl['obs_e_ph'])

        # and remaining uppper limits
        data['filter'] = p.filters[p.upperlim]
        data['wave'] = p.mean_wavelength()[p.upperlim]
        data['flux'] = p.fnujy[p.upperlim]
        data['err'] = p.e_fnujy[p.upperlim]
        pldata = ColumnDataSource(data=data)
        
        # plot non-detections
        fig.inverted_triangle('wave','flux',source=pldata,name='phot',
                              **cfg.pl['obs_ph'])


def add_obs_spec(fig,r):
    """Add observed spectra to an SED plot
        
    If the object passed is not spectrum.ObsSpectrum,
    don't do anything
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
    """Add model photometry to an SED plot."""
    
    # where filter and colours are (spectra have None for filter)
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


def add_model_spec(fig,r):
    """Add model spectra to an SED plot.
        
    If object components are not model.SpecModel,
    don't do anything
    """
    
    # create wavelengths to plot with
    wav_pl = np.power(10,np.arange(-1,4,0.01))

    i0 = 0
    i = 0
    for comp in r.pl_models:
        if not isinstance(comp,tuple):
            comp = (comp,)
        nparam = len(comp[0].parameters)+1
        for m in comp:
            if not isinstance(m,model.SpecModel):
                continue
            data = {}
            data['wave'] = wav_pl
            m_pl = m.copy()
            m_pl.interp_to_wavelengths(wav_pl)
            data['flux'] = m_pl.fnujy(r.best_params[i0:i0+nparam])
            try:
                totflux += data['flux']
            except NameError:
                totflux = np.zeros(len(data['flux']))
                totflux += data['flux']
            pldata = ColumnDataSource(data=data)
            
            fig.line('wave','flux',source=pldata,**cfg.pl['mod_sp'][i+1])

            i += 1
        i0 += nparam

    # sum of all models
    if i > 1 and i0 > 0:
        data = {'wave':wav_pl,'flux':totflux}
        pldata = ColumnDataSource(data=data)
        fig.line('wave','flux',source=pldata,**cfg.pl['mod_sp'][0])


def add_res(fig,r):
    """Add residuals to a plot."""

#    res,wav,filt = fitting.residual(param,o,m)
    res = r.residuals
    wav = r.wavelengths
    filt = r.filters
    data = {}
    data['filter'] = filt
    data['wave'] = wav
    data['res'] = res
    pldata = ColumnDataSource(data=data)
    fig.circle('wave','res',source=pldata,name='resid',**cfg.pl['obs_ph'])


def sed_limits(results):
    """Figure out plotting limits given observations
        
    For now just using photometry.
    """

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

    xlims = [xlims[0]/pad,  xlims[1]*pad]
    ylims = [ylims[0]/padlg,ylims[1]*pad]
    return xlims,ylims


def calibration(file):
    """Diagnostic plot showing quality of photometric calibration."""

    try:
        cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                      password=cfg.mysql['passwd'],
                                      host=cfg.mysql['host'],
                                      database=cfg.mysql['db_results'])
        cursor = cnx.cursor(buffered=True)
        print("Photometric calibration")
            
    except mysql.connector.InterfaceError:
        print("Can't connect to {} at {}".format(cfg.mysql['db_results'],
                                                 cfg.mysql['host']) )
        return

    output_file(file,mode='cdn')

    # get a wavelength-sorted list of filters. TODO, use 'order' to
    # also sort by filter name
    cursor.execute("SELECT DISTINCT filter FROM "+cfg.mysql['phot_table'])
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

        if wav[i] > 4 or wav[i] < 0.017:
            continue

        # grab the data for this filter
        data = {'id':[],'chi':[],'R':[],'Teff':[]}
        col = np.array([])
        stmt = ("SELECT id,chi,IF(R BETWEEN -100 and 100,R,100),parameters, "
                "chisq/IF(dof<1,1,dof) as cdof FROM "
                +cfg.mysql['model_table']+" ""LEFT JOIN "
                +cfg.mysql['phot_table']+" USING (id) "
                "WHERE filter='"+f+"' AND obs_upperlim=0")
        cursor.execute(stmt)
        for (id,chi,R,par,cdof) in cursor.fetchall():
            data['id'].append(id)
            data['chi'].append(chi)
            data['R'].append(R)
            pars = ast.literal_eval(par)
            data['Teff'].append(pars[0])
            col = np.append(col, cdof )
        
        print("  ",f,":",len(col))

        # set colour range, clipped at chisq/dof=10
        col = np.clip(255 * col / 10.,0,255)
        data['col'] = np.array(palettes.Viridis256)[col.astype(int)]

        pldata = ColumnDataSource(data=data)

        # set range for flux ratios
        std = np.percentile(data['R'],[5,95])

        # flux ratio plot
        hover = HoverTool(names=['pl'],tooltips=[('id',"@id")])
        tools = ['wheel_zoom,pan,save,reset']
        flux.append( figure(x_axis_label='Teff / K',y_axis_label=f,
                            y_range=std.tolist(),
                            tools=tools+[hover],active_scroll='wheel_zoom',
                            width=800,height=200) )

        center = (1 if std[0] > 0.5 else 0)
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

    # link all the x ranges
    for i in range(len(flux)-1):
        flux[i].x_range = flux[-1].x_range

    # add a title
    flux[0].title.text = 'Photometric calibration'

    pl = []
    for i in range(len(flux)):
        pl.append([flux[i],rhist[i],chist[i]])

    grid = gridplot(pl,toolbar_location='above')

    save(grid)
