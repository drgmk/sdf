"""Generate HTML pages to browse database.

Uses non-standard version of astropy to ensure html anchors retained in
jsviewer tables. This is a virtualenv created into which the modified
version of astropy is installed.

"""

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
from bokeh.resources import CDN
from bokeh.plotting import figure,ColumnDataSource
import bokeh.palettes
from bokeh.models.widgets import Panel,Tabs
from bokeh.models import HoverTool,OpenURL,TapTool
from bokeh.layouts import gridplot,layout
from bokeh.embed import components

from sdf import plotting
from sdf import templates
from sdf import utils
from sdf import config as cfg


def cleanup_sample_dirs():
    """Remove dirs for samples that no longer exist."""
    # TODO: add equivalent for calibration output

    dirs = glob.glob(cfg.www['root']+'samples/*/')
    samp = get_samples()
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
    
    fs = glob.glob(cfg.www['root']+'calibration/*')
    for f in fs:
        fname = basename(f.rstrip('.html'))
        if fname not in cfg.www['cal_samples']:
            print("  {} removed ".format(f))
            remove(f)
        else:
            print("  {} ok".format(f))


def get_samples():
    """Get a list of samples.
    
    Get a list of samples from the database. Add "everything" and
    "public" samples, "public" might not show everything in the list,
    but "everything" will (but may not be visible to anyone).
    
    Samples with no sdbid column, i.e. those that haven't been imported
    yet, will be excluded.
    """
    
    cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                  password=cfg.mysql['passwd'],
                                  host=cfg.mysql['host'],
                                  database=cfg.mysql['db_samples'])
    cursor = cnx.cursor(buffered=True)
    cursor.execute("SHOW TABLES;")
    samples_tmp = cursor.fetchall() # a list of tuples
    samples = []
    for s_tuple in samples_tmp:
        s = s_tuple[0]
        cursor.execute("SHOW COLUMNS FROM {} LIKE 'sdbid'".format(s))
        if cursor.rowcount == 1:
            samples.append(s)
#    samples = [i[0] for i in samples]
    cursor.close()
    cnx.close()
    return( samples + ['public','everything'] )


def create_dir(wwwroot,sample):
    """Create sample directories and .htaccess if necessary."""

    # make dir and .htaccess if dir doesn't exist
    if not isdir(wwwroot+sample):
        mkdir(wwwroot+sample)

    # make .htaccess if needed, don't put one in "public" or those
    # ending with "_" so stuff in those directories remains visible
    # to those not logged in
    if  sample[-1] != '_' and sample != 'public':
        fd = open(wwwroot+sample+'/.htaccess','w')
        fd.write('AuthName "Must login"\n')
        fd.write('AuthType Basic\n')
        fd.write('AuthUserFile '+cfg.www['root']+'.htpasswd\n')
        fd.write('AuthGroupFile '+cfg.www['root']+'.htgroup\n')
        fd.write('require group admin '+sample+'\n')
        fd.close()


def sample_tables():
    """Generate tables for all samples."""

    # set up connection
    cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                  password=cfg.mysql['passwd'],
                                  host=cfg.mysql['host'],
                                  database=cfg.mysql['db_sdb'])
    cursor = cnx.cursor(buffered=True)

    # get a list of samples and generate their pages
    samples = get_samples()
    for sample in samples:
        print("  sample:",sample)
        sample_table_www(cursor,sample)
        sample_table_votable(cursor,sample)

    cursor.close()
    cnx.close()


def sample_table_www(cursor,sample):
    """Generate an HTML page with a sample table.

    Extract the necessary information from the database and create HTML
    pages with the desired tables, one for each sample. These are 
    generated using astropy's XMLWriter the jsviewer,which makes tables
    that are searchable and sortable.
    """

    wwwroot = cfg.www['root']+'samples/'

    # create dir and .htaccess if neeeded
    create_dir(wwwroot,sample)

    # create temporary tables we want to join on
    sample_table_temp_tables(cursor)

    sel = ("SELECT "
           "CONCAT('<a target=\"_blank\" href=\"../../seds/masters/',sdbid,'/public\">',COALESCE(main_id,hd.xid,hip.xid,gj.xid,tmass.xid),'</a>') as id,"
           "hd.xid as HD,"
           "hip.xid as HIP,"
           "gj.xid as GJ,"
           "ROUND(Vmag,1) as Vmag,"
           "ROUND(raj2000/15.,1) as `RA/h`,"
           "ROUND(dej2000,1) as `Dec`,"
           "sp_type as SpType,"
           "ROUND(teff,0) as Teff,"
           "ROUND(log10(lstar),2) as `LogL*`,"
           "ROUND(1/COALESCE(star.plx_arcsec),1) AS Dist,"
           "ROUND(log10(SUM(ldisk_lstar)),1) as Log_f,"
           "GROUP_CONCAT(ROUND(temp,1)) as T_disk")
        
    # here we decide which samples get all targets, for now "everything" and "public"
    # get everything, but this could be changed so that "public" is some subset of
    # "everything"
    if sample == 'everything' or sample == 'public':
        sel += " FROM sdb_pm"
    else:
        sel += " FROM "+cfg.mysql['db_samples']+"."+sample+" LEFT JOIN sdb_pm USING (sdbid)"
        
    sel += (" LEFT JOIN simbad USING (sdbid)"
            " LEFT JOIN sdb_results.star on sdbid=star.id"
            " LEFT JOIN sdb_results.disk_r on sdbid=disk_r.id"
            " LEFT JOIN hd USING (sdbid)"
            " LEFT JOIN hip USING (sdbid)"
            " LEFT JOIN gj USING (sdbid)"
            " LEFT JOIN tmass USING (sdbid)"
            " LEFT JOIN phot USING (sdbid)"
            " WHERE sdb_pm.sdbid IS NOT NULL"
            " GROUP BY sdbid"
            " ORDER by raj2000")
    # limit table sizes
    if sample != 'everything':
        sel += " LIMIT "+str(cfg.www['tablemax'])+";"

    cursor.execute(sel)
    tsamp = Table(names=cursor.column_names,
                  dtype=('S1000','S50','S50','S50',
                         'S4','S8','S8','S10','S5','S4','S6','S4','S12'))
    for row in cursor:
        tsamp.add_row(row)

    for n in tsamp.colnames:
        none = np.where(tsamp[n] == b'None')
        tsamp[n][none] = '-'

    print("    got ",len(tsamp)," rows for html table")

    # get the table as xml
    s = io.StringIO()
    w = xml.writer.XMLWriter(s)
    with w.xml_cleaning_method('bleach_clean'):
        with w.tag('table',attrib={'class':'display compact','id':sample}):
            with w.tag('thead',attrib={'class':'datatable_header'}):
                with w.tag('tr'):
                    for col in tsamp.colnames:
                        w.element('td',text=col)
            for i in range(len(tsamp)):
                with w.tag('tr'):
                    for txt in tsamp[i]:
                        w.element('td',text=txt.decode())

    # write the table out to html
    template = Template(templates.datatable)
    html = template.render(css=templates.css,name=sample,table=s.getvalue(),
                           creation_time=datetime.utcnow().strftime("%d/%m/%y %X"))

    with io.open(wwwroot+sample+'/index.html',
                 mode='w', encoding='utf-8') as f:
        f.write(html)


def sample_table_votable(cursor,sample):
    """Generate a votable of the results."""

    wwwroot = cfg.www['root']+'samples/'

    # create dir and .htaccess if neeeded
    create_dir(wwwroot,sample)

    # generate the mysql statement
    sel = "SELECT *"

    if sample == 'everything' or sample == 'public':
        sel += " FROM sdb_pm"
    else:
        sel += " FROM "+cfg.mysql['db_samples']+"."+sample+" LEFT JOIN sdb_pm USING (sdbid)"
        
    sel += (" LEFT JOIN simbad USING (sdbid)"
            " LEFT JOIN sdb_results.star on sdbid=star.id"
            " LEFT JOIN sdb_results.disk_r using (id)"
            " WHERE sdb_pm.sdbid IS NOT NULL"
            " ORDER by raj2000")
    # limit table sizes
    if sample != 'everything':
        sel += " LIMIT "+str(cfg.www['votmax'])+";"

    cursor.execute(sel)
    rows = cursor.fetchall()
    tsamp = Table(rows=rows,names=cursor.column_names)

    print("    got ",len(tsamp)," rows for votable")

    tsamp.write(wwwroot+sample+'/'+sample+'.xml',
                format='votable',overwrite=True)


def sample_table_temp_tables(cursor):
    """Create temporary tables for creating sample tables."""

    cursor.execute("DROP TABLE IF EXISTS tmass;")
    cursor.execute("CREATE TEMPORARY TABLE tmass SELECT sdbid,GROUP_CONCAT(xid) as xid"
                   " FROM sdb_pm LEFT JOIN xids USING (sdbid) WHERE xid REGEXP('^2MASS')"
                   " GROUP BY sdbid;")
    cursor.execute("ALTER TABLE tmass ADD INDEX sdbid_tmass (sdbid);")
    cursor.execute("DROP TABLE IF EXISTS hd;")
    cursor.execute("CREATE TEMPORARY TABLE hd SELECT sdbid,GROUP_CONCAT(xid) as xid"
                   " FROM sdb_pm LEFT JOIN xids USING (sdbid) WHERE xid REGEXP('^HD')"
                   " GROUP BY sdbid;")
    cursor.execute("ALTER TABLE hd ADD INDEX sdbid_hd (sdbid);")
    cursor.execute("DROP TABLE IF EXISTS hip;")
    cursor.execute("CREATE TEMPORARY TABLE hip SELECT sdbid,GROUP_CONCAT(xid) as xid"
                   " FROM sdb_pm LEFT JOIN xids USING (sdbid) WHERE xid REGEXP('^HIP')"
                   " GROUP BY sdbid;")
    cursor.execute("ALTER TABLE hip ADD INDEX sdbid_hip (sdbid);")
    cursor.execute("DROP TABLE IF EXISTS gj;")
    cursor.execute("CREATE TEMPORARY TABLE gj SELECT sdbid,GROUP_CONCAT(xid) as xid"
                   " FROM sdb_pm LEFT JOIN xids USING (sdbid) WHERE xid REGEXP('^GJ')"
                   " GROUP BY sdbid;")
    cursor.execute("ALTER TABLE gj ADD INDEX sdbid_gj (sdbid);")
    cursor.execute("DROP TABLE IF EXISTS phot;")
    cursor.execute("CREATE TEMPORARY TABLE phot SELECT"
                   " id as sdbid,ROUND(-2.5*log10(ANY_VALUE(model_jy)/3882.37),1) as Vmag"
                   " FROM sdb_results.phot WHERE filter='VJ' GROUP BY id;")
    cursor.execute("ALTER TABLE phot ADD INDEX sdbid_phot (sdbid);")


def sample_plots():

    # set up connection
    cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                  password=cfg.mysql['passwd'],
                                  host=cfg.mysql['host'],
                                  database=cfg.mysql['db_sdb'])
    cursor = cnx.cursor(buffered=True)

    template = Template(templates.sample_plot_wide)
    bokeh_js = CDN.render_js()
    bokeh_css = CDN.render_css()

    # get a list of samples and generate their pages
    samples = get_samples()
    for sample in samples:
        
        print("  sample:",sample)
        
        wwwroot = cfg.www['root']+'samples/'

        # create dir and .htaccess if neeeded
        create_dir(wwwroot,sample)
        file = wwwroot+sample+"/hr.html"
        script,div = sample_plot(cursor,sample)

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


def sample_plot(cursor,sample):
    """Generate HTML sample plot.

    Extract the necessary information from the database and plot
    using bokeh.
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
    url = "/~grant/sdb/seds/masters/@sdbid/public"
    taptool = hr.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    taptool = ft.select(type=TapTool)
    taptool.callback = OpenURL(url=url)

    return components(p)


def flux_size_plots():

    # set up connection
    cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                  password=cfg.mysql['passwd'],
                                  host=cfg.mysql['host'],
                                  database=cfg.mysql['db_sdb'])
    cursor = cnx.cursor(buffered=True)


    template = Template(templates.sample_plot)
    bokeh_js = CDN.render_js()
    bokeh_css = CDN.render_css()

    # get a list of samples and generate their pages
    samples = get_samples()
    for sample in samples:
        
        print("  sample:",sample)
        
        wwwroot = cfg.www['root']+'samples/'

        # create dir and .htaccess if neeeded
        create_dir(wwwroot,sample)
        file = wwwroot+sample+"/fnuvsr.html"
        out = flux_size_plot(cursor,sample)
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


def flux_size_plot(cursor,sample):
    """Show disk fluxes at various bands vs. their size."""

    # bands to show disk fluxes at
    filters = ['WISE3P4','AKARI9','WISE12','AKARI18','WISE22','PACS70']

    # sensitivities from Smith & Wyatt 2010
    miri_10_sens = np.array([[1000.,0.15],[100,0.15],[1,0.2],[0.01,0.3],
                             [0.001,0.7],[0.001,1.5],[0.01,10],[0.1,50]])
    miri_18_sens = np.array([[1000,0.25],[100.,0.3],[10,0.4],[1,0.7],
                             [0.5,1],[0.1,2],[0.05,4],[0.05,10],[0.05,40]])
    miri_25_sens = np.array([[1000.,0.3],[100,0.3],[10,0.4],[3,0.5],[2,1],
                             [1,1.5],[0.5,2],[0.2,4],[0.2,10],[0.2,40]])
    sens = [None,miri_10_sens,miri_10_sens,miri_18_sens,miri_25_sens,None]

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
    return components(tab)


def colours_for_list(values_in,palette,log=False):
    """Return colours for an array of values."""
    
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
