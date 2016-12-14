#!/Users/grant/.virtualenvs/astropy-dev/bin/python3
#!/usr/bin/env python3

"""Generate HTML pages to browse database.

Uses non-standard version of astropy to ensure html anchors in
jsviewer tables. This is a virtualenv created into which the modified
version of astropy is installed.

"""

from os.path import isdir,isfile
from os import mkdir,remove,write
import argparse

import numpy as np
from astropy.table import Table,jsviewer
import mysql.connector

from bokeh.plotting import figure,output_file,save,ColumnDataSource
import bokeh.palettes
from bokeh.models.widgets import Panel,Tabs
from bokeh.models import HoverTool,OpenURL,TapTool
from bokeh.layouts import gridplot,layout

from sdf import config as cfg


def get_samples():
    """Get a list of samples.
    
    Get a list of samples from the database. Add "everything" and
    "public" samples, "public" might not show everything in the list,
    but "everything" will (but may not be visible to anyone).
    """
    
    cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                  password=cfg.mysql['passwd'],
                                  host=cfg.mysql['host'],
                                  database=cfg.mysql['db_results'])
    cursor = cnx.cursor(buffered=True)
    cursor.execute("SELECT DISTINCT project FROM "
                   +cfg.mysql['db_sdb']+".projects;")
    samples = cursor.fetchall() # a list of tuples
    samples = [i[0] for i in samples]
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
        sample_table(cursor,sample)

    cursor.close()
    cnx.close()


def sample_table(cursor,sample):
    """Generate an HTML page with a sample table.

    Extract the necessary information from the database and create HTML
    pages with the desired tables, one for each sample. These are 
    generated using astropy's HTML table writer and the jsviewer,which
    makes tables that are searchable and sortable.
    """

    wwwroot = cfg.www['root']+'samples/'

    # create dir and .htaccess if neeeded
    create_dir(wwwroot,sample)

    # grab table we want to display
    cursor.execute("DROP TABLE IF EXISTS hd;")
    cursor.execute("CREATE TEMPORARY TABLE hd SELECT sdbid,GROUP_CONCAT(xid) as xid"
                   " FROM sdb_pm LEFT JOIN xids USING (sdbid) WHERE xid REGEXP('^HD')"
                   " GROUP BY sdbid;")
    cursor.execute("DROP TABLE IF EXISTS hip;")
    cursor.execute("CREATE TEMPORARY TABLE hip SELECT sdbid,GROUP_CONCAT(xid) as xid"
                   " FROM sdb_pm LEFT JOIN xids USING (sdbid) WHERE xid REGEXP('^HIP')"
                   " GROUP BY sdbid;")
    cursor.execute("DROP TABLE IF EXISTS gj;")
    cursor.execute("CREATE TEMPORARY TABLE gj SELECT sdbid,GROUP_CONCAT(xid) as xid"
                   " FROM sdb_pm LEFT JOIN xids USING (sdbid) WHERE xid REGEXP('^GJ')"
                   " GROUP BY sdbid;")
    cursor.execute("DROP TABLE IF EXISTS phot;")
    cursor.execute("CREATE TEMPORARY TABLE phot SELECT"
                   " id as sdbid,ROUND(-2.5*log10(ANY_VALUE(model_jy)/3882.37),1) as Vmag"
                   " FROM sdb_results.phot WHERE filter='VJ' GROUP BY id;")
    sel = ("SELECT "
           "CONCAT('<a target=\"_blank\" href=\"../../seds/masters/',sdbid,'/public/',sdbid,'-sed.html\">',sdbid,'</a>') as sdbid,"
           "CONCAT('<a href=\"http://simbad.u-strasbg.fr/simbad/sim-basic?submit=SIMBAD+search&Ident=',main_id,'\" target=\"_blank\">',main_id,'</a>') as Simbad,"
           "hd.xid as HD,"
           "hip.xid as HIP,"
           "gj.xid as GJ,"
           "concat('<a href=\"http://irsa.ipac.caltech.edu/applications/finderchart/#id=Hydra_finderchart_finder_chart&RequestClass=ServerRequest&DoSearch=true&subsize=0.08333333400000001&thumbnail_size=medium&sources=DSS,SDSS,twomass,WISE,IRIS&overlay_catalog=true&catalog_by_radius=true&iras_radius=240&sdss_radius=5&twomass_radius=5&wise_radius=5&one_to_one=_none_&dss_bands=poss1_blue,poss1_red,poss2ukstu_blue,poss2ukstu_red,poss2ukstu_ir&SDSS_bands=u,g,r,i,z&twomass_bands=j,h,k&wise_bands=1,2,3,4&UserTargetWorldPt=',raj2000,';',dej2000,';EQ_J2000&projectId=finderchart&searchName=finder_chart&shortDesc=Finder%20Chart&isBookmarkAble=true&isDrillDownRoot=true&isSearchResult=true\" target=\"_blank\">images</a>') as Finder,"
           "Vmag,"
           "raj2000 as RA,"
           "dej2000 as `Dec`,"
           "sp_type as SpType,"
           "teff as Teff,"
           "ROUND(log10(lstar),2) as LogLstar,"
           "1e3/COALESCE(tgas.plx,simbad.plx_value) AS Dist,"
           "ROUND(log10(ldisk_lstar),1) as Log_f,"
           "tdisk as T_disk")
        
    # here we decide which samples get all targets, for now "everything" and "public"
    # get everything, but this could be changed so that "public" is some subset of
    # "everything"
    if sample == 'everything' or sample == 'public':
        sel += " FROM sdb_pm"
    else:
        sel += " FROM "+cfg.mysql['db_samples']+"."+sample+" LEFT JOIN sdb_pm USING (sdbid)"
        
    sel += (" LEFT JOIN simbad USING (sdbid)"
            " LEFT JOIN tyc2 USING (sdbid)"
            " LEFT JOIN photometry.tgas ON COALESCE(-tyc2.hip,tyc2.tyc2id)=tgas.tyc2hip"
            " LEFT JOIN sdb_results.star on sdbid=star.id"
            " LEFT JOIN sdb_results.disk_r on sdbid=disk_r.id"
            " LEFT JOIN hd USING (sdbid)"
            " LEFT JOIN hip USING (sdbid)"
            " LEFT JOIN gj USING (sdbid)"
            " LEFT JOIN phot USING (sdbid)"
            " ORDER by RA")
    # limit table sizes
    if sample != 'everything':
        sel += " LIMIT "+str(cfg.www['tablemax'])+";"

#        print(sel)
    cursor.execute(sel)
    tsamp = Table(names=cursor.column_names,
                  dtype=('S200','S200','S50','S50','S50','S1000',
                         'f','f','f','S10','f','f','f','f','f'))
    for row in cursor:
        tsamp.add_row(row)
    print("    got ",len(tsamp)," rows")

    # write html page with interactive table, astropy 1.2.1 doesn't allow all of the
    # the htmldict contents to be passed to write_table_jsviewer so links are
    # bleached out. can use jsviewer with my modifications to that function...
    fd = open(wwwroot+sample+'/table.html','w')
#        tsamp.write(fd,format='ascii.html',htmldict={'raw_html_cols':['sdbid','Simbad','Finder'],
#                                                     'raw_html_clean_kwargs':{'attributes':{'a':['href','target']}} })
    jsviewer.write_table_jsviewer(tsamp,fd,max_lines=10000,table_id=sample,
                                  table_class="display compact",
                                  jskwargs={'display_length':25},
                                  raw_html_cols=['sdbid','Simbad','Finder'],
                                  raw_html_clean_kwargs={'attributes':{'a':['href','target']
                                                            }} )
    fd.close()


def sample_plots():

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
        sample_plot(cursor,sample)
            
    cursor.close()
    cnx.close()


def sample_plot(cursor,sample):
    """Generate HTML sample plot.

    Extract the necessary information from the database and plot
    using bokeh.
    """

    wwwroot = cfg.www['root']+'samples/'

    # create dir and .htaccess if neeeded
    create_dir(wwwroot,sample)

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

    selall = ("SELECT sdbid,main_id,teff,lstar,IFNULL(ldisk_lstar,-1) as "
              "ldisklstar,IFNULL(tdisk,-1) as tdisk") + sel

    # limit table sizes
    if sample != 'everything':
        selall += " LIMIT "+str(cfg.www['tablemax'])+";"

    # number we could have plotted (more if some were nan)
    selnum = "SELECT COUNT(*)" + sel
    cursor.execute(selnum)
    ntot = cursor.fetchall()[0][0]

    cursor.execute(selall)
    t = {}
    allsql = cursor.fetchall()
    ngot = len(allsql)
    print("    got ",ngot," rows")
    l = list(zip(*allsql))
    keys = cursor.column_names
    dtypes = [None,None,float,float,float,float]
    for i in range(len(keys)):
        col = np.array(l[i],dtype=dtypes[i])
        t[keys[i]] = col

    # set up colour scale, grey for nans
    col = np.repeat('#969696',ngot)
    if 'ldisklstar' in t:
        ok = t['ldisklstar'] > 0
        if np.sum(ok) > 0:
            cr = np.array([np.nanmin(np.log(t['ldisklstar'][ok])),
                           np.nanmax(np.log(t['ldisklstar'][ok]))])
            # ensure top is below 1 for indexing
            ci = 0.999*(np.log(t['ldisklstar'])-cr[0])/(cr[1]-cr[0])
            ok = np.isfinite(ci)
            col = np.empty(ngot,dtype='U7')
            col[ok] = np.array(bokeh.palettes.plasma(100))[np.floor(90*ci[ok]).astype(int)]
            col[col==''] = '#969696'

    t['col'] = col
    data = ColumnDataSource(data=t)

    # remove the plot file to avoid overwrite warnings
    plfile = wwwroot+sample+"/hr.html"
    if isfile(plfile):
        remove(plfile)
    output_file(plfile,mode='cdn')

    # TODO: hover in one highlights in the other
    hover1 = HoverTool(tooltips=[("name","@main_id")])
    hover2 = HoverTool(tooltips=[("name","@main_id")])
    tools1 = ['wheel_zoom,box_zoom,box_select,tap,save,reset',hover1]
    tools2 = ['wheel_zoom,box_zoom,box_select,tap,save,reset',hover2]

    # hr diagram
    hr = figure(title="diagrams for "+sample+" ("+str(ngot)+" of "+str(ntot)+")",
                tools=tools1,active_scroll='wheel_zoom',
                x_axis_label='Effective temperature / K',y_axis_label='Stellar luminosity / Solar',
                y_axis_type="log",y_range=(0.5*min(t['lstar']),max(t['lstar'])*2),
                x_range=(300+max(t['teff']),min(t['teff'])-300) )
    hr.circle('teff','lstar',source=data,size=10,fill_color='col',
              fill_alpha=0.6,line_color=None)

    # f vs temp (if we have any)
    if np.max(t['ldisklstar']) > 0:
        ft = figure(tools=tools2,active_scroll='wheel_zoom',
                    x_axis_label='Disk temperature / K',y_axis_label='Disk fractional luminosity',
                    y_axis_type="log",y_range=(0.5*np.exp(cr[0]),2*np.exp(cr[1])),
                    x_axis_type="log",x_range=(0.5*min(t['tdisk']),max(t['tdisk'])*2) )
        ft.circle('tdisk','ldisklstar',source=data,size=10,fill_color='col',
                  fill_alpha=0.6,line_color=None)
    else:
        ft = figure(title='no IR excesses')
            
    p = gridplot([[hr,ft]],sizing_mode='stretch_both',
                 toolbar_location='above')
                 
    url = "/~grant/sdb/seds/masters/@sdbid/public/@sdbid"+"-sed.html"
    taptool = hr.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    taptool = ft.select(type=TapTool)
    taptool.callback = OpenURL(url=url)

    save(p)


def flux_size_plots():

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
        flux_size_plot(cursor,sample)
            
    cursor.close()
    cnx.close()


def flux_size_plot(cursor,sample):
    """Show disk fluxes at various bands vs. their size."""

    # bands to show disk fluxes at
    filters = ['2MKS','WISE3P4','AKARI9','WISE12','WISE22','PACS70']

    wwwroot = cfg.www['root']+'samples/'
    plfile = wwwroot+sample+"/fnuvsr.html"
    if isfile(plfile):
        remove(plfile)
    output_file(plfile,mode='cdn')

    tabs = []
    for f in filters:

        # get the results
        stmt = "SELECT coalesce(main_id,id),sdbid,chisq,1e3*disk_jy,rdisk*plx_arcsec "
        if sample == 'everything' or sample == 'public':
            stmt += "FROM sdb_pm "
        else:
            stmt += "FROM "+cfg.mysql['db_samples']+"."+sample+" "
         
        stmt += ("LEFT JOIN "+cfg.mysql['db_results']+".model ON sdbid=id "
                 "LEFT JOIN "+cfg.mysql['db_results']+".phot USING (id) "
                 "LEFT JOIN "+cfg.mysql['db_results']+".star USING (id) "
                 "LEFT JOIN "+cfg.mysql['db_results']+".disk_r USING (id) "
                 "LEFT JOIN "+cfg.mysql['db_sdb']+".simbad USING (sdbid) "
                 "WHERE filter = %s AND disk_jy > 0 AND id IS NOT NULL")
        # limit table sizes
        if sample != 'everything':
            stmt += " LIMIT "+str(cfg.www['tablemax'])+";"

        cursor.execute(stmt,(str(f),))

        # fill a table with fluxes
        t = {'id':[],'sdbid':[],'chisq':[],'flux':[],'rdisk':[]}
        ntot = 0
        for (id,sdbid,chisq,flux,rdisk) in cursor:
            ntot += 1
            if flux is not None:
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

        data = ColumnDataSource(data=t)

        hover = HoverTool(tooltips=[("name","@id")])
        tools = ['wheel_zoom,box_zoom,box_select,tap,save,reset',hover]

        pl = figure(title="disk flux vs radius for "+sample+" ("+str(ngot)+" of "+str(ntot)+")",
                    tools=tools,active_scroll='wheel_zoom',
                    x_axis_label='Disk black body radius / arcsec',
                    y_axis_label='Disk flux / mJy',
                    y_axis_type="log",y_range=(0.5*min(t['flux']),max(t['flux'])*2),
                    x_axis_type="log",x_range=(0.5*min(t['rdisk']),max(t['rdisk'])*2),
                    width=800,height=500)
        pl.circle('rdisk','flux',source=data,size=10,fill_color='#969696',
                  fill_alpha=0.6,line_color=None)

        url = "/~grant/sdb/seds/masters/@sdbid/public/@sdbid"+"-sed.html"
        taptool = pl.select(type=TapTool)
        taptool.callback = OpenURL(url=url)

        tabs.append( Panel(child=pl, title=f) )
    
    tab = Tabs(tabs=tabs)
    save(tab)


# run from the command line
if __name__ == "__main__":

    # inputs
    parser = argparse.ArgumentParser(description='Update web pages')
    parser.add_argument('--tables','-t',action='store_true',help='Update sample tables')
    parser.add_argument('--plots','-p',action='store_true',help='Update sample plots')
    args = parser.parse_args()

    if args.tables:
        print("Updating sample tables")
        sample_tables()

    if args.plots:
        print("Updating sample plots")
        flux_size_plots()
        sample_plots()



