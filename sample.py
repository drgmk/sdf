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
import bokeh.resources

from sdf import plotting
from sdf import templates
from sdf import utils
from sdf import config as cfg


def cleanup_sample_dirs():
    """Remove dirs for samples that no longer exist."""

    dirs = glob.glob(cfg.file['www_root']+'samples/*/')
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
    
    fs = glob.glob(cfg.file['www_root']+'calibration/*')
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
        fd.write('AuthUserFile '+cfg.file['www_root']+'.htpasswd\n')
        fd.write('AuthGroupFile '+cfg.file['www_root']+'.htgroup\n')
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

    wwwroot = cfg.file['www_root']+'samples/'

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

    wwwroot = cfg.file['www_root']+'samples/'

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
    bokeh_js = bokeh.resources.CDN.render_js()
    bokeh_css = bokeh.resources.CDN.render_css()

    # get a list of samples and generate their pages
    samples = get_samples()
    for sample in samples:
        
        print("  sample:",sample)
        
        wwwroot = cfg.file['www_root']+'samples/'

        # create dir and .htaccess if neeeded
        create_dir(wwwroot,sample)
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
    samples = get_samples()
    for sample in samples:
        
        print("  sample:",sample)
        
        wwwroot = cfg.file['www_root']+'samples/'

        # create dir and .htaccess if neeeded
        create_dir(wwwroot,sample)
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
