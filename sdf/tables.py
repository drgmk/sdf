"""Generate www tables."""

import io
import ast
from datetime import datetime
from os.path import isdir,isfile,basename
from os import mkdir,remove,write,rmdir

import numpy as np
from astropy.table import Table
from astropy.utils import xml
import mysql.connector

import jinja2
import bokeh.resources

from . import db
from . import www
from . import utils
from . import config as cfg


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
        sample_table_www(cursor,sample)
        sample_table_votable(cursor,sample)

    cursor.close()
    cnx.close()


def sample_table_www(cursor,sample,file='index.html',
                     absolute_paths=True,rel_loc=None):
    """Generate an HTML page with a sample table.

    Extract the necessary information from the database and create HTML
    pages with the desired tables, one for each sample. These are 
    generated using astropy's XMLWriter the jsviewer, which makes tables
    that are searchable and sortable.
    
    Parameters
    ----------
    cursor : mysql.connector.connect.cursor
        Cursor used to execute database query
    file : str, optional
        File name of html table
    absolute_paths : bool, optional
        Use absolute paths to sdf files, set by cfg.www and structure
        set by sdb_getphot.py. Otherwise files are in relative location
        set by rel_loc keyword.
    rel_loc : str, optional
        Relative location of sdf files to the created table, used only
        if absolute_paths is False. This string goes in the middle of an
        sql CONCAT statement, so could have sql in it, in which case
        double and single-quotes must be used
        (e.g. "folder1/',sql_column,'/file.html")
    """

    wwwroot = cfg.www['site_root']
    sample_root = cfg.file['www_root']+'samples/'

    # create temporary tables we want to join on
    sample_table_temp_tables(cursor)

    # get link for sed and create dir and .htaccess if
    # needed, absolute is used for main sample pages, otherwise these
    # are for special cases
    if absolute_paths:
        www.create_dir(sample_root,sample)
        url_str = wwwroot+"seds/masters/',sdbid,'/public"
        file = sample_root+sample+'/'+file
    else:
        if rel_loc is None:
            url_str = "',sdbid,'.html"
        else:
            url_str = rel_loc

    sel = ("SELECT "
           "CONCAT('<a target=\"_blank\" href=\""+url_str+"\">',"
              "COALESCE(main_id,hd.xid,hip.xid,gj.xid,tmass.xid,sdbid),"
              "'<span><img src=\""+url_str+"/',sdbid,'_thumb.png\"></span></a>') as id,"
           "hd.xid as HD,"
           "hip.xid as HIP,"
           "gj.xid as GJ,"
           "ROUND(Vmag,1) as Vmag,"
           "ROUND(raj2000/15.,1) as `RA/h`,"
           "ROUND(dej2000,1) as `Dec`,"
           "sp_type as SpType,"
           "ROUND(sdb_results.star.teff,0) as Teff,"
           "ROUND(log10(lstar),2) as `LogL*`,"
           "ROUND(1/COALESCE(star.plx_arcsec),1) AS Dist,"
           "ROUND(log10(SUM(ldisk_lstar)),1) as Log_f,"
           "GROUP_CONCAT(ROUND(temp,1)) as T_disk")
        
    # here we decide which samples get all targets, for now "everything"
    # and "public" get everything, but this could be changed so that
    # "public" is some subset of "everything"
    if sample == 'everything' or sample == 'public':
        sel += " FROM sdb_pm"
    else:
        sel += (" FROM "+cfg.mysql['db_samples']+"."+sample+" "
                "LEFT JOIN sdb_pm USING (sdbid)")
        
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
    with w.xml_cleaning_method('bleach_clean',tags=['a','img','span'],
                               attributes=['src','href','target']):
        with w.tag('table',attrib={'class':'display compact','id':sample}):
            with w.tag('thead',attrib={'class':'datatable_header'}):
                with w.tag('tr'):
                    for col in tsamp.colnames:
                        w.element('td',text=col)
            for i in range(len(tsamp)):
                with w.tag('tr'):
                    for j,txt in enumerate(tsamp[i]):
                        if j == 0:
                            w.element('td',text=txt,
                                      attrib={'class':'td_img_hover'})
                        else:
                            w.element('td',text=txt)

    # write the table out to html
    env = jinja2.Environment(autoescape=False,
         loader=jinja2.PackageLoader('sdf',package_path='www/templates'))
    template = env.get_template("sample_table.html")

    html = template.render(title=sample,table=s.getvalue(),
                 creation_time=datetime.utcnow().strftime("%d/%m/%y %X"))

    with io.open(file, mode='w', encoding='utf-8') as f:
        f.write(html)


def sample_table_votable(cursor,sample):
    """Generate a votable of the results."""

    wwwroot = cfg.file['www_root']+'samples/'

    # create dir and .htaccess if neeeded
    www.create_dir(wwwroot,sample)

    # generate the mysql statement
    sel = "SELECT *"

    if sample == 'everything' or sample == 'public':
        sel += " FROM sdb_pm"
    else:
        sel += (" FROM "+cfg.mysql['db_samples']+"."+sample+" "
                "LEFT JOIN sdb_pm USING (sdbid)")
        
    sel += (" LEFT JOIN simbad USING (sdbid)"
            " LEFT JOIN sdb_results.star ON sdbid=star.id"
            " LEFT JOIN sdb_results.disk_r USING (id)"
            " LEFT JOIN sdb_results.model USING (id)"
            " WHERE sdb_pm.sdbid IS NOT NULL"
            " ORDER by raj2000")
    # limit table sizes
    if sample != 'everything':
        sel += " LIMIT "+str(cfg.www['votmax'])+";"

    cursor.execute(sel)
    rows = cursor.fetchall()
    tsamp = Table(rows=rows,names=cursor.column_names)

    # add some url columns with links
    tsamp['url'] = np.core.defchararray.add(
                     np.core.defchararray.add(
        np.repeat(cfg.www['site_url']+'seds/masters/',len(tsamp)),tsamp['sdbid']
                                              ),
                        np.repeat('/public/',len(tsamp))
                                             )
    # to photometry file
    tsamp['phot_url'] = np.core.defchararray.add(
                         np.core.defchararray.add(
                            tsamp['url'],tsamp['sdbid']
                                                  ),
                    np.repeat('-rawphot.txt',len(tsamp))
                                                 )
    # to best fit model json
    tsamp['model_url'] = np.empty(len(tsamp),dtype='U150')
    for i,comps in enumerate(tsamp['model_comps']):
        if comps is not None:
            try:
                tsamp['model_url'][i] = (
                     cfg.www['site_url'] + 'seds/masters/' +
                     tsamp['sdbid'][i] + '/public/' +
                     tsamp['sdbid'][i] + cfg.fitting['pmn_dir_suffix'] + '/' +
                     cfg.fitting['model_join'].join(ast.literal_eval(comps)) +
                     cfg.fitting['pmn_model_suffix'] + '.json'
                                         )
            except ValueError:
                raise utils.SdfError("{}".format(comps))

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
