"""Generate www tables."""

import io
import ast
from datetime import datetime

import numpy as np
from astropy.table import Table
from astropy.utils import xml

import jinja2

from . import db
from . import www
from . import utils
from . import config as cfg


def sample_tables(samples=None):
    """Generate tables for all samples."""

    # set up connection
    cnx, cursor = db.get_cnx(cfg.db['db_sdb'])

    # get a list of samples and generate their pages
    if samples is None:
        samples = db.get_samples()
    
    for sample in samples:
        print("  sample:", sample)
        sample_table_www(cursor, sample)
        sample_table_votable(cursor, sample)
        sample_table_photometry(cursor, sample)

    cursor.close()
    cnx.close()


def sample_table_www(cursor, sample, file='index.html',
                     absolute_paths=True, rel_loc=None):
    """Generate an HTML page with a sample table.

    Extract the necessary information from the database and create HTML
    pages with the desired tables, one for each sample. These are 
    generated using astropy's XMLWriter the jsviewer, which makes tables
    that are searchable and sortable.
    
    Change sdf.config.mysql to specify db tables to look in.

    Parameters
    ----------
    cursor : mysql.connector.connect.cursor
        Cursor used to execute database query
    sample : str
        Name of sample.
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
        (e.g. "folder1/', sql_column, '/file.html")
    """

    sample_root = cfg.file['www_root']+'samples/'

    # create temporary tables we want to join on
    sample_table_temp_tables(cursor)

    # get link for sed and create dir and .htaccess if
    # needed, absolute is used for main sample pages, otherwise these
    # are for special cases
    if absolute_paths:
        www.create_dir(sample_root, sample)
        url_str = cfg.www['sdb_path'] + "/seds/masters/' || sdbid || '/public"
        file = sample_root+sample+'/'+file
    else:
        if rel_loc is None:
            url_str = "' || sdbid || '.html"
        else:
            url_str = rel_loc

    sel = ("SELECT "
           "'<a target=\"_blank\" href=\""+url_str+"\">' || "
           "COALESCE(main_id, hd.xid, hip.xid, gj.xid, tmass.xid, sdbid) ||"
           "'<span><img src=\""+url_str+"/' || sdbid || '_thumb.png\"></span></a>' as id, "
           "hd.xid as HD, "
           "hip.xid as HIP, "
           "gj.xid as GJ, "
           "ROUND(Vmag, 1) as Vmag, "
           "ROUND(raj2000/15., 1) as `RA/h`, "
           "ROUND(dej2000, 1) as `Dec`, "
           "sp_type as SpType, "
           "ROUND(star_tmp.teff, 0) as Teff, "
           "ROUND(log10(star_tmp.lstar), 2) as `LogL*`, "
           "ROUND(1/star_tmp.plx_arcsec, 1) AS Dist, "
           "ROUND(log10(disk_r_tmp.ldisk_lstar), 1) as Log_f, "
           "disk_r_tmp.temp as T_disk")

    # here we decide which samples get all targets, for now "everything"
    # and "public" get everything, but this could be changed so that
    # "public" is some subset of "everything"
    if sample == 'everything' or sample == 'public':
        sel += " FROM sdb_pm"
    else:
        sel += (" FROM "+cfg.db['db_samples']+"."+sample+" "
                "LEFT JOIN sdb_pm USING (sdbid)")
        
    sel += (" LEFT JOIN simbad USING (sdbid)"
            " LEFT JOIN star_tmp on sdbid=star_tmp.id"
            " LEFT JOIN disk_r_tmp on sdbid=disk_r_tmp.id"
            " LEFT JOIN hd USING (sdbid)"
            " LEFT JOIN hip USING (sdbid)"
            " LEFT JOIN gj USING (sdbid)"
            " LEFT JOIN tmass USING (sdbid)"
            " LEFT JOIN phot_tmp USING (sdbid)"
            " WHERE sdb_pm.sdbid IS NOT NULL"
            " ORDER by raj2000")
    # limit table sizes
    if sample != 'everything':
        sel += " LIMIT "+str(cfg.www['tablemax'])+";"

    cursor.execute(sel)
    tsamp = Table(names=[desc[0] for desc in cursor.description],
                  dtype=('S1000', 'S50', 'S50', 'S50',
                         'S4', 'S8', 'S8', 'S10', 'S5', 'S4', 'S6', 'S4', 'S12'))
    for row in cursor:
        add = [str(x) for x in row]
        tsamp.add_row(add)

    for n in tsamp.colnames:
        none = np.where(tsamp[n] == 'None')
        tsamp[n][none] = '-'

    print("    got ", len(tsamp), " rows for html table")

    # get the table as xml
    s = io.StringIO()
    w = xml.writer.XMLWriter(s)
    with w.xml_cleaning_method('bleach_clean', tags=['a', 'img', 'span'],
                               attributes=['src', 'href', 'target']):
        with w.tag('table', attrib={'class': 'display compact', 'id': sample}):
            with w.tag('thead', attrib={'class': 'datatable_header'}):
                with w.tag('tr'):
                    for col in tsamp.colnames:
                        w.element('td', text=col)
            for i in range(len(tsamp)):
                with w.tag('tr'):
                    for j, txt in enumerate(tsamp[i]):
                        if j == 0:
                            w.element('td', text=txt,
                                      attrib={'class': 'td_img_hover'})
                        else:
                            w.element('td', text=txt)

    # write the table out to html
    env = jinja2.Environment(autoescape=False,
                             loader=jinja2.PackageLoader('sdf', package_path='www/templates'))
    template = env.get_template("sample_table.html")

    html = template.render(title=sample, table=s.getvalue(),
                           sdb_url=cfg.www['sdb_url'],
                           creation_time=datetime.utcnow().strftime("%d/%m/%y %X"))

    with io.open(file, mode='w', encoding='utf-8') as f:
        f.write(html)


def sample_table_votable(cursor, sample, file_path=None):
    """Generate a votable of the results.

    Change sdf.config.mysql to specify db tables to look in.

    Parameters
    ----------
    cursor : mysql.connector.connect.cursor
        Cursor pointing to main sdb database.
    sample : str
        Name of sample.
    file_path : str, optional
        Where to put the file, defaults set by www config.
    """

    # create dir and .htaccess if neeeded
    if file_path is None:
        wwwroot = cfg.file['www_root']+'samples/'
        www.create_dir(wwwroot, sample)
        file_path = wwwroot+sample+'/'

    # generate the mysql statement
    sel = "SELECT *"

    if sample == 'everything' or sample == 'public':
        sel += " FROM sdb_pm"
    else:
        sel += (" FROM "+cfg.db['db_samples']+"."+sample+" "
                "LEFT JOIN sdb_pm USING (sdbid)")
        
    sel += (" LEFT JOIN simbad USING (sdbid)"
            " LEFT JOIN "+cfg.db['db_results']+".star ON sdbid=star.id"
            " LEFT JOIN "+cfg.db['db_results']+".disk_r USING (id)"
            " LEFT JOIN "+cfg.db['db_results']+".model USING (id)"
            " WHERE sdb_pm.sdbid IS NOT NULL"
            " ORDER by raj2000, disk_r.temp")
    # limit table sizes
    if sample != 'everything':
        sel += " LIMIT "+str(cfg.www['votmax'])+";"

    cursor.execute(sel)
    rows = cursor.fetchall()
    tsamp = Table(rows=rows, names=[desc[0] for desc in cursor.description])

    # add some url columns with links
    tsamp['url'] = np.core.defchararray.add(
                     np.core.defchararray.add(
        np.repeat(cfg.www['sdb_url']+'/seds/masters/', len(tsamp)), tsamp['sdbid']
                                              ),
                        np.repeat('/public/', len(tsamp))
                                             )
    # to photometry file
    tsamp['phot_url'] = np.core.defchararray.add(
                         np.core.defchararray.add(
                            tsamp['url'], tsamp['sdbid']
                                                  ),
                    np.repeat('-rawphot.txt', len(tsamp))
                                                 )
    # to mnest folder
    tsamp['mnest_url'] = np.core.defchararray.add(
                          np.core.defchararray.add(
                            tsamp['url'], tsamp['sdbid']
                                                  ),
                    np.repeat(cfg.fitting['pmn_dir_suffix'], len(tsamp))
                                                 )

    # to best fit model json
    tsamp['model_url'] = np.empty(len(tsamp), dtype='U250')
    for i, comps in enumerate(tsamp['model_comps']):
        if comps is not None:
            try:
                tsamp['model_url'][i] = (
                     cfg.www['sdb_url'] + '/seds/masters/' +
                     tsamp['sdbid'][i] + '/public/' +
                     tsamp['sdbid'][i] + cfg.fitting['pmn_dir_suffix'] + '/' +
                     cfg.fitting['model_join'].join(ast.literal_eval(comps)) +
                     cfg.fitting['pmn_model_suffix'] + '.json'
                                         )
            except ValueError:
                raise utils.SdfError("{}".format(comps))

    print("    got ", len(tsamp), " rows for votable")

    # this may get written to the votable in the future...
    tsamp.meta = {'updated': datetime.utcnow().strftime("%d/%m/%y %X")}

    tsamp.write(file_path+sample+'.xml',
                format='votable', overwrite=True)


def sample_table_photometry(cursor, sample, file_path=None):
    """Generate a table of the photometry.
        
    Seems that votable can't have a blank entry where there needn't be
    one (e.g. here an observed flux when none was observed). So write
    as a csv.

    Change sdf.config.mysql to specify db tables to look in.

    Parameters
    ----------
    cursor : mysql.connector.connect.cursor
        Cursor pointing to main sdb database.
    sample : str
        Name of sample.
    file_path : str, optional
        Where to put the file, defaults set by www config.
    """

    # create dir and .htaccess if neeeded
    if file_path is None:
        wwwroot = cfg.file['www_root']+'samples/'
        www.create_dir(wwwroot, sample)
        file_path = wwwroot+sample+'/'

    # generate the mysql statement
    sel = "SELECT name, phot.*"

    if sample == 'everything' or sample == 'public':
        sel += " FROM sdb_pm"
    else:
        sel += (" FROM "+cfg.db['db_samples']+"."+sample+" "
                "LEFT JOIN "+cfg.db['db_results']+".phot ON sdbid=id"
                " WHERE comp_no=-1 ORDER BY filter")

    # these are large, so don't limit table sizes
#    if sample != 'everything':
#        sel += " LIMIT "+str(cfg.www['votmax'])+";"

    cursor.execute(sel)
    rows = cursor.fetchall()
    tsamp = Table(rows=rows, names=[desc[0] for desc in cursor.description], masked=True)

    print("    got ", len(tsamp), " rows for photometry table")

    # mask blank entries
    for n in tsamp.colnames:
        no = tsamp[n] is None
        if np.any(no):
            tsamp[n].mask = no

    tsamp.write(file_path+sample+'_photometry.csv',
                format='csv', overwrite=True)


def sample_table_temp_tables(cursor):
    """Create temporary tables for creating sample tables."""

    cursor.execute("DROP TABLE IF EXISTS tmass;")
    cursor.execute("CREATE TEMPORARY TABLE tmass AS SELECT sdbid, GROUP_CONCAT(xid) as xid"
                   " FROM sdb_pm LEFT JOIN xids USING (sdbid) WHERE xid LIKE '2MASS%'"
                   " GROUP BY sdbid;")
    cursor.execute("CREATE INDEX sdbid_tmass ON tmass (sdbid);")

    cursor.execute("DROP TABLE IF EXISTS hd;")
    cursor.execute("CREATE TEMPORARY TABLE hd AS SELECT sdbid, GROUP_CONCAT(xid) as xid"
                   " FROM sdb_pm LEFT JOIN xids USING (sdbid) WHERE xid LIKE 'HD%'"
                   " GROUP BY sdbid;")
    cursor.execute("CREATE INDEX sdbid_hd ON hd (sdbid);")

    cursor.execute("DROP TABLE IF EXISTS hip;")
    cursor.execute("CREATE TEMPORARY TABLE hip AS SELECT sdbid, GROUP_CONCAT(xid) as xid"
                   " FROM sdb_pm LEFT JOIN xids USING (sdbid) WHERE xid LIKE 'HIP%'"
                   " GROUP BY sdbid;")
    cursor.execute("CREATE INDEX sdbid_hip ON hip (sdbid);")

    cursor.execute("DROP TABLE IF EXISTS gj;")
    cursor.execute("CREATE TEMPORARY TABLE gj AS SELECT sdbid, GROUP_CONCAT(xid) as xid"
                   " FROM sdb_pm LEFT JOIN xids USING (sdbid) WHERE xid LIKE 'GJ%'"
                   " GROUP BY sdbid;")
    cursor.execute("CREATE INDEX sdbid_gj ON gj (sdbid);")

    cursor.execute("DROP TABLE IF EXISTS phot_tmp;")
    cursor.execute("CREATE TEMPORARY TABLE phot_tmp AS SELECT"
                   " id as sdbid, ROUND(-2.5*log10(model_jy/3882.37), 1) as Vmag"
                   " FROM "+cfg.db['db_results']+".phot WHERE filter='VJ' AND comp_no=-1;")
    cursor.execute("CREATE INDEX sdbid_phot ON phot_tmp (sdbid);")

    cursor.execute("DROP TABLE IF EXISTS star_tmp;")
    cursor.execute("CREATE TEMPORARY TABLE star_tmp AS SELECT"
                   " id, ROUND("+cfg.db['db_results']+".star.teff, 0) as teff, "
                   " plx_arcsec as plx_arcsec, lstar"
                   " from "+cfg.db['db_results']+".star"
                   " WHERE star_comp_no=0;")
    cursor.execute("CREATE INDEX id_star ON star_tmp (id);")

    cursor.execute("DROP TABLE IF EXISTS disk_r_tmp;")
    cursor.execute("CREATE TEMPORARY TABLE disk_r_tmp AS SELECT"
                   " id, GROUP_CONCAT(ROUND(temp, 1)) as temp, SUM(ldisk_lstar) as ldisk_lstar"
                   " from "+cfg.db['db_results']+".disk_r"
                   " GROUP BY id;")
    cursor.execute("CREATE INDEX id_dr ON disk_r_tmp (id);")
