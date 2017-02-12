import numpy as np

import mysql.connector
import astropy.units as u

from . import filter
from . import utils
from . import config as cfg


def write_all(r,update=False):
    """Write sdf results to db."""

    print(" Database")

    # set up connection
    try:
        cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                      password=cfg.mysql['passwd'],
                                      host=cfg.mysql['host'],
                                      database=cfg.mysql['db_results'])
        cursor = cnx.cursor(buffered=True)

    except mysql.connector.InterfaceError:
        print("   Can't connect to {} at {}".format(cfg.mysql['db_results'],
                                                    cfg.mysql['host']) )
        return

    # write to the tables
    if update or update_needed(cursor,cfg.mysql['model_table'],r):
        
        print("   writing")
        
        stmt = ("DELETE FROM "+cfg.mysql['phot_table']+" WHERE id = %(a)s;")
        cursor.execute(stmt,{'a':str(r.id)})
        write_phot(cursor,r)

        stmt = ("DELETE FROM "+cfg.mysql['model_table']+" WHERE id = %(a)s;")
        cursor.execute(stmt,{'a':str(r.id)})
        write_model(cursor,r)

        cursor.execute("DELETE FROM "+cfg.mysql['star_table']+" "
                       "WHERE id = '{}'".format(r.id))
        write_star(cursor,r)

        cursor.execute("DELETE FROM "+cfg.mysql['disk_r_table']+" "
                       "WHERE id = '{}'".format(r.id))
        write_disk_r(cursor,r)
    
    else:
        print("   no update needed")

    # commit and close
    cnx.commit()
    cursor.close()
    cnx.close()


def update_needed(cursor,table,r):
    """Check mtime against result, deleting old entries."""

    # if entries in table, get file modification time
    stmt = ("SELECT MIN(model_mtime) FROM "
            +table+" WHERE id = %(a)s;")
    cursor.execute(stmt,{'a':str(r.id)})
    db_mtime = cursor.fetchall()[0][0]

    # False if sql up to date
    if db_mtime is not None:
        if r.mtime <= db_mtime:
            return False

    return True


def write_phot(cursor,r):
    """Write photometry from a fit to db.
        
    Assumes rows have been removed already (if necessary).
    """

    # write to table
    for i,filt in enumerate(r.all_filters):

        # filters with observations
        if filt in r.filters:
            
            j = np.where(filt == r.filters)[0][0]

            stmt = ("INSERT INTO "+cfg.mysql['phot_table']+" "
                    "(id,filter,obs_jy,e_obs_jy,obs_upperlim,bibcode,"
                    "model_jy,chi,R) "
                    "VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)")
                    
            # compute R properly, / for flux and - for colurs
            ratio = r.obs_fnujy[j] / r.model_fnujy[j]
            if filter.iscolour(r.filters[j]):
                ratio = r.obs_fnujy[j] - r.model_fnujy[j]
            
            values = (str(r.id),str(r.filters[j]),str(r.obs_fnujy[j]),
                      str(r.obs_e_fnujy[j]),str(r.obs_upperlim[j].astype(int)),
                      str(r.obs_bibcode[j]),str(r.model_fnujy[j]),
                      str(r.residuals[j]),str(ratio))
            
            cursor.execute(stmt,values)

        # the other filters
        else:

            stmt = ("INSERT INTO "+cfg.mysql['phot_table']+" "
                    "(id,filter,model_jy) "
                    "VALUES (%s,%s,%s)")
            values = (str(r.id),str(filt),
                      str(r.all_phot[i]) )
            cursor.execute(stmt,values)

        # add these on to both
        if r.star_phot is not None:
            cursor.execute("UPDATE "+cfg.mysql['phot_table']+" "
                           "SET star_jy = {:e} WHERE id = '{}' "
                           "AND filter = '{}'".format(r.star_phot[i],
                                                    r.id,filt) )

        if r.disk_phot is not None:
            cursor.execute("UPDATE "+cfg.mysql['phot_table']+" "
                           "SET disk_jy = {:e} WHERE id = '{}' "
                           "AND filter = '{}'".format(r.disk_phot[i],
                                                    r.id,filt) )


def write_model(cursor,r):
    """Write model results to db.
    
    Assumes rows have been removed already (if necessary).
    """

    stmt = ("INSERT INTO "+cfg.mysql['model_table']+" "
            "(id,model_comps,parameters,evidence,chisq,dof,model_mtime) "
            "VALUES (%s,%s,%s,%s,%s,%s,%s)")

    values = (str(r.id),str(r.model_comps),
              '['+','.join("{:e}".format(p) for p in r.best_params)+']',
              str(r.evidence),str(r.chisq),
              str(r.dof),str(r.mtime))

    cursor.execute(stmt,values)


def write_star(cursor,r):
    """Write stellar properties to db.

    Assumes rows have been removed already (if necessary).
    """

    # loop over tuple of dicts of star results
    for star in r.star:

        # start a row
        cursor.execute("INSERT INTO "+cfg.mysql['star_table']+" "
                       "(id,teff,e_teff) VALUES "
                       "('{}',{:e},{:e})".format(r.id,star['Teff'],
                                                 star['e_Teff']))

        # and add the rest of the results
        for key in star.keys():
            if key.find('e_') == 0:
                continue
            cursor.execute("UPDATE "+cfg.mysql['star_table']+" "
                           "SET {} = {:e}, {} = {:e} WHERE "
                           "id = '{}'".format(key,star[key],
                                              'e_'+key,star['e_'+key],
                                              r.id) )


def write_disk_r(cursor,r):
    """Write narrow disk properties to db.

    Assumes rows have been removed already (if necessary).
    """

    # loop over tuple of dicts of disk_r results
    for i,disk_r in enumerate(r.disk_r):
        
        cursor.execute("INSERT INTO "+cfg.mysql['disk_r_table']+" "
                       "(id,comp_no) VALUES "
                       "('{}',{})".format(r.id,comp_no) )

        # and add the rest of the results
        for key in disk_r.keys():
            if key.find('e_') == 0:
                continue
            cursor.execute("UPDATE "+cfg.mysql['disk_r_table']+" "
                           "SET {} = {:e}, {} = {:e} WHERE "
                           "id = '{}'".format(key,disk_r[key],
                                              'e_'+key,disk_r['e_'+key],
                                              r.id) )


def sample_targets(sample,db='sdb_samples'):
    """Return list of sdbids of targets in some sample."""

    # set up connection
    try:
        cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                      password=cfg.mysql['passwd'],
                                      host=cfg.mysql['host'],
                                      database=db)
        cursor = cnx.cursor(buffered=True)

    except mysql.connector.InterfaceError:
        print("Can't connect to {} at {}".format(cfg.mysql[db],
                                                 cfg.mysql['host']) )
        return

    cursor.execute("SELECT sdbid FROM "+sample+" WHERE sdbid IS NOT NULL")
    ids = []
    for (id,) in cursor:
        ids.append(id)

    return ids


def sdb_info(id):
    """Get info for a given sdb id."""

    # set up connection
    try:
        cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                      password=cfg.mysql['passwd'],
                                      host=cfg.mysql['host'],
                                      database=cfg.mysql['db_sdb'])
        cursor = cnx.cursor(buffered=True)

    except mysql.connector.InterfaceError:
        print("Can't connect to {} at {}".format(cfg.mysql['db_sdb'],
                                                 cfg.mysql['host']) )
        return

    cursor.execute("SELECT sdbid,COALESCE(main_id,sdbid), "
                   "raj2000,dej2000, "
                   "COALESCE(gaia.pmra,sdb_pm.pmra), "
                   "COALESCE(gaia.pmde,sdb_pm.pmde), "
                   "1e3/COALESCE(gaia.plx,simbad.plx_value)"
                   "FROM sdb_pm LEFT JOIN simbad USING (sdbid) "
                   "LEFT JOIN gaia USING (sdbid) "
                   "WHERE sdbid = '{}'".format(id))
    out = cursor.fetchall()
    if len(out) > 0:
        sdbid,main_id,ra,dec,pmra,pmde,dist = out[0]
    else:
        return

    cursor.execute("SELECT xid FROM xids WHERE sdbid = '{}'".format(id))
    xids = []
    for (id,) in cursor:
        xids.append(id)

    return sdbid,main_id,xids,ra,dec,dist
