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
        
        cursor.execute("DELETE FROM {} "
                       "WHERE id = '{}'".format(cfg.mysql['phot_table'],r.id))
        write_phot(cursor,r)

        cursor.execute("DELETE FROM {} "
                       "WHERE id = '{}'".format(cfg.mysql['model_table'],r.id))
        write_model(cursor,r)

        cursor.execute("DELETE FROM {} "
                       "WHERE id = '{}'".format(cfg.mysql['star_table'],r.id))
        write_star(cursor,r)

        cursor.execute("DELETE FROM {} "
                       "WHERE id = '{}'".format(cfg.mysql['disk_r_table'],r.id))
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
    cursor.execute("SELECT MIN(model_mtime) "
                   "FROM {} WHERE id = '{}';".format(table,r.id))
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

            # compute R properly, / for flux and - for colurs
            ratio = r.obs_fnujy[j] / r.model_fnujy[j]
            if filter.iscolour(r.filters[j]):
                ratio = r.obs_fnujy[j] - r.model_fnujy[j]
            
            cursor.execute("INSERT INTO {} "
                    "(id,filter,obs_jy,e_obs_jy,obs_upperlim,bibcode,"
                    "model_jy,e_model_jy_lo,e_model_jy_hi,chi,R) "
                    "VALUES ('{}','{}',{:e}, "
                    "{:e},{},'{}',{:e},{:e},{:e}, "
                    "{:e},{:e})".format(cfg.mysql['phot_table'],
                                        r.id,r.filters[j],r.obs_fnujy[j],
                                        r.obs_e_fnujy[j],r.obs_upperlim[j].astype(int),
                                        r.obs_bibcode[j],r.model_fnujy[j],
                                        r.model_fnujy_1sig_lo[j],r.model_fnujy_1sig_hi[j],
                                        r.residuals[j],ratio))
                    
        # the other filters
        else:

            cursor.execute("INSERT INTO {} "
                    "(id,filter,model_jy,e_model_jy_lo,e_model_jy_hi) "
                    "VALUES ('{}','{}',{:e},{:e}, "
                    "{:e})".format(cfg.mysql['phot_table'],
                                   r.id,filt,r.all_phot[i],
                                   r.all_phot_1sig_lo[i],
                                   r.all_phot_1sig_hi[i]))

        # add these on to both
        if r.star_phot is not None:
            cursor.execute("UPDATE "+cfg.mysql['phot_table']+" "
                           "SET star_jy = {:e}, "
                           "e_star_jy_lo = {:e}, e_star_jy_hi = {:e} "
                           "WHERE id = '{}' "
                           "AND filter = '{}'".format(r.star_phot[i],
                                                      r.star_phot_1sig_lo[i],
                                                      r.star_phot_1sig_hi[i],
                                                      r.id,filt) )

        if r.disk_phot is not None:
            cursor.execute("UPDATE "+cfg.mysql['phot_table']+" "
                           "SET disk_jy = {:e}, "
                           "e_disk_jy_lo = {:e}, e_disk_jy_hi = {:e} "
                           "WHERE id = '{}' "
                           "AND filter = '{}'".format(r.disk_phot[i],
                                                      r.disk_phot_1sig_lo[i],
                                                      r.disk_phot_1sig_hi[i],
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
        cursor.execute("INSERT INTO {} (id,teff,e_teff) VALUES "
                       "('{}',{:e},{:e})".format(cfg.mysql['star_table'],
                                                 r.id,star['Teff'],
                                                 star['e_Teff']))

        # and add the rest of the results
        for key in star.keys():
            if key.find('e_') == 0:
                continue
            cursor.execute("UPDATE {} SET {} = {:e}, {} = {:e} WHERE "
                           "id = '{}'".format(cfg.mysql['star_table'],
                                              key,star[key],
                                              'e_'+key,star['e_'+key],
                                              r.id) )


def write_disk_r(cursor,r):
    """Write narrow disk properties to db.

    Assumes rows have been removed already (if necessary).
    """

    # loop over tuple of dicts of disk_r results
    for i,disk_r in enumerate(r.disk_r):
        
        cursor.execute("INSERT INTO {} (id,comp_no) VALUES "
                       "('{}',{})".format(cfg.mysql['disk_r_table'],r.id,i) )

        # and add the rest of the results
        for key in disk_r.keys():
            if key.find('e_') == 0:
                continue
            cursor.execute("UPDATE {} SET {} = {:e}, {} = {:e} WHERE id = '{}' "
                           "AND comp_no = {}".format(cfg.mysql['disk_r_table'],
                                                     key,disk_r[key],
                                                     'e_'+key,disk_r['e_'+key],
                                                     r.id,i) )


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


def sdb_xids(id):
    """Get selected xids for a given sdb id."""

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

    cursor.execute("SELECT xid FROM xids WHERE sdbid = '{}' AND xid "
                   "REGEXP('^HD|^HR|^HIP|^GJ|^TYC|^NAME|^\\\\* ')".format(id))
    xids = []
    for (id,) in cursor:
        xids.append(id)

    if len(xids) == 0:
        return
    else:
        return xids
