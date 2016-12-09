import numpy as np
import mysql.connector

from . import filter
from . import config as cfg


def write_all(r,update=False):
    """Write sdf results to db."""

    # set up connection
    try:
        cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                      password=cfg.mysql['passwd'],
                                      host=cfg.mysql['host'],
                                      database=cfg.mysql['db_results'])
        cursor = cnx.cursor(buffered=True)
        print(" Writing to db")

    except mysql.connector.InterfaceError:
        print("Can't connect to {} at {}".format(cfg.mysql['db_results'],
                                                 cfg.mysql['host']) )
        return

    # write to the tables
    if update or update_needed(cursor,cfg.mysql['model_table'],r):
        
        stmt = ("DELETE FROM "+cfg.mysql['phot_table']+" WHERE id = %(a)s;")
        cursor.execute(stmt,{'a':str(r.id)})
        write_phot(cursor,r)

        stmt = ("DELETE FROM "+cfg.mysql['model_table']+" WHERE id = %(a)s;")
        cursor.execute(stmt,{'a':str(r.id)})
        write_model(cursor,r)

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
        
    Assumes rows have been removed already by update_needed (if
    necessary).
    """

    # write to table
    for i in range(len(r.filters)):

        # skip spectra, for which filters[i] is None
        if r.filters[i] is None:
            continue
        
        stmt = ("INSERT INTO "+cfg.mysql['phot_table']+" "
                "(id,filter,obs_jy,e_obs_jy,obs_upperlim,bibcode,"
                "model_jy,comps_jy,chi,R) "
                "VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)")
                
        # compute R properly, / for flux and - for colurs
        ratio = r.obs_fnujy[i] / r.model_fnujy[i]
        if filter.iscolour(r.filters[i]):
            ratio = r.obs_fnujy[i] - r.model_fnujy[i]
        
        comp_str = '['+','.join("{:e}".format(p) \
                                for p in r.model_comp_fnujy[:,i])+']'
        
        values = (str(r.id),str(r.filters[i]),str(r.obs_fnujy[i]),
                  str(r.obs_e_fnujy[i]),str(r.obs_upperlim[i].astype(int)),
                  str(r.obs_bibcode[i]),str(r.model_fnujy[i]),
                  comp_str,str(r.residuals[i]),str(ratio))
        
        cursor.execute(stmt,values)


def write_model(cursor,r):
    """Write model results to db.
    
    Assumes rows have been removed already by update_needed (if
    necessary).
    """

    stmt = ("INSERT INTO "+cfg.mysql['model_table']+" "
            "(id,model_comps,parameters,evidence,chisq,dof,model_mtime) "
            "VALUES (%s,%s,%s,%s,%s,%s,%s)")

    values = (str(r.id),str(r.model_comps),
              '['+','.join("{:e}".format(p) for p in r.best_params)+']',
              str(r.evidence),str(r.chisq),
              str(r.dof),str(r.mtime))

    cursor.execute(stmt,values)


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

    cursor.execute("SELECT sdbid FROM "+sample)
    ids = []
    for (id,) in cursor:
        ids.append(id)

    return ids
