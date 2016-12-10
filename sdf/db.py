import numpy as np

import mysql.connector
import astropy.units as u

from . import filter
from . import utils
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

        stmt = ("DELETE FROM "+cfg.mysql['star_table']+" WHERE id = %(a)s;")
        cursor.execute(stmt,{'a':str(r.id)})
        write_star(cursor,r)

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
    
    TODO: allow for multiple stellar components.
    """

    # loop over plotting models, only one SpecModel per component
    for i,comp in enumerate(r.pl_models):
        if 'kurucz' in r.model_comps[i] or 'phoenix' in r.model_comps[i]:

            cursor.execute("INSERT INTO "+cfg.mysql['star_table']+" "
                           "(id) VALUES (%s)",(str(r.id),))

            # parameters directly from model
            for j,par in enumerate(r.parameters):
                if par == 'norm' or par == 'spec_norm':
                    continue
                cursor.execute("UPDATE "+cfg.mysql['star_table']+" "
                               "SET {} = {:e} WHERE id = '{}'".format(
                               par,r.best_params[j],r.id) )

            # derived parameters
            lstar_1pc = r.comp_spectra[i].irradiance \
                        * 4 * np.pi * (u.pc.to(u.m))**2 / u.L_sun.to(u.W)

            cursor.execute("UPDATE "+cfg.mysql['star_table']+" "
                           "SET lstar_1pc = {:e} WHERE id = '{}'".format(
                           lstar_1pc,r.id) )

            # distance-dependent params
            if r.obs_keywords['plx_value'] is not None:
                plx_arcsec = r.obs_keywords['plx_value'] / 1e3
                cursor.execute("UPDATE "+cfg.mysql['star_table']+" "
                               "SET lstar = {:e} WHERE id = '{}'".format(
                               lstar_1pc / plx_arcsec**2,r.id) )
                               
                cursor.execute("UPDATE "+cfg.mysql['star_table']+" "
                               "SET plx_arcsec = {} WHERE id = '{}'".format(
                               plx_arcsec,r.id) )

                rstar = np.sqrt(cfg.ssr * 10**r.comp_best_params[i][-1]/np.pi) \
                        * u.pc.to(u.m) / plx_arcsec / u.R_sun.to(u.m)

                cursor.execute("UPDATE "+cfg.mysql['star_table']+" "
                               "SET rstar = {:e} WHERE id = '{}'".format(
                               rstar,r.id) )


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
