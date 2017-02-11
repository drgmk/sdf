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
    
    For uncertainty propagation see examples at:
        https://en.wikipedia.org/wiki/Propagation_of_uncertainty
    """
    # TODO: put all the derived stuff in a function elsewhere

    # loop over plotting models, only one SpecModel per component
    for i,comp in enumerate(r.model_comps):
        if comp in cfg.models['star']:

            # find the temperature for this component
            for j,par in enumerate(r.parameters):
                if par == 'Teff':
                    teff = r.best_params[j]
                    e_teff = r.best_params_1sig[j]
        
            cursor.execute("INSERT INTO "+cfg.mysql['star_table']+" "
                           "(id,teff,e_teff) VALUES "
                           "('{}',{:e},{:e})".format(r.id,teff,e_teff))

            # parameters directly from model
            for j,par in enumerate(r.comp_parameters[i]):
                if par == 'Teff' or par == 'norm' or par == 'spec_norm':
                    continue
                cursor.execute("UPDATE "+cfg.mysql['star_table']+" "
                               "SET {} = {:e}, e_{} = {:e} WHERE "
                               "id = '{}'".format(par,r.comp_best_params[i][j],
                                                  par,r.comp_best_params_1sig[i][j],
                                                  r.id) )

            # stellar luminosity at 1pc, uncertainty is normalisation
            frac_norm = np.log(10) * r.comp_best_params_1sig[i][-1]
            lstar_1pc = r.comp_spectra[i].irradiance \
                        * 4 * np.pi * (u.pc.to(u.m))**2 / u.L_sun.to(u.W)
            e_lstar_1pc = lstar_1pc * frac_norm

            cursor.execute("UPDATE "+cfg.mysql['star_table']+" "
                           "SET lstar_1pc = {:e},e_lstar_1pc = {:e} "
                           "WHERE id = '{}'".format(lstar_1pc,
                                                    e_lstar_1pc,r.id) )

            # distance-dependent params
            if r.obs_keywords['plx_value'] is not None:
                if r.obs_keywords['plx_value'] > 0:
                    plx_arcsec = r.obs_keywords['plx_value'] / 1e3
                    
                    if r.obs_keywords['plx_err'] is not None:
                        e_plx_arcsec = r.obs_keywords['plx_err'] / 1e3
                    else:
                        e_plx_arcsec = plx_arcsec / 3.
                    
                    lstar = lstar_1pc / plx_arcsec**2
                    e_lstar = lstar * np.sqrt( frac_norm**2
                                              + (2*e_plx_arcsec/plx_arcsec)**2 )
                                      
                    cursor.execute("UPDATE "+cfg.mysql['star_table']+" "
                                   "SET lstar = {:e},e_lstar = {:e} WHERE "
                                   "id = '{}'".format(lstar,e_lstar,r.id) )
                                   
                    cursor.execute("UPDATE "+cfg.mysql['star_table']+" "
                                   "SET plx_arcsec = {:e},e_plx_arcsec = {:e} WHERE "
                                   "id = '{}'".format(plx_arcsec,e_plx_arcsec,r.id) )

                    rstar = np.sqrt(cfg.ssr * 10**r.comp_best_params[i][-1]/np.pi) \
                            * u.pc.to(u.m) / plx_arcsec / u.R_sun.to(u.m)
                    e_rstar = rstar * ( np.sqrt( frac_norm**2
                                                + (2*e_plx_arcsec/plx_arcsec)**2) )

                    cursor.execute("UPDATE "+cfg.mysql['star_table']+" "
                                   "SET rstar = {:e},e_rstar = {:e} WHERE "
                                   "id = '{}'".format(rstar,e_rstar,r.id) )


def write_disk_r(cursor,r):
    """Write narrow disk properties to db.

    Assumes rows have been removed already (if necessary).

    For uncertainty propagation see examples at:
        https://en.wikipedia.org/wiki/Propagation_of_uncertainty
    """
    # TODO: put all the derived stuff in a function elsewhere

    # loop over plotting models, only one SpecModel per component
    comp_no = 0
    for i,comp in enumerate(r.model_comps):
        if comp in cfg.models['disk_r']:
            
            cursor.execute("INSERT INTO "+cfg.mysql['disk_r_table']+" "
                           "(id,comp_no) VALUES "
                           "('{}',{})".format(r.id,comp_no) )

            # parameters directly from model, keep disk temp and e_temp
            for j,par in enumerate(r.comp_parameters[i]):
                if par == 'norm' or par == 'spec_norm':
                    continue
                elif 'log_' in par:
                    par_in = par.replace('log_','')
                    val = 10**r.comp_best_params[i][j]
                    e_val = 10**r.comp_best_params_1sig[i][j]
                    if par == 'log_Temp':
                        tdisk = val
                        e_tdisk = e_val
                else:
                    par_in = par
                    val = r.comp_best_params[i][j]
                    e_val = r.comp_best_params_1sig[i][j]

                cursor.execute("UPDATE "+cfg.mysql['disk_r_table']+" "
                               "SET {} = {:e}, e_{} = {:e} WHERE "
                               "id = '{}' AND comp_no = {}".\
                               format(par_in,val,par_in,e_val,r.id,comp_no) )

            # disk and fractional luminosity
            frac_norm = np.log(10) * r.comp_best_params_1sig[i][-1]
            ldisk_1pc = r.comp_spectra[i].irradiance \
                        * 4 * np.pi * (u.pc.to(u.m))**2 / u.L_sun.to(u.W)
            e_ldisk_1pc = ldisk_1pc * frac_norm

            cursor.execute("SELECT SUM(lstar_1pc) FROM "+cfg.mysql['star_table']+" "
                           "WHERE id = '{}'".format(r.id) )
            lstar_1pc = cursor.fetchall()[0][0]
            
            cursor.execute("SELECT SUM(e_lstar_1pc) FROM "
                           +cfg.mysql['star_table']+" "
                           "WHERE id = '{}'".format(r.id) )
            e_lstar_1pc = cursor.fetchall()[0][0]
            frac_lstar_1pc = e_lstar_1pc / lstar_1pc

            ldls = ldisk_1pc/lstar_1pc
            e_ldls = ldls * np.sqrt(frac_norm**2 + frac_lstar_1pc**2)

            cursor.execute("UPDATE "+cfg.mysql['disk_r_table']+" "
                           "SET ldisk_1pc = {:e},e_ldisk_1pc = {:e}, "
                           "ldisk_lstar = {:e},e_ldisk_lstar = {:e} WHERE "
                           "id = '{}' AND comp_no = {}".\
                           format(ldisk_1pc,e_ldisk_1pc,ldls,e_ldls,r.id,comp_no) )

            # distance-dependent params
            if r.obs_keywords['plx_value'] is not None:
                plx_arcsec = r.obs_keywords['plx_value'] / 1e3
                rdisk = (lstar_1pc/plx_arcsec**2)**0.5 * (278.3/tdisk)**2
                e_rdisk = rdisk * np.sqrt( (0.5*frac_lstar_1pc)**2
                                          + (2*(e_tdisk/tdisk))**2 )
                cursor.execute("UPDATE "+cfg.mysql['disk_r_table']+" "
                               "SET rdisk = {:e},e_rdisk = {:e} WHERE "
                               "id = '{}' AND comp_no = {}".\
                               format(rdisk,e_rdisk,r.id,comp_no) )

            comp_no += 1


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
