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
        if r.pickle_time == db_mtime:
            return False

    return True


def write_phot(cursor,r):
    """Write photometry from a fit to db.
        
    Photometry is split by model component number, with -1 being total
    fluxes, and observed, "star", and "disk" photometry appearing only
    in these rows.

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
                    "(id,comp_no,filter,obs_jy,e_obs_jy,obs_upperlim,bibcode,"
                    "model_jy,e_model_jy_lo,e_model_jy_hi,residual,R) "
                    "VALUES ('{}',{:d},'{}',{:e}, "
                    "{:e},{},'{}',{:e},{:e},{:e}, "
                    "{:e},{:e})".format(cfg.mysql['phot_table'],
                                        r.id,-1,r.filters[j],r.obs_fnujy[j],
                                        r.obs_e_fnujy[j],r.obs_upperlim[j].astype(int),
                                        r.obs_bibcode[j],r.model_fnujy[j],
                                        r.model_fnujy_1sig_lo[j],r.model_fnujy_1sig_hi[j],
                                        r.residuals[j],ratio))

            if r.all_star_phot is not None:
                ratio = r.obs_fnujy[j]/r.all_star_phot[i]
                chi = (r.obs_fnujy[j] - r.all_star_phot[i]) / \
                            np.sqrt(r.obs_e_fnujy[j]**2 + \
                np.mean([r.all_star_phot_1sig_lo[i],r.all_star_phot_1sig_hi[i]])**2)

                if filter.iscolour(r.filters[j]):
                    ratio = r.obs_fnujy[j] - r.all_star_phot[i]

                cursor.execute("UPDATE {} "
                               "SET R_star = {:e}, chi_star = {:e}"
                               "WHERE id = '{}' AND comp_no = {:d} "
                               "AND filter = '{}'".format(cfg.mysql['phot_table'],
                                                          ratio,chi,r.id,-1,filt) )

        # the other filters
        else:
            cursor.execute("INSERT INTO {} "
                    "(id,comp_no,filter,model_jy,e_model_jy_lo,e_model_jy_hi) "
                    "VALUES ('{}',{:d},'{}',{:e},{:e}, "
                    "{:e})".format(cfg.mysql['phot_table'],
                                   r.id,-1,filt,r.all_phot[i],
                                   r.all_phot_1sig_lo[i],
                                   r.all_phot_1sig_hi[i]))

        # add these on to both
        if r.all_star_phot is not None:
            cursor.execute("UPDATE "+cfg.mysql['phot_table']+" "
                           "SET star_jy = {:e}, "
                           "e_star_jy_lo = {:e}, e_star_jy_hi = {:e} "
                           "WHERE id = '{}' AND comp_no = {:d} "
                           "AND filter = '{}'".format(r.all_star_phot[i],
                                                      r.all_star_phot_1sig_lo[i],
                                                      r.all_star_phot_1sig_hi[i],
                                                      r.id,-1,filt) )

        if r.all_disk_phot is not None:
            cursor.execute("UPDATE "+cfg.mysql['phot_table']+" "
                           "SET disk_jy = {:e}, "
                           "e_disk_jy_lo = {:e}, e_disk_jy_hi = {:e} "
                           "WHERE id = '{}' AND comp_no = {:d} "
                           "AND filter = '{}'".format(r.all_disk_phot[i],
                                                      r.all_disk_phot_1sig_lo[i],
                                                      r.all_disk_phot_1sig_hi[i],
                                                      r.id,-1,filt) )

        # fluxes for each model component
        for j in range(r.n_comps):
            cursor.execute("INSERT INTO {} "
                    "(id,comp_no,filter,model_jy,e_model_jy_lo,e_model_jy_hi) "
                    "VALUES ('{}',{:d},'{}',{:e},{:e}, "
                    "{:e})".format(cfg.mysql['phot_table'],r.id,j,filt,
                                   r.all_comp_phot[j,i],
                                   r.all_comp_phot_1sig_lo[j,i],
                                   r.all_comp_phot_1sig_hi[j,i]))



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
              str(r.dof),str(r.pickle_time))

    cursor.execute(stmt,values)


def write_star(cursor,r):
    """Write stellar properties to db.

    Assumes rows have been removed already (if necessary).
    """

    # loop over tuple of dicts of star results
    for i,star in enumerate(r.star):

        # start a row
        cursor.execute("INSERT INTO {} (id,star_comp_no) VALUES "
                       "('{}',{:d})".format(cfg.mysql['star_table'],
                                            r.id,star['comp_no']))

        # and add the rest of the results, only keys starting 'e_' and
        # without _lo or _hi on the end
        for key in star.keys():
            if key[:2] != 'e_':
                continue
            if key[-3:] == '_hi' or key[-3:] == '_lo':
                continue
            cursor.execute(
               "UPDATE {} SET {} = {:e}, {} = {:e} WHERE id = '{}' "
               "AND star_comp_no = {}".format(cfg.mysql['star_table'],
                                              key[2:],star[key[2:]],
                                              key,star[key],
                                              r.id,star['comp_no'])
                           )


def write_disk_r(cursor,r):
    """Write narrow disk properties to db.

    Assumes rows have been removed already (if necessary).
    """

    # loop over tuple of dicts of disk_r results
    for i,disk_r in enumerate(r.disk_r):
        
        cursor.execute("INSERT INTO {} (id,disk_r_comp_no) VALUES "
                       "('{}',{:d})".format(cfg.mysql['disk_r_table'],
                                            r.id,disk_r['comp_no']))

        # and add the rest of the results, only keys starting 'e_' and
        # without _lo or _hi on the end
        for key in disk_r.keys():
            if key[:2] != 'e_':
                continue
            if key[-3:] == '_hi' or key[-3:] == '_lo':
                continue
            cursor.execute(
               "UPDATE {} SET {} = {:e}, {} = {:e} WHERE id = '{}' "
               "AND disk_r_comp_no = {}".format(cfg.mysql['disk_r_table'],
                                         key[2:],disk_r[key[2:]],
                                         key,disk_r[key],
                                         r.id,disk_r['comp_no'])
                           )


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


def get_sdbids(ids):
    """Return a list of sdbids given a list of ids.

    Parameters
    ----------
    ids : list
        List of target ids.
    """

    if not isinstance(ids,list):
        raise utils.SdfError('please pass a list, not {}'.format(type(ids)))

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

    sdbids = []
    for id in ids:

        cursor.execute("SELECT sdbid FROM xids WHERE xid = '{}';".format(id))
        if cursor.rowcount == 0:
            sdbids.append( None )
        else:
            for (sdbid,) in cursor:
                sdbids.append( sdbid )

    return sdbids
