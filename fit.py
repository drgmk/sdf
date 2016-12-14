import os
import glob
import pickle
import argparse
from multiprocessing import Pool

import numpy as np
import corner
import pymultinest as pmn
import filelock

from sdf import result
from sdf import plotting
from sdf import db
from sdf import config as cfg


def fit_results(file,update=False,sort=True,nospec=False):
    """Return a list of fitting results.
        
    If necessary, the fitting will be done, otherwise the results are
    just loaded and returned.
    
    By default sort the list by evidence.
    """

    print(" Fitting")

    # load models, fit is a list of model combos
    fit = cfg.fitting['models']

    results = []
    for f in fit:

        print("  ",f)
        
        res = result.Result(file,f,update=update,nospec=nospec)

        if hasattr(res,'obs'):
            results.append(res)
        else:
            print("  no photometry = no results")
            return None

        # if evidence didn't go up with second model component, stop
        if len(f) > 1:
            evs = [r.evidence for r in results]
            if np.max(evs) != evs[-1]:
                break

    if sort:
        results = [results[i] for i in result.sort_results(results)]

    return results


def plot_seds(results,update=False):
    """Create SED plots plot of fitting results.
        
    Collects all model fits and puts them in tabs in a plot.
    """
    
    print(" Plotting")

    # see whether the sed needs updating (unless update is enforced)
    if os.path.exists(results[0].sed_plot):
        pkltime = []
        for r in results:
            pkltime.append( os.path.getmtime(r.pickle) )

        if os.path.getmtime(r.sed_plot) > np.max(pkltime):
            if not update:
                print("  no update needed")
                return
        print("  updating sed")
    else:
        print("  generating sed")

    plotting.sed(results,file=results[0].sed_plot)


def everything(f,up_res,nospec,up_plot,up_db):
    """Wrapper to do everything for parallelism."""

    print(f)

    # evidence-sorted list of results
    results = fit_results( os.path.abspath(f),update=up_res,
                           nospec=nospec)
    if results is None:
        return

    if args.plot:
        plot_seds(results,update=up_plot)

    # write best model to db
    if args.dbwrite:
        db.write_all(results[0],update=up_db)


# command line
if __name__ == '__main__':
    """Run it, with options."""

    # inputs
    parser1 = argparse.ArgumentParser(description='Fit SED models to sdb    \
                                                   rawphot files')
    parser = parser1.add_mutually_exclusive_group(required=True)
    parser.add_argument('--file','-f',nargs='+',action='append',
                        help='Fit SED to file or files')
    parser.add_argument('--dir','-d',nargs='+',action='append',
                        help='Fit SED to *-rawphot.txt files in path(s)')
    parser.add_argument('--sample','-s',nargs='+',
                        help='Restrict to target sample')

    parser1.add_argument('--subset',nargs='+',default='*',
                         help='Restrict to subset of targets (e.g. public)')

    parser1.add_argument('--plot','-p',action='store_true',
                         help='Plot SEDs')
    parser1.add_argument('--update-plot',action='store_true',
                         help='Force update of SEDs')
                         
    parser1.add_argument('--dbwrite','-w',action='store_true',
                         help='Write results to db')
    parser1.add_argument('--update-db',action='store_true',
                         help='Force udpate of db')

    parser1.add_argument('--no-spectra',action='store_true',
                         help='Exclude spectra from fitting')

    parser1.add_argument('--update-all','-u',action='store_true',
                         help='Force udpate of everything')

    args = parser1.parse_args()
    
    # collect the files
    if args.file is not None:
        files = args.file[0]

    elif args.dir is not None:
        for d in args.dir[0]:
            d.rstrip('/')
            files = glob.glob(os.path.abspath(d)+'/**/'\
                              +args.subset[0]+'/*-rawphot.txt',
                              recursive=True)

    elif args.sample is not None:
        ids = []
        for s in args.sample:
            ids += db.sample_targets(s)
        files = []
        for id in ids:
            files += glob.glob( cfg.file['sdb_root']+'masters/'\
                               +id+'/**/*-rawphot.txt' )

    # one at a time, locking before we start
    for f in files:
        
        lock = filelock.FileLock(os.path.dirname(f)
                                 +'/.sdf_lock-'
                                 +os.path.basename(f))
        try:
            with lock.acquire(timeout = 0):
                
                print(f)

                # evidence-sorted list of results
                results = fit_results(os.path.abspath(f),
                                      update=args.update_all,
                                      nospec=args.no_spectra)
                if results is None:
                    continue

                if args.plot:
                    plot_seds(results,update=args.update_plot)

                # write best model to db
                if args.dbwrite:
                    db.write_all(results[0],update=args.update_db)

        except filelock.Timeout:
            pass

# parallel, this looks like a big fail due to something forking
#    args = [(f,args.update_all,args.no_spectra,
#             args.update_plot,args.update_db) for f in files]
#    print(args)
#    with Pool(2) as p:
#        p.starmap(everything,args)

