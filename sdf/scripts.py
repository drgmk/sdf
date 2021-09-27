import os
import glob
import argparse

#from pympler import summary, muppy, tracker

import filelock

from sdf import fitting
from sdf import plotting
from sdf import tables
from sdf import db
from sdf import www
from sdf import config as cfg


def sdf_fit():
    """Run fitting with sdf."""

#    tr = tracker.SummaryTracker()

    # inputs
    parser1 = argparse.ArgumentParser(description='Run sdf')
    parser = parser1.add_mutually_exclusive_group(required=True)
    parser.add_argument('--file','-f',nargs='+',action='append',
                        help='Fit SED to file or files')
    parser.add_argument('--dir','-d',nargs='+',action='append',
                        help='Fit SED to *-rawphot.txt files in path(s)')
    parser.add_argument('--samples','-s',nargs='+',
                        help='Restrict to target samples')

    parser1.add_argument('--subset',nargs='+',default='*',
                         help='Restrict to subset of targets (e.g. public)')

    parser1.add_argument('--www','-w',action='store_true',
                         help='www material',default=False)
    parser1.add_argument('--update-www',action='store_true',
                         help='Force update of www')
                         
    parser1.add_argument('--dbwrite','-b',action='store_true',
                         help='Write results to db',default=False)
    parser1.add_argument('--update-db',action='store_true',
                         help='Force udpate of db')

    parser1.add_argument('--no-spectra',action='store_true',
                         help='Exclude spectra from fitting')

    parser1.add_argument('--update-all','-u',action='store_true',
                         help='Force udpate of everything')

    parser1.add_argument('--update-analysis','-a',action='store_true',
                         help='Force udpate of post-multinest analysis')

    parser1.add_argument('--update-json','-j',action='store_true',
                         help='Force udpate of json file')

    parser1.add_argument('--update-thumb','-t',action='store_true',
                         help='Force udpate of SED thumbnail image')

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

    elif args.samples is not None:
        ids = []
        for s in args.samples:
            ids += db.sample_targets(s)
        files = []
        for id in ids:
            files += glob.glob( cfg.file['sdb_root']+'masters/'\
                               +id+'/**/*-rawphot.txt' )

    # one at a time, locking before we start
    for f in files:
        
        lock = filelock.FileLock(os.path.dirname(os.path.abspath(f))
                                 +'/.sdf_lock-'
                                 +os.path.basename(f))
        try:
            with lock.acquire(timeout = 0):
                
                print(f)

                # evidence-sorted list of results
                results = fitting.fit_results(os.path.abspath(f),
                                              update_mn=args.update_all,
                                              update_an=args.update_analysis,
                                              update_json=args.update_json,
                                              update_thumb=args.update_thumb,
                                              nospec=args.no_spectra)
                if results is None:
                    continue

                if args.www or args.update_www:
                    www.www_all(results,update=args.update_www)

                # write best model to db
                if args.dbwrite or args.update_db:
                    db.write_all(results[0],update=args.update_db)

        except filelock.Timeout:
            pass

#        tr.print_diff()


def sdf_sample():
    """Generate HTML sample pages to browse database."""
    
    # inputs
    parser = argparse.ArgumentParser(description='Update sample www pages')
    parser.add_argument('--samples','-s',nargs='+',
                        help='Restrict to target samples')
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
        tables.sample_tables(args.samples)

    if args.plots:
        print("Updating sample plots")
        plotting.flux_size_plots(args.samples)
        plotting.sample_plots(args.samples)

    if args.calibration:
        print("Updating calibration plots")
        for sample in cfg.www['cal_samples']:
            plotting.calibration(sample=sample)

    if args.cleanup:
        print("Cleaning up")
        www.cleanup_sample_dirs()
        www.cleanup_calibration_dirs()


def sdf_cleanup():
    """Cleanup."""

    parser = argparse.ArgumentParser(description='Remove multinest files')
    parser.add_argument('--sample','-s',nargs='+',
                        help='Sample(s) to clean up')
    args = parser.parse_args()

    if args.sample is not None:
        for s in args.sample:

            ids = db.sample_targets(s)
            fs = []
            for id in ids:
                fs += glob.glob( cfg.file['sdb_root']+'masters/'\
                                +id+'/*/*mnest/*' )

            removed = [os.unlink(f) for f in fs]
