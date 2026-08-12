import os
import glob
import argparse

# from pympler import summary, muppy, tracker

import filelock

from sdf import fitting
from sdf import plotting
from sdf import tables
from sdf import db
from sdf import www
from sdf import config as cfg
from sdf.dataset import discover_fit_sources
from sdf.coordinator import fit_package


def sdf_fit():
    """Run fitting with sdf."""

#    tr = tracker.SummaryTracker()

    # inputs
    parser1 = argparse.ArgumentParser(description='Run sdf')
    parser = parser1.add_mutually_exclusive_group(required=True)
    parser.add_argument('--file', '-f', nargs='+', action='append',
                        help='Fit SED to file or files')
    parser.add_argument('--dir', '-d', nargs='+', action='append',
                        help='Fit SED to *-rawphot.txt files in path(s)')
    parser.add_argument('--samples', '-s', nargs='+',
                        help='Restrict to target samples')

    parser1.add_argument('--subset', nargs='+', default='*',
                         help='Restrict to subset of targets (e.g. public)')

    parser1.add_argument('--www', '-w', action='store_true',
                         help='www material', default=False)
    parser1.add_argument('--update-www', action='store_true',
                         help='Force update of www')
                         
    parser1.add_argument('--dbwrite', '-b', action='store_true',
                         help='Write results to db', default=False)
    parser1.add_argument('--update-db', action='store_true',
                         help='Force udpate of db')

    parser1.add_argument('--no-spectra', action='store_true',
                         help='Exclude spectra from fitting')

    parser1.add_argument('--update-all', '-u', action='store_true',
                         help='Force udpate of everything')

    parser1.add_argument('--update-analysis', '-a', action='store_true',
                         help='Force udpate of post-multinest analysis')

    parser1.add_argument('--update-json', '-j', action='store_true',
                         help='Force udpate of json file')

    parser1.add_argument('--update-thumb', '-t', action='store_true',
                         help='Force udpate of SED thumbnail image')

    args = parser1.parse_args()
    
    # collect the files
    if args.file is not None:
        files = [path for group in args.file for path in group]

    elif args.dir is not None:
        files = []
        for d in (path for group in args.dir for path in group):
            d.rstrip('/')
            files += glob.glob(os.path.abspath(d) + '/**/'
                               + args.subset[0] + '/*-rawphot.txt',
                               recursive=True)

    elif args.samples is not None:
        ids = []
        for s in args.samples:
            ids += db.sample_targets(s)
        files = []
        for id_ in ids:
            files += glob.glob(cfg.file['sdb_root']+'masters/'
                               + id_ + '/**/*-rawphot.txt')

    # Joint-package member files collapse to one package-level fit. Ordinary
    # rawphot files remain independent and keep the existing Result workflow.
    sources = discover_fit_sources(files)
    for source in sources:
        source = os.fspath(source)
        is_joint = os.path.isdir(source)
        if is_joint:
            lock_path = os.path.join(source, '.sdf_lock')
        else:
            lock_path = os.path.dirname(os.path.abspath(source)) \
                + '/.sdf_lock-' + os.path.basename(source)
        lock = filelock.FileLock(lock_path)
        try:
            with lock.acquire(timeout=0):
                print(source)

                if is_joint:
                    if args.www or args.update_www or \
                            args.dbwrite or args.update_db:
                        raise NotImplementedError(
                            'joint-fit web/database output will follow the '
                            'hierarchy-aware analysis layer'
                        )
                    fit_package(
                        source,
                        force=args.update_all,
                        nospec=args.no_spectra,
                    )
                    continue

                # evidence-sorted list of results
                results = fitting.fit_results(os.path.abspath(source),
                                              update_mn=args.update_all,
                                              update_an=args.update_analysis,
                                              update_json=args.update_json,
                                              update_thumb=args.update_thumb,
                                              nospec=args.no_spectra)
                if results is None:
                    continue

                if args.www or args.update_www:
                    www.www_all(results, update=args.update_www)

                # write best model to db
                if args.dbwrite or args.update_db:
                    db.write_all(results[0], update=args.update_db)

        except filelock.Timeout:
            pass

#        tr.print_diff()


def sdf_sample():
    """Generate HTML sample pages to browse database."""
    
    # inputs
    parser = argparse.ArgumentParser(description='Update sample www pages')
    parser.add_argument('--samples', '-s', nargs='+',
                        help='Restrict to target samples')
    parser.add_argument('--tables', '-t', action='store_true',
                        help='Update sample tables')
    parser.add_argument('--plots', '-p', action='store_true',
                        help='Update sample plots')
    parser.add_argument('--calibration', '-l', action='store_true',
                        help='Update calibration plots')
    parser.add_argument('--cleanup', '-c', action='store_true',
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
            plotting.calibration(sample=sample,
                                 min_wav=cfg.www['min_cal_wav'],
                                 max_wav=cfg.www['max_cal_wav'])

    if args.cleanup:
        print("Cleaning up")
        www.cleanup_sample_dirs()
        www.cleanup_calibration_dirs()


def sdf_cleanup():
    """Cleanup."""

    parser = argparse.ArgumentParser(description='Remove multinest files')
    parser.add_argument('--sample', '-s', nargs='+',
                        help='Sample(s) to clean up')
    args = parser.parse_args()

    if args.sample is not None:
        for s in args.sample:

            ids = db.sample_targets(s)
            fs = []
            for id_ in ids:
                fs += glob.glob(cfg.file['sdb_root'] + 'masters/'
                                + id_ + '/*/*mnest/*')

            removed = [os.unlink(f) for f in fs]
