import os
import glob
import argparse
from dataclasses import asdict
import json
from pathlib import Path

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
from sdf.fit_management import FitManager
from sdf.joint_plotting import write_joint_results_page
from sdf.result_index import ResultIndex
from sdf.operator_jobs import run_next_job, run_worker
from sdf.operator_state import OperatorState
from sdf.catalog_server import serve_catalog
from sdf.sample_products import write_sample_products


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
                    if args.dbwrite or args.update_db:
                        raise NotImplementedError(
                            'joint-fit database output requires the result '
                            'index and is not implemented yet'
                        )
                    fit_package(
                        source,
                        force=args.update_all,
                        nospec=args.no_spectra,
                    )
                    write_joint_results_page(source)
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
    parser.add_argument(
        '--index',
        help='Use retained joint-fit index instead of legacy result tables',
    )
    parser.add_argument(
        '--output-root',
        help='Sample output root (default: configured www_root/samples)',
    )
    parser.add_argument(
        '--sed-base-url',
        help='Base URL containing package pages (default: configured sdb_url)',
    )
    args = parser.parse_args()

    if args.index and (args.tables or args.plots):
        sample_names = args.samples or db.get_samples()
        output_root = Path(
            args.output_root or (cfg.file['www_root'] + 'samples')
        )
        index = ResultIndex(args.index)
        for sample in sample_names:
            print("Updating indexed sample {}".format(sample))
            write_sample_products(
                index, sample, db.sample_targets(sample),
                output_root / sample,
                table=args.tables, plots=args.plots,
                sed_base_url=args.sed_base_url or cfg.www['sdb_url'],
            )
        return

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


def sdf_index():
    """Build or verify the retained-fit query index."""

    parser = argparse.ArgumentParser(
        description='Manage the rebuildable SDF result index'
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    for command in ('rebuild', 'verify'):
        subparser = subparsers.add_parser(command)
        subparser.add_argument('root', help='Root containing exported packages')
        subparser.add_argument(
            '--database', help='Index path (default: ROOT/sdf-results.sqlite)'
        )
    args = parser.parse_args()
    root = os.path.abspath(args.root)
    database = args.database or os.path.join(root, 'sdf-results.sqlite')
    index = ResultIndex(database)
    if args.command == 'rebuild':
        print(json.dumps(asdict(index.rebuild(root)), sort_keys=True))
        return
    problems = index.verify(root)
    if problems:
        parser.exit(1, '\n'.join(problems) + '\n')
    print(json.dumps(asdict(index.summary()), sort_keys=True))


def sdf_fits():
    """Inspect and select retained package fits."""

    parser = argparse.ArgumentParser(
        description='Manage retained fits and the SDF export catalog'
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    list_parser = subparsers.add_parser('list')
    list_parser.add_argument('package')
    list_parser.add_argument('--json', action='store_true')
    select_parser = subparsers.add_parser('select')
    select_parser.add_argument('package')
    select_parser.add_argument('--default')
    select_parser.add_argument('--display', nargs='+')
    render_parser = subparsers.add_parser('render')
    render_parser.add_argument('package')
    render_parser.add_argument('--output')
    render_parser.add_argument(
        '--operator', action='store_true',
        help='Render the shared page with the local operator drawer',
    )
    render_parser.add_argument(
        '--fit', action='append', dest='fits',
        help='Render only this retained fit (repeat for several tabs)',
    )
    fit_parser = subparsers.add_parser('fit')
    fit_parser.add_argument('package')
    fit_parser.add_argument(
        '--stellar', action='append', required=True, type=_model_assignment,
        metavar='TARGET=MODEL',
    )
    fit_parser.add_argument(
        '--disk', action='append', default=[], type=_model_assignment,
        metavar='TARGET=MODEL',
    )
    fit_parser.add_argument('--force', action='store_true')
    fit_parser.add_argument('--no-spectra', action='store_true')
    refit_parser = subparsers.add_parser('refit')
    refit_parser.add_argument('package')
    refit_parser.add_argument('fit')
    refit_parser.add_argument('--no-spectra', action='store_true')
    serve_parser = subparsers.add_parser('serve')
    serve_parser.add_argument(
        'package', metavar='EXPORT_ROOT',
        help='SDB export root containing all fitting packages',
    )
    serve_parser.add_argument('--state')
    serve_parser.add_argument('--actor')
    serve_parser.add_argument('--host', default='127.0.0.1')
    serve_parser.add_argument('--port', type=int, default=8765)
    serve_parser.add_argument('--index')
    serve_parser.add_argument('--rebuild-index', action='store_true')
    serve_parser.add_argument('--no-worker', action='store_true')
    serve_parser.add_argument('--no-watch', action='store_true')
    serve_parser.add_argument('--autofit-existing', action='store_true')
    serve_parser.add_argument('--poll-seconds', type=float, default=5.0)
    args = parser.parse_args()

    if args.command == 'serve':
        serve_catalog(
            args.package, index_path=args.index, state_path=args.state,
            actor=args.actor, rebuild=args.rebuild_index,
            host=args.host, port=args.port,
            worker=not args.no_worker, watch=not args.no_watch,
            autofit_existing=args.autofit_existing,
            poll_seconds=args.poll_seconds,
        )
        return

    manager = FitManager(args.package)

    if args.command == 'list':
        available = manager.available()
        if args.json:
            print(json.dumps([asdict(fit) for fit in available], indent=2))
            return
        for fit in available:
            flags = []
            if fit.is_default:
                flags.append('default')
            if fit.display_position is not None:
                flags.append('display {}'.format(fit.display_position + 1))
            evidence = (
                'n/a' if fit.evidence is None else '{:.3f}'.format(fit.evidence)
            )
            print('{}\t{}\t{}\t{}'.format(
                fit.reference, evidence, ', '.join(flags) or '-',
                fit.display_name,
            ))
        return

    if args.command == 'select':
        selection = manager.select(
            default=args.default, displayed=args.display
        )
        manager.render()
        print(json.dumps(asdict(selection), sort_keys=True))
        return

    if args.command == 'fit':
        result = manager.fit_exact(
            args.stellar, args.disk, force=args.force,
            nospec=args.no_spectra,
        )
        manager.render()
        print(result.path.name)
        return

    if args.command == 'refit':
        result = manager.refit(args.fit, nospec=args.no_spectra)
        manager.render()
        print(result.path.name)
        return

    print(manager.render(
        args.output, fits=args.fits, operator=args.operator,
    ))


def sdf_jobs():
    """Inspect or run the asynchronous operator fit queue."""

    parser = argparse.ArgumentParser(description='Manage queued SDF fit jobs')
    parser.add_argument('--state', required=True, help='Operator-state SQLite')
    subparsers = parser.add_subparsers(dest='command', required=True)
    list_parser = subparsers.add_parser('list')
    list_parser.add_argument('--status')
    run_parser = subparsers.add_parser('run')
    run_parser.add_argument('--once', action='store_true')
    run_parser.add_argument('--poll-seconds', type=float, default=5.0)
    args = parser.parse_args()
    state = OperatorState(args.state).initialize()
    if args.command == 'list':
        print(json.dumps([
            asdict(job) for job in state.jobs(status=args.status)
        ], indent=2))
        return
    if args.once:
        job = run_next_job(state, raise_errors=True)
        print('no queued jobs' if job is None else 'completed {}'.format(job.id))
        return
    run_worker(state, poll_seconds=args.poll_seconds)


def _model_assignment(value):
    owner, separator, model_name = value.rpartition('=')
    if not separator or not owner.strip() or not model_name.strip():
        raise argparse.ArgumentTypeError(
            'model assignment must be TARGET=MODEL'
        )
    return owner.strip(), model_name.strip()


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
