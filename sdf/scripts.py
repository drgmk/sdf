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
from sdf.coordinator import fit_package, reprioritize_package
from sdf.fit_management import FitManager
from sdf.joint_plotting import write_joint_results_page
from sdf.result_index import ResultIndex, discover_packages
from sdf.operator_jobs import run_next_job, run_worker
from sdf.operator_state import OperatorState
from sdf.catalog_server import serve_catalog
from sdf.catalog_auth import read_session_secret
from sdf.candidate_library import CandidateLibrary
from sdf.sample_products import write_sample_products
from sdf.sdb_import import SdbImportConfig


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
        description=(
            'Build or check the disposable SQLite index of retained SDF fits.\n\n'
            'The index makes catalog and sample queries fast by projecting '
            'package-local fit results below an SDB export root. It contains '
            'no administrative policy or job state and can always be rebuilt.'
        ),
        epilog=(
            'Examples:\n'
            '  sdf-index rebuild /data/sdb/export\n'
            '  sdf-index verify /data/sdb/export\n'
            '  sdf-index rebuild /data/sdb/export --database /data/index.sqlite'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    commands = (
        (
            'rebuild',
            'Rebuild the index from current package-local retained fits.',
            'Existing rebuildable index contents are replaced. Retained fit '
            'artifacts and the separate administrative-state database are '
            'unchanged.',
        ),
        (
            'verify',
            'Check the index against its export root and retained fits.',
            'Print the current index summary when it is valid; otherwise print '
            'each problem and exit with status 1.',
        ),
    )
    for command, summary, detail in commands:
        subparser = subparsers.add_parser(
            command,
            help=summary,
            description=summary + '\n\n' + detail,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        subparser.add_argument(
            'root', metavar='EXPORT_ROOT',
            help='SDB export root containing fitting-package directories',
        )
        subparser.add_argument(
            '--database', metavar='PATH',
            help='index SQLite path (default: EXPORT_ROOT/sdf-results.sqlite)',
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
        description=(
            'Fit, inspect, present, or serve retained SDF fitting packages.\n\n'
            'Package commands operate on one directory exported by SDB. The '
            'serve command operates on the complete export root and provides '
            'the catalog-wide administrative UI, fit queue, and export watcher.'
        ),
        epilog=(
            'Examples:\n'
            '  sdf-fits serve /data/sdb/export --actor Grant\n'
            '  sdf-fits list /data/sdb/export/sdbid-v3-...\n'
            '  sdf-fits render /data/sdb/export/sdbid-v3-... --operator\n'
            '  sdf-fits select PACKAGE --default FIT --display FIT OTHER_FIT\n\n'
            'Run `sdf-fits COMMAND --help` for command-specific details.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    list_parser = _sdf_fits_command(
        subparsers, 'list',
        'List all retained candidate fits for one package.',
        'Shows the fit reference, evidence, default/display state, and model '
        'summary. Fit references may be full IDs, storage names, or '
        'unambiguous prefixes.',
    )
    _package_argument(list_parser)
    list_parser.add_argument(
        '--json', action='store_true',
        help='write the complete candidate records as indented JSON',
    )
    select_parser = _sdf_fits_command(
        subparsers, 'select',
        'Choose the default fit and the fits shown on the public page.',
        'Updates fits/selection.json, regenerates the package page, and '
        'refreshes a discoverable result index. This does not run a fit.',
    )
    _package_argument(select_parser)
    select_parser.add_argument(
        '--default', metavar='FIT',
        help='fit reference to make the default; retain the current default if omitted',
    )
    select_parser.add_argument(
        '--display', nargs='+', metavar='FIT',
        help='replace the displayed fit-tab list; retain the current list if omitted',
    )
    render_parser = _sdf_fits_command(
        subparsers, 'render',
        'Regenerate or preview the HTML SED page for one package.',
        'Uses the retained compact fit artifacts; MultiNest is not rerun. By '
        'default the package selection controls the tabs.',
    )
    _package_argument(render_parser)
    render_parser.add_argument(
        '--output', metavar='PATH',
        help='HTML output path (default: PACKAGE/index.html or operator.html)',
    )
    render_parser.add_argument(
        '--operator', action='store_true',
        help='Render the shared page with the local operator drawer',
    )
    render_parser.add_argument(
        '--fit', action='append', dest='fits',
        metavar='FIT',
        help='render only this retained fit; repeat for several tabs',
    )
    fit_parser = _sdf_fits_command(
        subparsers, 'fit',
        'Synchronously fit one explicit model combination.',
        'Assign exactly one stellar model to every physical contributor. '
        'Repeat --disk to add multiple disk components to an owner. An '
        'identical current candidate is reused unless --force is supplied.',
        examples=(
            "sdf-fits fit PACKAGE --stellar 'HD 113766A=phoenix_sol' "
            "--stellar 'HD 113766B=phoenix_sol'",
            "sdf-fits fit PACKAGE --stellar 'HD 113766A=phoenix_sol' "
            "--disk 'HD 113766A=modbb_disk_r'",
        ),
    )
    _package_argument(fit_parser)
    fit_parser.add_argument(
        '--stellar', action='append', required=True, type=_model_assignment,
        metavar='TARGET=MODEL',
        help='stellar model owned by TARGET; repeat once per physical contributor',
    )
    fit_parser.add_argument(
        '--disk', action='append', default=[], type=_model_assignment,
        metavar='TARGET=MODEL',
        help='disk model owned by TARGET; repeat for each requested disk component',
    )
    fit_parser.add_argument(
        '--force', action='store_true',
        help='rerun the sampler even when this exact current candidate exists',
    )
    fit_parser.add_argument(
        '--no-spectra', action='store_true',
        help='exclude spectra from the fitted dataset',
    )
    refit_parser = _sdf_fits_command(
        subparsers, 'refit',
        'Rerun the sampler for one retained candidate.',
        'Reuses the retained candidate model definition, replaces that exact '
        'candidate artifact, and regenerates the package page.',
    )
    _package_argument(refit_parser)
    refit_parser.add_argument(
        'fit', metavar='FIT',
        help='full candidate ID, storage name, display name, or unambiguous prefix',
    )
    refit_parser.add_argument(
        '--no-spectra', action='store_true',
        help='exclude spectra from the fitted dataset',
    )
    reprioritize_parser = _sdf_fits_command(
        subparsers, 'reprioritize',
        'Reapply the evidence-priority rules to retained fits.',
        'Only selection metadata and the disposable result index are changed; '
        'this command never invokes a sampler. Pass one PACKAGE, or use '
        '--all with an EXPORT_ROOT to update every fitted package.',
        examples=(
            'sdf-fits reprioritize PACKAGE',
            'sdf-fits reprioritize /data/sdb/export --all',
        ),
    )
    reprioritize_parser.add_argument(
        'package', metavar='PACKAGE_OR_EXPORT_ROOT',
        help='one fitting package, or an export root when --all is used',
    )
    reprioritize_parser.add_argument(
        '--all', action='store_true',
        help='reprioritize every fitted package below EXPORT_ROOT',
    )
    reprioritize_parser.add_argument(
        '--state', metavar='PATH',
        help=(
            'administrative policy SQLite path (default with --all: '
            'EXPORT_ROOT/sdf-operator.sqlite)'
        ),
    )
    reprioritize_parser.add_argument(
        '--threshold', type=float, metavar='LOG_EVIDENCE',
        help='significance threshold (default: configured ev_threshold)',
    )
    refresh_parser = _sdf_fits_command(
        subparsers, 'refresh-artifacts',
        'Refresh retained posterior and SED data files.',
        'Splits older combined NPZ files and makes their fitted, derived, '
        'and model component axes self-describing. This operation never '
        'invokes a sampler. Pass one PACKAGE, or use --all with an '
        'EXPORT_ROOT.',
        examples=(
            'sdf-fits refresh-artifacts PACKAGE',
            'sdf-fits refresh-artifacts /data/sdb/export --all',
        ),
    )
    refresh_parser.add_argument(
        'package', metavar='PACKAGE_OR_EXPORT_ROOT',
        help='one fitting package, or an export root when --all is used',
    )
    refresh_parser.add_argument(
        '--all', action='store_true',
        help='refresh every fitted package below EXPORT_ROOT',
    )
    serve_parser = _sdf_fits_command(
        subparsers, 'serve',
        'Serve the catalog-wide SDF administrative UI and fitting loop.',
        'Starts a localhost web UI over every package in EXPORT_ROOT. By '
        'default it also supervises one fit worker and watches for new or '
        'changed SDB exports. Existing unfit packages are not queued unless '
        '--autofit-existing is supplied. Use --workers to fit independent '
        'systems concurrently; one package is never run by two workers.',
        examples=(
            'sdf-fits serve /data/sdb/export --actor Grant',
            'sdf-fits serve /data/sdb/export --workers 4 --autofit-existing',
            'sdf-fits serve /data/sdb/export --auth-user-header X-SDF-User '
            '--auth-database /srv/sdb/sdb.sqlite '
            '--auth-session-secret-file /srv/sdb/session-secret',
            'sdf-fits serve /data/sdb/export --no-worker --no-watch',
        ),
    )
    serve_parser.add_argument(
        'package', metavar='EXPORT_ROOT',
        help='SDB export root containing all fitting packages',
    )
    serve_parser.add_argument(
        '--state', metavar='PATH',
        help='durable policy/job SQLite path (default: EXPORT_ROOT/sdf-operator.sqlite)',
    )
    serve_parser.add_argument(
        '--actor',
        help=(
            'local-development operator; proxy-authenticated use records '
            'the authenticated username'
        ),
    )
    serve_parser.add_argument(
        '--host', default='127.0.0.1',
        help='bind host (default: 127.0.0.1)',
    )
    serve_parser.add_argument(
        '--port', type=int, default=8765,
        help='TCP port (default: 8765)',
    )
    serve_parser.add_argument(
        '--index', metavar='PATH',
        help='rebuildable result-index SQLite path (default: EXPORT_ROOT/sdf-results.sqlite)',
    )
    serve_parser.add_argument(
        '--rebuild-index', action='store_true',
        help='rebuild the result index from retained fits before serving',
    )
    serve_parser.add_argument(
        '--no-worker', action='store_true',
        help='do not start the supervised asynchronous fit worker',
    )
    serve_parser.add_argument(
        '--workers', type=_positive_int, default=1, metavar='N',
        help='fit up to N independent systems concurrently (default: 1)',
    )
    serve_parser.add_argument(
        '--no-watch', action='store_true',
        help='do not watch the export root for new or changed packages',
    )
    serve_parser.add_argument(
        '--autofit-existing', action='store_true',
        help='queue the existing unfit backlog when the export watcher starts',
    )
    serve_parser.add_argument(
        '--poll-seconds', type=float, default=5.0, metavar='SECONDS',
        help='worker and export-watcher polling interval (default: 5)',
    )
    serve_parser.add_argument(
        '--auth-user-header', metavar='HEADER',
        help='trust an authenticated username from this reverse-proxy header',
    )
    serve_parser.add_argument(
        '--auth-database', metavar='PATH',
        help=(
            'SDB SQLite database containing users and sample access policy; '
            'defaults to --sdb-database'
        ),
    )
    serve_parser.add_argument(
        '--auth-session-secret-file', metavar='PATH',
        help='shared file containing at least 32 bytes for signed login sessions',
    )
    serve_parser.add_argument(
        '--auth-cookie-secure', action='store_true',
        help='mark login-session cookies HTTPS-only (recommended in production)',
    )
    serve_parser.add_argument(
        '--review-base-url', metavar='URL',
        help='SDB review mount URL used for target-level Review data links',
    )
    serve_parser.add_argument(
        '--access-log', metavar='PATH',
        help=(
            'Apache access log used for aggregate admin activity statistics; '
            'client addresses are never displayed or copied'
        ),
    )
    serve_parser.add_argument(
        '--activity-timezone', default='UTC', metavar='ZONE',
        help=(
            'IANA timezone used for activity dates and displayed times '
            '(default: UTC; for example Pacific/Auckland)'
        ),
    )
    serve_parser.add_argument(
        '--sdb-database', metavar='PATH',
        help=(
            'enable signed-in target import using this authoritative SDB '
            'SQLite database'
        ),
    )
    serve_parser.add_argument(
        '--sdb-executable', default='sdb', metavar='PATH',
        help='SDB CLI executable used for target import (default: sdb on PATH)',
    )
    serve_parser.add_argument(
        '--sdb-config', metavar='PATH',
        help=(
            'SDB TOML configuration used during import, provider refresh, '
            'and mirror-aware SIMBAD/VizieR browser links'
        ),
    )
    serve_parser.add_argument(
        '--sdb-reference-database', metavar='PATH',
        help='SDB reference-snapshot SQLite database used during import',
    )
    serve_parser.add_argument(
        '--sdb-cache-database', metavar='PATH',
        help='SDB provider cache SQLite database used during import',
    )
    serve_parser.add_argument(
        '--sdb-import-timeout', type=float, default=1800.0,
        metavar='SECONDS',
        help='maximum time for each SDB command during import (default: 1800)',
    )
    serve_parser.add_argument(
        '--sdb-sample-import-timeout', type=float, default=43200.0,
        metavar='SECONDS',
        help=(
            'final wall-clock guard for a resumable sample import '
            '(default: 43200)'
        ),
    )
    args = parser.parse_args()

    if args.command == 'serve':
        import_options = (
            args.sdb_config, args.sdb_reference_database,
            args.sdb_cache_database,
        )
        if args.sdb_database is None and any(import_options):
            parser.error('--sdb-database is required with SDB import options')
        if args.sdb_import_timeout <= 0:
            parser.error('--sdb-import-timeout must be positive')
        if args.sdb_sample_import_timeout <= 0:
            parser.error('--sdb-sample-import-timeout must be positive')
        auth_database = args.auth_database or args.sdb_database
        if args.auth_user_header is not None and auth_database is None:
            parser.error(
                '--auth-database or --sdb-database is required with '
                '--auth-user-header'
            )
        auth_session_secret = (
            None if args.auth_session_secret_file is None
            else read_session_secret(args.auth_session_secret_file)
        )
        sdb_import_config = None
        if args.sdb_database is not None:
            sdb_import_config = SdbImportConfig.create(
                args.sdb_database, executable=args.sdb_executable,
                config=args.sdb_config,
                reference_database=args.sdb_reference_database,
                cache_database=args.sdb_cache_database,
                timeout_seconds=args.sdb_import_timeout,
                sample_timeout_seconds=args.sdb_sample_import_timeout,
            )
        serve_catalog(
            args.package, index_path=args.index, state_path=args.state,
            actor=args.actor, rebuild=args.rebuild_index,
            host=args.host, port=args.port,
            worker=not args.no_worker, workers=args.workers,
            watch=not args.no_watch,
            autofit_existing=args.autofit_existing,
            poll_seconds=args.poll_seconds,
            auth_user_header=args.auth_user_header,
            auth_database=auth_database,
            auth_session_secret=auth_session_secret,
            auth_cookie_secure=args.auth_cookie_secure,
            sdb_import_config=sdb_import_config,
            review_base_url=args.review_base_url,
            access_log=args.access_log,
            activity_timezone=args.activity_timezone,
        )
        return

    if args.command == 'reprioritize':
        root = Path(args.package).expanduser().resolve()
        packages = discover_packages(root) if args.all else (root,)
        state_path = args.state
        if state_path is None and args.all:
            default_state = root / 'sdf-operator.sqlite'
            state_path = default_state if default_state.is_file() else None
        state = (
            OperatorState(state_path).initialize()
            if state_path is not None else None
        )
        changed = 0
        skipped = 0
        for package in packages:
            selection = reprioritize_package(
                package, threshold=args.threshold, policy_store=state,
                result_index=False,
            )
            if selection is None:
                skipped += 1
            else:
                changed += 1
        index_path = root / 'sdf-results.sqlite' if args.all else None
        if args.all and (changed or index_path.is_file()):
            ResultIndex(index_path).rebuild(root)
        print(json.dumps({
            'reprioritized': changed, 'skipped': skipped,
            'packages': len(packages),
        }, sort_keys=True))
        return

    if args.command == 'refresh-artifacts':
        root = Path(args.package).expanduser().resolve()
        packages = discover_packages(root) if args.all else (root,)
        refreshed = sum(
            CandidateLibrary(package).refresh_artifacts()
            for package in packages
        )
        print(json.dumps({
            'refreshed': refreshed, 'packages': len(packages),
        }, sort_keys=True))
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
            if fit.is_recommended:
                flags.append('recommended')
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
    """Inspect or run the asynchronous administrative fit queue."""

    parser = argparse.ArgumentParser(
        description=(
            'Inspect or run the asynchronous SDF administrative fit queue.\n\n'
            'This is mainly useful with `sdf-fits serve --no-worker`, for '
            'headless workers, and for queue diagnostics. Only one worker '
            'should consume a given administrative-state database.'
        ),
        epilog=(
            'Examples:\n'
            '  sdf-jobs --state /data/sdb/export/sdf-operator.sqlite list\n'
            '  sdf-jobs --state /data/sdb/export/sdf-operator.sqlite list --status queued\n'
            '  sdf-jobs --state /data/sdb/export/sdf-operator.sqlite run --once\n'
            '  sdf-jobs --state /data/sdb/export/sdf-operator.sqlite run'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--state', required=True, metavar='PATH',
        help='durable administrative policy and fit-job SQLite database',
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    list_parser = subparsers.add_parser(
        'list',
        help='List queued and completed fit jobs as JSON.',
        description=(
            'List fit jobs newest first as indented JSON. Filter by an exact '
            'status when diagnosing the queue.'
        ),
    )
    list_parser.add_argument(
        '--status', choices=('queued', 'running', 'completed', 'failed'),
        help='show only jobs with this status',
    )
    run_parser = subparsers.add_parser(
        'run',
        help='Run fit jobs from the queue.',
        description=(
            'Claim and run queued jobs serially. On startup, jobs left running '
            'by an interrupted sole worker are returned to the queue.'
        ),
    )
    run_parser.add_argument(
        '--once', action='store_true',
        help='run at most one queued job and exit',
    )
    run_parser.add_argument(
        '--poll-seconds', type=float, default=5.0, metavar='SECONDS',
        help='queue polling interval for the continuous worker (default: 5)',
    )
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


def _positive_int(value):
    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError('must be a positive integer') from error
    if parsed < 1:
        raise argparse.ArgumentTypeError('must be a positive integer')
    return parsed


def _sdf_fits_command(subparsers, name, summary, detail, examples=()):
    epilog = None
    if examples:
        epilog = 'Examples:\n  ' + '\n  '.join(examples)
    return subparsers.add_parser(
        name,
        help=summary,
        description=summary + '\n\n' + detail,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )


def _package_argument(parser):
    parser.add_argument(
        'package', metavar='PACKAGE',
        help='fitting-package directory containing exported rawphot inputs',
    )


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
