"""Run disposable old/new comparisons for real single-star rawphot files."""

from __future__ import annotations

import argparse
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile


HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]
SCIENCE_REGRESSION = HERE / 'science_regression.py'


def run_cases(inputs, *, model='phoenix_sol', legacy_revision='master'):
    """Fit every input in temporary storage and return compact comparisons."""

    inputs = tuple(Path(value).expanduser().resolve() for value in inputs)
    missing = tuple(path for path in inputs if not path.is_file())
    if missing:
        raise FileNotFoundError(
            'missing rawphot input(s): {}'.format(
                ', '.join(str(path) for path in missing)
            )
        )

    current_revision = _git_output('rev-parse', 'HEAD')
    current_dirty = bool(_git_output(
        'status', '--porcelain', '--untracked-files=no'
    ))
    legacy_commit = _git_output('rev-parse', legacy_revision)
    rows = []
    with tempfile.TemporaryDirectory(
            prefix='sdf-real-single-regression-') as temporary:
        temporary = Path(temporary)
        legacy_source = temporary / 'legacy-source'
        _extract_revision(legacy_revision, legacy_source)
        for index, source in enumerate(inputs, start=1):
            label = _main_id(source) or source.name.removesuffix(
                '-rawphot.txt'
            )
            print('[{}/{}] {}'.format(index, len(inputs), label), flush=True)
            case = temporary / 'case-{:02d}'.format(index)
            case.mkdir()
            clean_input = case / source.name
            _copy_without_package_markers(source, clean_input)

            legacy_input = case / 'legacy' / source.name
            legacy_input.parent.mkdir()
            legacy_input.write_bytes(clean_input.read_bytes())
            legacy_summary = case / 'legacy.json'
            _run_quietly(
                [
                    sys.executable, str(SCIENCE_REGRESSION), 'legacy',
                    str(legacy_input), '--model', model, '--force',
                    '--output', str(legacy_summary),
                ],
                env=_pythonpath_environment(legacy_source),
                label='legacy fit for {}'.format(label),
            )

            current_package = case / 'current'
            _run_quietly(
                [
                    sys.executable, '-c',
                    'from sdf.scripts import sdf_fit; sdf_fit()',
                    str(clean_input), '--output', str(current_package),
                    '--stellar', model, '--no-spectra', '--force',
                ],
                env=_pythonpath_environment(REPOSITORY),
                label='current fit for {}'.format(label),
            )
            current_summary = case / 'current.json'
            _run_quietly(
                [
                    sys.executable, str(SCIENCE_REGRESSION), 'candidate',
                    str(current_package), '--output', str(current_summary),
                ],
                env=_pythonpath_environment(REPOSITORY),
                label='current summary for {}'.format(label),
            )
            comparison_path = case / 'comparison.json'
            _run_quietly(
                [
                    sys.executable, str(SCIENCE_REGRESSION), 'compare',
                    str(legacy_summary), str(current_summary),
                    '--output', str(comparison_path),
                ],
                env=_pythonpath_environment(REPOSITORY),
                label='comparison for {}'.format(label),
            )
            comparison = json.loads(comparison_path.read_text())
            rows.append(_summary_row(label, source, comparison))

    return {
        'schema': 'sdf-real-single-regression-1',
        'legacy_revision': legacy_commit,
        'current_revision': current_revision,
        'current_dirty': current_dirty,
        'model': model,
        'cases': rows,
    }


def _summary_row(label, source, comparison):
    quantities = {
        row['quantity']: row for row in comparison['quantities']
    }
    selected = {
        name: quantities[name]['z']
        for name in (
            'stellar.Teff', 'stellar.logg', 'parameters.norm',
            'stellar.lstar', 'stellar.rstar',
        )
        if name in quantities
    }
    finite_z = tuple(abs(value) for value in selected.values()
                     if value is not None)
    return {
        'target': label,
        'input': str(source),
        'z': selected,
        'max_abs_z': max(finite_z) if finite_z else None,
        'evidence_delta': comparison['evidence']['difference'],
        'evidence_z': comparison['evidence']['z'],
        'chi_squared_delta': comparison['chi_squared']['difference'],
    }


def _extract_revision(revision, destination):
    archive = subprocess.run(
        ['git', 'archive', '--format=tar', revision],
        cwd=REPOSITORY, check=True, capture_output=True,
    ).stdout
    destination.mkdir(parents=True)
    with tarfile.open(fileobj=io.BytesIO(archive), mode='r:') as handle:
        handle.extractall(destination, filter='data')


def _git_output(*arguments):
    return subprocess.run(
        ['git', *arguments], cwd=REPOSITORY, check=True,
        capture_output=True, text=True,
    ).stdout.strip()


def _copy_without_package_markers(source, destination):
    """Make an unpackaged photometry-only copy of one exported rawphot file."""

    lines = []
    for line in source.read_text(encoding='utf-8').splitlines(keepends=True):
        if line.startswith('\\'):
            name, separator, _ = line[1:].partition('=')
            if separator and name.strip() in {
                'sdf_package', 'spectra_manifest'
            }:
                continue
        lines.append(line)
    destination.write_text(''.join(lines), encoding='utf-8')


def _main_id(path):
    for line in path.read_text(encoding='utf-8').splitlines():
        if not line.startswith('\\'):
            break
        name, separator, value = line[1:].partition('=')
        if separator and name.strip() == 'main_id':
            return value.strip().strip("'\"")
    return None


def _case_list(path, export_root):
    export_root = Path(export_root).expanduser().resolve()
    lines = Path(path).read_text(encoding='utf-8').splitlines()
    identifiers = tuple(
        line.strip() for line in lines
        if line.strip() and not line.lstrip().startswith('#')
    )
    return tuple(
        export_root / identifier / '{}-rawphot.txt'.format(identifier)
        for identifier in identifiers
    )


def _pythonpath_environment(source):
    environment = os.environ.copy()
    environment['PYTHONPATH'] = str(source)
    return environment


def _run_quietly(command, *, env, label):
    completed = subprocess.run(
        command, cwd=Path(tempfile.gettempdir()), env=env,
        capture_output=True, text=True,
    )
    if completed.returncode:
        detail = '\n'.join(
            (completed.stdout + '\n' + completed.stderr).splitlines()[-40:]
        )
        raise RuntimeError('{} failed:\n{}'.format(label, detail))


def _print_summary(report):
    print('\n{:<22} {:>8} {:>8} {:>8} {:>8} {:>8} {:>9}'.format(
        'target', 'Teff', 'logg', 'norm', 'Lstar', 'Rstar', 'max |z|'
    ))
    for case in report['cases']:
        values = case['z']
        print('{:<22} {:>8} {:>8} {:>8} {:>8} {:>8} {:>9}'.format(
            case['target'][:22],
            _z_text(values.get('stellar.Teff')),
            _z_text(values.get('stellar.logg')),
            _z_text(values.get('parameters.norm')),
            _z_text(values.get('stellar.lstar')),
            _z_text(values.get('stellar.rstar')),
            _z_text(case['max_abs_z']),
        ))
    print('\nValues are (current - legacy) / combined posterior sigma.')


def _z_text(value):
    return 'n/a' if value is None else '{:+.3f}'.format(value)


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Run disposable pre-branch/current fits for real single-star '
            'rawphot files and retain only a compact comparison summary.'
        )
    )
    parser.add_argument('rawphot', nargs='*')
    parser.add_argument('--model', default='phoenix_sol')
    parser.add_argument('--legacy-revision', default='master')
    parser.add_argument(
        '--export-root',
        help='resolve IDs from --case-list below this SDB export root',
    )
    parser.add_argument(
        '--case-list', default=str(HERE / 'real-single-cases.txt'),
        help='newline-separated SDB IDs used with --export-root',
    )
    parser.add_argument('--output', help='optional compact summary JSON')
    args = parser.parse_args()

    inputs = tuple(args.rawphot)
    if args.export_root is not None:
        if inputs:
            parser.error('use rawphot paths or --export-root, not both')
        inputs = _case_list(args.case_list, args.export_root)
    if not inputs:
        parser.error('provide rawphot paths or --export-root')

    report = run_cases(
        inputs, model=args.model,
        legacy_revision=args.legacy_revision,
    )
    _print_summary(report)
    if args.output is not None:
        Path(args.output).write_text(
            json.dumps(report, indent=2, sort_keys=True) + '\n',
            encoding='utf-8',
        )


if __name__ == '__main__':
    main()
