"""Summarize and compare legacy and retained-candidate stellar fits.

Run the ``legacy`` command with ``PYTHONPATH`` pointing at the pre-branch SDF
checkout. Run ``candidate`` with the current checkout. The resulting small JSON
files have one common schema and can be compared without loading old pickles.
"""

from __future__ import annotations

import argparse
import json
from math import hypot
from pathlib import Path


SCHEMA = 'sdf-stellar-science-regression-1'


def legacy_summary(rawphot, model, force=False):
    """Run and summarize the established single-file Result workflow."""

    from sdf.result import Result

    result = Result.get(
        str(Path(rawphot).expanduser().resolve()), (model,),
        update_mn=force, update_an=force, update_json=force,
        nospec=True,
    )
    if result is None or not hasattr(result, 'parameters'):
        raise RuntimeError('legacy fit produced no result')
    stats = result.analyzer.get_stats()
    stellar = next(
        value for value in result.main_results
        if 'Teff' in value
    )
    return {
        'schema': SCHEMA,
        'kind': 'legacy-result',
        'model': model,
        'parameters': {
            name: _symmetric(value, uncertainty)
            for name, value, uncertainty in zip(
                result.parameters, result.best_params,
                result.best_params_1sig,
            )
        },
        'stellar': {
            'Teff': _symmetric(stellar['Teff'], stellar['e_Teff']),
            'logg': _symmetric(stellar['logg'], stellar['e_logg']),
            'lstar': _asymmetric(
                stellar['lstar'], stellar['e_lstar_lo'],
                stellar['e_lstar_hi'],
            ),
            'rstar': _asymmetric(
                stellar['rstar'], stellar['e_rstar_lo'],
                stellar['e_rstar_hi'],
            ),
        },
        'evidence': float(stats['global evidence']),
        'evidence_error': float(stats['global evidence error']),
        'evidence_sampling_error': float(stats.get(
            'nested sampling global log-evidence error',
            stats['global evidence error'],
        )),
        'chi_squared': float(result.chisq),
    }


def candidate_summary(package, reference=None, owner=None):
    """Summarize one current retained candidate in the common schema."""

    from sdf.candidate_library import CandidateLibrary
    from sdf.fit_management import FitManager

    package = Path(package).expanduser().resolve()
    manager = FitManager(package)
    library = CandidateLibrary(package)
    if reference is None:
        try:
            candidate_id = library.read_selection().default
        except FileNotFoundError:
            candidates = library.candidates()
            if len(candidates) != 1:
                raise RuntimeError(
                    'select a candidate when the package does not have one '
                    'unambiguous default'
                )
            candidate_id = candidates[0].candidate_id
    else:
        candidate_id = manager.resolve(reference)
    result = library.find(candidate_id)
    components = result.metadata['analysis']['components']
    stellar_components = tuple(
        component for component in components
        if component['role'] == 'stellar'
    )
    if owner is None:
        if len(stellar_components) != 1:
            raise RuntimeError(
                'select --owner for a fit with multiple stellar components'
            )
        stellar_component = stellar_components[0]
    else:
        owner = ' '.join(owner.split())
        target_matches = {
            key for key, value in result.metadata['analysis']['targets'].items()
            if owner in {
                key,
                ' '.join(str(value.get('main_id') or '').split()),
                str(value.get('sdbid') or ''),
            }
        }
        matches = tuple(
            component for component in stellar_components
            if component['owner_key'] in target_matches
        )
        if len(matches) != 1:
            raise RuntimeError(
                'owner {!r} matched {} stellar components'.format(
                    owner, len(matches)
                )
            )
        stellar_component = matches[0]
    owner_key = stellar_component['owner_key']
    target = result.metadata['analysis']['targets'][owner_key]
    target_component = next(
        component for component in target['stellar_components']
        if component['component_id'] == stellar_component['component_id']
    )
    parameters = {
        name.rsplit('.', 1)[-1]: dict(values)
        for name, values in result.metadata['analysis']['parameters'].items()
        if name.startswith(stellar_component['component_id'] + '.')
    }
    return {
        'schema': SCHEMA,
        'kind': 'retained-candidate',
        'model': stellar_component['model_name'],
        'parameters': parameters,
        'stellar': {
            'Teff': dict(stellar_component['parameters']['Teff']),
            'logg': dict(stellar_component['parameters']['logg']),
            'lstar': dict(target_component['luminosity_lsun']),
            'rstar': dict(target_component['radius_rsun']),
        },
        'evidence': result.posterior.evidence,
        'evidence_error': result.posterior.evidence_error,
        'evidence_sampling_error': result.posterior.diagnostics.get(
            'stats', {}
        ).get(
            'nested sampling global log-evidence error',
            result.posterior.evidence_error,
        ),
        'chi_squared': result.metadata['analysis']['chi_squared'],
    }


def compare_summaries(first, second, include_fit_statistics=True):
    """Return normalized differences for compatible science summaries."""

    if first.get('schema') != SCHEMA or second.get('schema') != SCHEMA:
        raise ValueError('science summaries use an unsupported schema')
    if first['model'] != second['model']:
        raise ValueError(
            'cannot compare {} with {}'.format(
                first['model'], second['model']
            )
        )
    rows = []
    for section in ('parameters', 'stellar'):
        common = sorted(set(first[section]) & set(second[section]))
        for name in common:
            left = first[section][name]
            right = second[section][name]
            left_sigma = 0.5 * (left['minus'] + left['plus'])
            right_sigma = 0.5 * (right['minus'] + right['plus'])
            combined = hypot(left_sigma, right_sigma)
            difference = right['median'] - left['median']
            rows.append({
                'quantity': '{}.{}'.format(section, name),
                'first': left['median'],
                'second': right['median'],
                'difference': difference,
                'combined_sigma': combined,
                'z': None if combined == 0 else difference / combined,
            })
    comparison = {
        'schema': SCHEMA,
        'model': first['model'],
        'quantities': rows,
    }
    if include_fit_statistics:
        evidence_delta = second['evidence'] - first['evidence']
        evidence_sigma = hypot(
            first.get('evidence_sampling_error', first.get('evidence_error'))
            or 0.0,
            second.get('evidence_sampling_error', second.get('evidence_error'))
            or 0.0,
        )
        comparison['evidence'] = {
            'first': first['evidence'],
            'second': second['evidence'],
            'difference': evidence_delta,
            'combined_sigma': evidence_sigma,
            'error_kind': 'nested-sampling run uncertainty',
            'z': (
                None if evidence_sigma == 0 else
                evidence_delta / evidence_sigma
            ),
        }
        comparison['chi_squared'] = {
            'first': first['chi_squared'],
            'second': second['chi_squared'],
            'difference': second['chi_squared'] - first['chi_squared'],
        }
    return comparison


def _symmetric(median, uncertainty):
    return _asymmetric(median, uncertainty, uncertainty)


def _asymmetric(median, minus, plus):
    return {
        'median': float(median),
        'minus': float(minus),
        'plus': float(plus),
    }


def _write(payload, output):
    text = json.dumps(payload, indent=2, sort_keys=True) + '\n'
    if output is None:
        print(text, end='')
    else:
        Path(output).write_text(text, encoding='utf-8')


def _read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def _print_comparison(comparison):
    print('{:<24} {:>14} {:>14} {:>10}'.format(
        'quantity', 'first', 'second', 'delta/sigma'
    ))
    for row in comparison['quantities']:
        z = 'n/a' if row['z'] is None else '{:.3f}'.format(row['z'])
        print('{:<24} {:>14.6g} {:>14.6g} {:>10}'.format(
            row['quantity'], row['first'], row['second'], z
        ))
    if 'evidence' in comparison:
        evidence = comparison['evidence']
        z = 'n/a' if evidence['z'] is None else '{:.3f}'.format(evidence['z'])
        print('\nlog evidence: {:.6g} -> {:.6g} (delta {:.4g}, {} sigma)'.format(
            evidence['first'], evidence['second'], evidence['difference'], z,
        ))
        chi = comparison['chi_squared']
        print('chi squared:  {:.6g} -> {:.6g} (delta {:.4g})'.format(
            chi['first'], chi['second'], chi['difference']
        ))


def main():
    parser = argparse.ArgumentParser(
        description='Summarize and compare old/new SDF stellar fits.'
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    legacy = subparsers.add_parser('legacy')
    legacy.add_argument('rawphot')
    legacy.add_argument('--model', default='phoenix_sol')
    legacy.add_argument('--force', action='store_true')
    legacy.add_argument('--output')
    candidate = subparsers.add_parser('candidate')
    candidate.add_argument('package')
    candidate.add_argument('--fit')
    candidate.add_argument(
        '--owner', help='stellar owner key, SIMBAD name, or SDB ID',
    )
    candidate.add_argument('--output')
    compare = subparsers.add_parser('compare')
    compare.add_argument('first')
    compare.add_argument('second')
    compare.add_argument(
        '--parameters-only', action='store_true',
        help='omit evidence and chi-squared for different datasets',
    )
    compare.add_argument('--output')
    args = parser.parse_args()

    if args.command == 'legacy':
        _write(
            legacy_summary(args.rawphot, args.model, force=args.force),
            args.output,
        )
    elif args.command == 'candidate':
        _write(
            candidate_summary(
                args.package, reference=args.fit, owner=args.owner
            ),
            args.output,
        )
    else:
        comparison = compare_summaries(
            _read(args.first), _read(args.second),
            include_fit_statistics=not args.parameters_only,
        )
        if args.output is not None:
            _write(comparison, args.output)
        _print_comparison(comparison)


if __name__ == '__main__':
    main()
