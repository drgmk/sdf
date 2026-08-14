"""Shared public/operator page shell for retained SED fits."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode

import jinja2
from bokeh.resources import CDN, INLINE


def render_sed_page(*, title, plot_script, plot_div, targets, candidates,
                    observations, default_candidate, inline=True,
                    operator=False, available_fits=(), operator_api=None,
                    operator_token=None, operator_actor=None, queued_jobs=()):
    """Render one SED document, optionally extended with operator controls."""

    env = jinja2.Environment(
        autoescape=True,
        loader=jinja2.PackageLoader('sdf', package_path='www/templates'),
    )
    template = env.get_template('joint_sed_page.html')
    resources = INLINE if inline else CDN
    return template.render(
        title=title,
        plot_script=plot_script,
        plot_div=plot_div,
        targets=targets,
        candidates=candidates,
        observations=observations,
        default_candidate=default_candidate,
        operator=operator,
        available_fits=available_fits,
        operator_api=operator_api,
        operator_token=operator_token,
        operator_actor=operator_actor,
        queued_jobs=queued_jobs,
        bokeh_js=resources.render_js(),
        bokeh_css=resources.render_css(),
        creation_time=datetime.now(timezone.utc).strftime('%d/%m/%y %X'),
    )


def target_page_rows(dataset, package_directory):
    """Project physical contributors and their component-specific links."""

    package_directory = Path(package_directory)
    observation_by_target = {
        observation.target_key: observation
        for observation in dataset.observations
    }
    rows = []
    for target in dataset.targets:
        keywords = target.keywords
        ra = _number(keywords.get('raj2000'))
        dec = _number(keywords.get('dej2000'))
        parallax = _number(keywords.get('plx_value'))
        distance = (
            None if parallax is None or parallax <= 0 else 1000.0 / parallax
        )
        observation = observation_by_target.get(target.key)
        source = _target_source(package_directory, target, observation)
        rows.append({
            'key': target.key,
            'main_id': target.main_id,
            'sdbid': target.sdbid,
            'label': target.main_id or target.sdbid or target.key,
            'spectral_type': _clean(keywords.get('sp_type')),
            'distance_pc': distance,
            'ra_deg': ra,
            'dec_deg': dec,
            'input_file': (
                None if source is None else source.relative_to(
                    package_directory
                ).as_posix()
            ),
            'external_links': _external_links(target.main_id, ra, dec),
        })
    return rows


def candidate_page_rows(candidates, labels, default_candidate,
                        package_directory, dataset, displayed_ids=None):
    """Project retained fits and component-aware result links."""

    package_directory = Path(package_directory)
    rows = []
    displayed_ids = (
        {candidate.candidate_id for candidate in candidates}
        if displayed_ids is None else set(displayed_ids)
    )
    for index, candidate in enumerate(candidates, start=1):
        evidence = candidate.posterior.evidence
        analysis = candidate.metadata['analysis']
        components = []
        for component in candidate.metadata['components']:
            components.append({
                'owner': labels.get(
                    component['owner_key'], component['owner_key']
                ),
                'role': component['role'],
                'model': component['model_name'],
            })
        links = [{
            'label': 'result',
            'url': (candidate.path / 'result.json').relative_to(
                package_directory
            ).as_posix(),
        }, {
            'label': 'posterior',
            'url': (candidate.path / candidate.metadata['posterior']['file']).relative_to(
                package_directory
            ).as_posix(),
        }]
        for label, filename in (
            ('fit posteriors', 'corner.png'),
            ('derived distributions', 'distributions.png'),
        ):
            path = candidate.path / filename
            if path.is_file():
                links.append({
                    'label': label,
                    'url': path.relative_to(package_directory).as_posix(),
                })
        rows.append({
            'tab': index,
            'candidate_id': candidate.candidate_id,
            'reference': candidate.path.name,
            'display_name': candidate.metadata['display_name'],
            'evidence': evidence,
            'chi_squared': analysis.get('chi_squared'),
            'current': candidate.is_current(dataset),
            'default': candidate.candidate_id == default_candidate,
            'displayed': candidate.candidate_id in displayed_ids,
            'components': components,
            'links': links,
            'preview_url': '?' + urlencode({'fit': candidate.path.name}),
        })
    return rows


def add_operator_policy(targets, policies, registry, *, policy_store=None,
                        package_directory=None, actor=None):
    """Add editable policy state without changing the public projection."""

    policies = {policy.target_key: policy for policy in policies}
    for target in targets:
        policy = policies[target['key']]
        record = (
            None if policy_store is None else
            policy_store.policy(target['key'], package=package_directory)
        )
        if policy.origin in {'package', 'operator-package'}:
            save_scope = 'package'
        else:
            save_scope = 'global'
        target['policy'] = {
            'stellar_model': policy.stellar_models[0],
            'stellar_models': policy.stellar_models,
            'disk_model': policy.disk_models[0] if policy.disk_models else '',
            'disk_models': policy.disk_models,
            'max_disks': policy.max_disks,
            'search_mode': policy.search_mode,
            'reddening': policy.reddening,
            'origin': policy.origin,
            'save_scope': save_scope,
            'resettable': policy.origin in {'operator', 'operator-package'},
            'reason': (
                record.reason if record is not None else
                _suggested_reason(target['label'], policy)
            ),
        }
        target['stellar_options'] = registry.names('stellar')
        target['disk_options'] = registry.names('disk')
        target['operator_actor'] = actor or ''
    return targets


def _suggested_reason(label, policy):
    stellar = ' / '.join(policy.stellar_models)
    if policy.max_disks:
        disks = ' / '.join(policy.disk_models)
        return 'Use {} for {}; allow up to {} {} disk component(s)'.format(
            stellar, label, policy.max_disks, disks,
        )
    return 'Use {} for {}; do not fit a disk'.format(stellar, label)


def observation_page_rows(dataset, package_directory):
    """Project observed scopes without conflating them with model owners."""

    package_directory = Path(package_directory)
    target_labels = {
        target.key: target.main_id or target.sdbid or target.key
        for target in dataset.targets
    }
    return [{
        'observation_id': observation.observation_id,
        'contributors': [
            target_labels.get(key, key) for key in observation.contributor_keys
        ],
        'input_file': Path(observation.source).relative_to(
            package_directory
        ).as_posix(),
    } for observation in dataset.observations]


def _target_source(package_directory, target, observation):
    if observation is not None:
        return Path(observation.source)
    if target.sdbid is not None:
        candidate = package_directory / '{}-rawphot.txt'.format(target.sdbid)
        if candidate.is_file():
            return candidate
    return None


def _external_links(main_id, ra, dec):
    links = []
    if main_id:
        links.append({
            'label': 'SIMBAD',
            'url': 'https://simbad.cds.unistra.fr/simbad/sim-id?' + urlencode({
                'submit': 'submit id', 'Ident': main_id,
            }),
        })
    if ra is None or dec is None:
        return links
    coordinate = '{} {}'.format(ra, dec)
    comma_coordinate = '{},{}'.format(ra, dec)
    links.extend([
        {
            'label': 'CDS',
            'url': 'https://cdsportal.u-strasbg.fr/?' + urlencode({
                'target': coordinate,
            }),
        },
        {
            'label': 'CASSIS',
            'url': 'https://cassis.sirtf.com/atlas/cgi/radec.py?' + urlencode({
                'ra': ra, 'dec': dec, 'radius': 20,
            }),
        },
        {
            'label': 'Finder',
            'url': (
                'https://irsa.ipac.caltech.edu/applications/finderchart/'
                'servlet/api?' + urlencode({
                    'mode': 'getResult', 'locstr': comma_coordinate,
                })
            ),
        },
        {
            'label': 'Spitzer',
            'url': (
                'https://sha.ipac.caltech.edu/applications/Spitzer/SHA/?' +
                urlencode({
                    'api': 'search', 'searchoption': 'POSITION', 'sr': '180s',
                    'WorldPt': '{};{};EQ_J2000'.format(ra, dec),
                    'execute': 'true',
                })
            ),
        },
        {
            'label': 'MAST',
            'url': (
                'https://mast.stsci.edu/portal/Mashup/Clients/Mast/'
                'Portal.html?' + urlencode({'searchQuery': comma_coordinate})
            ),
        },
        {
            'label': 'ESASky',
            'url': 'https://sky.esa.int/?' + urlencode({
                'action': 'goto', 'fov': '0.25', 'cooframe': 'J2000',
                'sci': 'true', 'hips': 'AllWISE color', 'target': coordinate,
            }),
        },
    ])
    return links


def _clean(value):
    if value is None or str(value).strip() in ('', 'None'):
        return None
    return str(value).strip()


def _number(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value
