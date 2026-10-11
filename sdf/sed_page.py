"""Shared public/operator page shell for retained SED fits."""

from __future__ import annotations

from datetime import datetime, timezone
import math
from pathlib import Path
import re
from urllib.parse import urlencode

import jinja2
from bokeh.resources import CDN, INLINE

from . import utils


def render_sed_page(*, title, plot_script, plot_div, targets, candidates,
                    observations, systems, default_candidate, inline=True,
                    operator=False, available_fits=(), operator_api=None,
                    operator_token=None, operator_actor=None,
                    operator_actor_locked=False, queued_jobs=(),
                    navigation=None):
    """Render one SED document, optionally extended with operator controls."""

    env = jinja2.Environment(
        autoescape=True,
        loader=jinja2.PackageLoader('sdf', package_path='www/templates'),
    )
    template = env.get_template('joint_sed_page.html')
    resources = INLINE if inline else CDN
    navigation = dict(navigation or {})
    navigation.setdefault(
        'search_action', '/admin/search' if operator else '/search'
    )
    default_fit = next((
        candidate for candidate in candidates if candidate['default']
    ), candidates[0] if candidates else None)
    return template.render(
        title=title,
        plot_script=plot_script,
        plot_div=plot_div,
        targets=targets,
        candidates=candidates,
        observations=observations,
        systems=systems,
        default_candidate=default_candidate,
        default_fit=default_fit,
        operator=operator,
        available_fits=available_fits,
        operator_api=operator_api,
        operator_token=operator_token,
        operator_actor=operator_actor,
        operator_actor_locked=operator_actor_locked,
        queued_jobs=queued_jobs,
        bokeh_js=resources.render_js(),
        bokeh_css=resources.render_css(),
        creation_time=datetime.now(timezone.utc).strftime('%d/%m/%y %X'),
        admin=operator,
        **navigation,
    )


def target_page_rows(dataset, package_directory):
    """Project physical contributors and their component-specific links."""

    package_directory = Path(package_directory)
    observation_by_target = {
        observation.target_key: observation
        for observation in dataset.observations
        if observation.kind == 'photometry'
    }
    system_positions = [
        (_number(target.keywords.get('raj2000')),
         _number(target.keywords.get('dej2000')))
        for target in dataset.targets
    ]
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
        search_radius = _search_radius_arcsec(ra, dec, system_positions)
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
            'external_links': _external_links(
                target.main_id, ra, dec, radius_arcsec=search_radius,
            ),
        })
    return rows


def candidate_page_rows(candidates, labels, default_candidate,
                        package_directory, dataset, displayed_ids=None,
                        component_labels=None, owner_colours=None):
    """Project retained fits and component-aware result links."""

    package_directory = Path(package_directory)
    rows = []
    displayed_ids = (
        {candidate.candidate_id for candidate in candidates}
        if displayed_ids is None else set(displayed_ids)
    )
    component_labels = component_labels or {}
    owner_colours = owner_colours or {}
    for index, candidate in enumerate(candidates, start=1):
        evidence = candidate.posterior.evidence
        analysis = candidate.metadata['analysis']
        components = []
        for component in candidate.metadata['components']:
            components.append({
                'owner': labels.get(
                    component['owner_key'], component['owner_key']
                ),
                'owner_key': component['owner_key'],
                'component': component_labels.get(component['owner_key']),
                'colour': owner_colours.get(
                    component['owner_key'], '#64748b'
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
        shorthand = _fit_shorthand(components)
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
            'shorthand': shorthand,
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
    targets_by_key = {target.key: target for target in dataset.targets}
    target_positions = {
        target.key: (
            _number(target.keywords.get('raj2000')),
            _number(target.keywords.get('dej2000')),
        )
        for target in dataset.targets
    }
    rows = []
    for observation in dataset.observations:
        if observation.kind != 'photometry':
            continue
        target = targets_by_key.get(observation.target_key)
        if target is not None:
            keywords = target.keywords
        else:
            try:
                keywords = utils.get_sdb_keywords(observation.source)
            except (OSError, ValueError):
                keywords = {}
        main_id = _clean(keywords.get('main_id'))
        ra = _number(keywords.get('raj2000'))
        dec = _number(keywords.get('dej2000'))
        parallax = _number(keywords.get('plx_value'))
        search_radius = _search_radius_arcsec(
            ra, dec, target_positions.values()
        )
        rows.append({
            'observation_id': observation.observation_id,
            'target_key': observation.target_key,
            'contributor_keys': list(observation.contributor_keys),
            'contributors': [
                target_labels.get(key, key)
                for key in observation.contributor_keys
            ],
            'label': main_id or observation.observation_id,
            'ra_deg': ra,
            'dec_deg': dec,
            'search_radius_arcsec': search_radius,
            'spectral_type': _clean(keywords.get('sp_type')),
            'distance_pc': (
                None if parallax is None or parallax <= 0 else
                1000.0 / parallax
            ),
            'input_file': Path(observation.source).relative_to(
                package_directory
            ).as_posix(),
            'external_links': _external_links(
                main_id, ra, dec, radius_arcsec=search_radius,
            ),
        })
    return rows


def system_page_rows(targets, observations):
    """Build a compact hierarchy projection and owner component labels."""

    targets_by_key = {target['key']: target for target in targets}
    component_labels = {}
    systems = []
    grouped = set()
    for observation in observations:
        contributor_keys = observation['contributor_keys']
        if len(contributor_keys) < 2:
            continue
        children = [
            targets_by_key[key] for key in contributor_keys
            if key in targets_by_key
        ]
        if len(children) < 2:
            continue
        child_components = _component_tokens(
            observation['label'], [child['label'] for child in children]
        )
        child_rows = []
        for child, component in zip(children, child_components):
            child = dict(child)
            child['component'] = component
            component_labels[child['key']] = component
            child_rows.append(child)
            grouped.add(child['key'])
        root_component = ''.join(child_components)
        if not all(re.fullmatch(r'[A-Z]+', value or '')
                   for value in child_components):
            root_component = 'system'
        systems.append({
            'component': root_component,
            'label': observation['label'],
            'spectral_type': observation['spectral_type'],
            'distance_pc': observation['distance_pc'],
            'input_file': observation['input_file'],
            'external_links': observation['external_links'],
            'children': child_rows,
        })

    for target in targets:
        if target['key'] in grouped:
            continue
        row = dict(target)
        row['component'] = None
        systems.append(row | {'children': []})

    return systems, component_labels


def _component_tokens(parent, children):
    parent_key = re.sub(r'\s+', '', parent or '').casefold()
    tokens = []
    for child in children:
        child_key = re.sub(r'\s+', '', child or '')
        folded = child_key.casefold()
        suffix = child_key[len(parent_key):] if folded.startswith(parent_key) \
            else ''
        if not re.fullmatch(r'[A-Za-z]{1,3}', suffix):
            match = re.search(r'(?:\s|\d)([A-Z]{1,2})$', child or '')
            suffix = match.group(1) if match else ''
        tokens.append(suffix.upper())
    if any(not token for token in tokens) or len(set(tokens)) != len(tokens):
        return [None] * len(children)
    return tokens


def _fit_shorthand(components):
    grouped = {}
    order = []
    for component in components:
        key = component['owner_key']
        if key not in grouped:
            grouped[key] = {'stellar': 0, 'disk': 0, 'row': component}
            order.append(key)
        grouped[key][component['role']] += 1
    parts = []
    for key in order:
        values = grouped[key]
        owner = values['row']['component']
        if owner is None and len(order) > 1:
            owner = values['row']['owner']
        symbols = '★' * values['stellar'] + '○' * values['disk']
        parts.append({
            'owner': owner,
            'symbols': symbols,
            'colour': values['row']['colour'],
            'title': '{}: {} stellar, {} disk'.format(
                values['row']['owner'], values['stellar'], values['disk']
            ),
        })
    return parts


def _target_source(package_directory, target, observation):
    if observation is not None:
        return Path(observation.source)
    if target.sdbid is not None:
        candidate = package_directory / '{}-rawphot.txt'.format(target.sdbid)
        if candidate.is_file():
            return candidate
    return None


def _external_links(main_id, ra, dec, radius_arcsec=20.0):
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
    radius_arcsec = max(1.0, float(radius_arcsec))
    spitzer_radius = max(180.0, radius_arcsec)
    esasky_fov = max(0.25, 2.4 * radius_arcsec / 3600.0)
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
                'ra': ra, 'dec': dec, 'radius': round(radius_arcsec, 2),
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
                    'api': 'search', 'searchoption': 'POSITION',
                    'sr': '{:g}s'.format(spitzer_radius),
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
                'action': 'goto', 'fov': '{:g}'.format(esasky_fov),
                'cooframe': 'J2000',
                'sci': 'true', 'hips': 'AllWISE color', 'target': coordinate,
            }),
        },
    ])
    return links


def _angular_separation_arcsec(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = map(math.radians, (ra1, dec1, ra2, dec2))
    delta_ra = ra2 - ra1
    delta_dec = dec2 - dec1
    value = (
        math.sin(delta_dec / 2) ** 2 +
        math.cos(dec1) * math.cos(dec2) * math.sin(delta_ra / 2) ** 2
    )
    angle = 2 * math.asin(min(1.0, math.sqrt(value)))
    return math.degrees(angle) * 3600.0


def _search_radius_arcsec(ra, dec, positions):
    if ra is None or dec is None:
        return 20.0
    separations = [
        _angular_separation_arcsec(ra, dec, other_ra, other_dec)
        for other_ra, other_dec in positions
        if other_ra is not None and other_dec is not None
    ]
    return max(20.0, (max(separations) + 10.0) if separations else 20.0)


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
