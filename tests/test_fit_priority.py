from types import SimpleNamespace

from .context import sdf
import sdf.fit_priority
import sdf.model_policy
import sdf.model_selection


def _stored(owned_models, evidence, error=None, n_parameters=4):
    candidate = sdf.model_selection.ModelCandidate.from_owned_models(
        owned_models
    )
    return SimpleNamespace(
        candidate_id=candidate.canonical_id,
        metadata={'components': [
            {
                'component_id': component.component_id,
                'owner_key': component.owner_key,
                'role': component.role,
                'model_name': component.model_name,
                'ordinal': component.ordinal,
            }
            for component in candidate.components
        ]},
        posterior=SimpleNamespace(
            evidence=evidence,
            evidence_error=error,
            parameter_names=tuple(
                'p{}'.format(index) for index in range(n_parameters)
            ),
        ),
    )


def _policy(owner):
    return sdf.model_policy.TargetModelPolicy(
        target_key=owner,
        stellar_models=('kurucz-0.0',),
        disk_models=('modbb_disk_r',),
        max_disks=2,
    )


def test_recommendation_chooses_simplest_candidate_within_evidence_threshold():
    root = _stored((('A', 'kurucz-0.0'),), 0.0, n_parameters=2)
    one = _stored((
        ('A', 'kurucz-0.0'), ('A', 'modbb_disk_r'),
    ), 7.0, n_parameters=5)
    two = _stored((
        ('A', 'kurucz-0.0'),
        ('A', 'modbb_disk_r'), ('A', 'modbb_disk_r'),
    ), 10.0, n_parameters=8)

    recommendation = sdf.fit_priority.recommend(
        (root, one, two), policies=(_policy('A'),),
        owner_order=('A',), threshold=5.0,
    )

    assert recommendation.candidate_id == one.candidate_id
    assert set(recommendation.competitive_ids) == {
        one.candidate_id, two.candidate_id,
    }


def test_equivalent_disk_placements_prefer_concentration_then_primary():
    stellar = (('A', 'kurucz-0.0'), ('B', 'kurucz-0.0'))
    around_a = _stored(stellar + (
        ('A', 'modbb_disk_r'), ('A', 'modbb_disk_r'),
    ), 10.0)
    divided = _stored(stellar + (
        ('A', 'modbb_disk_r'), ('B', 'modbb_disk_r'),
    ), 14.0)
    around_b = _stored(stellar + (
        ('B', 'modbb_disk_r'), ('B', 'modbb_disk_r'),
    ), 13.0)

    recommendation = sdf.fit_priority.recommend(
        (divided, around_b, around_a),
        policies=(_policy('A'), _policy('B')),
        owner_order=('A', 'B'), threshold=5.0,
    )

    assert recommendation.candidate_id == around_a.candidate_id


def test_evidence_uncertainty_conservatively_keeps_simpler_fit_competitive():
    one = _stored((
        ('A', 'kurucz-0.0'), ('A', 'modbb_disk_r'),
    ), 4.8, error=0.2, n_parameters=5)
    two = _stored((
        ('A', 'kurucz-0.0'),
        ('A', 'modbb_disk_r'), ('A', 'modbb_disk_r'),
    ), 10.0, error=0.2, n_parameters=8)

    recommendation = sdf.fit_priority.recommend(
        (one, two), owner_order=('A',), threshold=5.0,
    )

    assert recommendation.candidate_id == one.candidate_id
