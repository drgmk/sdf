"""Evidence-aware recommendations for retained model candidates."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class FitRecommendation:
    """Automatic choice and the candidates statistically competitive with it."""

    candidate_id: str
    competitive_ids: tuple[str, ...]
    best_evidence: float | None
    threshold: float


def recommend(candidates, *, policies=(), owner_order=(), threshold=5.0):
    """Choose the simplest fit not significantly worse than the best fit.

    Evidence establishes a competitive set. Optional-component count and the
    deterministic topology preferences are considered only inside that set.
    """

    candidates = tuple(candidates)
    if not candidates:
        raise ValueError('at least one candidate is required')
    if threshold < 0:
        raise ValueError('evidence threshold cannot be negative')

    evidenced = tuple(
        candidate for candidate in candidates
        if _finite(candidate.posterior.evidence)
    )
    if evidenced:
        best = max(evidenced, key=lambda value: value.posterior.evidence)
        best_evidence = float(best.posterior.evidence)
        best_error = _evidence_error(best)
        competitive = tuple(
            candidate for candidate in evidenced
            if best_evidence - float(candidate.posterior.evidence) <=
            threshold + math.hypot(best_error, _evidence_error(candidate))
        )
    else:
        best_evidence = None
        competitive = candidates

    minimum_disks = min(_disk_count(candidate) for candidate in competitive)
    simplest = tuple(
        candidate for candidate in competitive
        if _disk_count(candidate) == minimum_disks
    )
    order = _complete_owner_order(owner_order, simplest)
    model_ranks = _model_ranks(policies)
    selected = min(
        simplest,
        key=lambda candidate: _preference_key(
            candidate, order, model_ranks
        ),
    )
    return FitRecommendation(
        candidate_id=selected.candidate_id,
        competitive_ids=tuple(sorted(
            candidate.candidate_id for candidate in competitive
        )),
        best_evidence=best_evidence,
        threshold=float(threshold),
    )


def dataset_owner_order(dataset):
    """Return physical contributors in their first topology appearance."""

    owners = []
    observations = tuple(dataset.observations)
    for observation in sorted(
            observations,
            key=lambda value: len(value.contributor_keys) <= 1):
        for owner in observation.contributor_keys:
            if owner not in owners:
                owners.append(owner)
    for target in dataset.targets:
        if target.key not in owners:
            owners.append(target.key)
    return tuple(owners)


def _preference_key(candidate, owner_order, model_ranks):
    components = tuple(candidate.metadata['components'])
    disks = tuple(
        component for component in components if component['role'] == 'disk'
    )
    counts = Counter(component['owner_key'] for component in disks)
    distribution = tuple(-counts[owner] for owner in owner_order)
    owner_positions = {owner: index for index, owner in enumerate(owner_order)}
    ordered_components = sorted(
        components,
        key=lambda component: (
            owner_positions.get(component['owner_key'], len(owner_positions)),
            component['owner_key'],
            0 if component['role'] == 'stellar' else 1,
            component.get('ordinal', 1),
            component['model_name'],
        ),
    )
    family_preference = tuple(
        model_ranks.get(
            (component['owner_key'], component['role'], component['model_name']),
            (10**6, component['model_name']),
        )
        for component in ordered_components
    )
    return (
        len(counts),
        distribution,
        family_preference,
        len(candidate.posterior.parameter_names),
        candidate.candidate_id,
    )


def _model_ranks(policies):
    ranks = {}
    for policy in policies:
        for role, names in (
                ('stellar', policy.stellar_models), ('disk', policy.disk_models)):
            for position, name in enumerate(names):
                ranks[(policy.target_key, role, name)] = (position, name)
    return ranks


def _complete_owner_order(owner_order, candidates):
    order = list(dict.fromkeys(owner_order))
    present = {
        component['owner_key']
        for candidate in candidates
        for component in candidate.metadata['components']
    }
    order.extend(sorted(present - set(order)))
    return tuple(order)


def _disk_count(candidate):
    return sum(
        component['role'] == 'disk'
        for component in candidate.metadata['components']
    )


def _evidence_error(candidate):
    value = candidate.posterior.evidence_error
    return float(value) if _finite(value) else 0.0


def _finite(value):
    return value is not None and math.isfinite(float(value))
