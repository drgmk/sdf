"""Deterministic model-candidate graph generation."""

from __future__ import annotations

from collections import Counter, deque
from itertools import product

from .model_selection import ModelCandidate


class CandidateSearch:
    """Generate roots and disk additions without fitting-order dependence."""

    def __init__(self, policies):
        self.policies = tuple(sorted(policies, key=lambda item: item.target_key))
        keys = [policy.target_key for policy in self.policies]
        if len(set(keys)) != len(keys):
            raise ValueError('target policies must have unique target keys')

    def roots(self):
        choices = [policy.stellar_models for policy in self.policies]
        return tuple(
            self._candidate(stellar_models, Counter())
            for stellar_models in product(*choices)
        )

    def children(self, candidate):
        stellar_models, disk_counts = self._decompose(candidate)
        children = {}
        for policy in self.policies:
            current = sum(
                count for (owner, _), count in disk_counts.items()
                if owner == policy.target_key
            )
            if current >= policy.max_disks:
                continue
            for disk_model in policy.disk_models:
                child_counts = disk_counts.copy()
                child_counts[(policy.target_key, disk_model)] += 1
                child = self._candidate(stellar_models, child_counts)
                children[child.canonical_id] = child
        return tuple(children[key] for key in sorted(children))

    def exhaustive(self):
        """Return every unique candidate permitted by the policies."""

        candidates = {}
        queue = deque(self.roots())
        while queue:
            candidate = queue.popleft()
            if candidate.canonical_id in candidates:
                continue
            candidates[candidate.canonical_id] = candidate
            queue.extend(self.children(candidate))
        return tuple(candidates[key] for key in sorted(candidates))

    def next_layer(self, candidates):
        """Return the complete deduplicated disk-count layer after candidates."""

        children = {}
        for candidate in candidates:
            for child in self.children(candidate):
                children[child.canonical_id] = child
        return tuple(children[key] for key in sorted(children))

    def _candidate(self, stellar_models, disk_counts):
        owned_models = [
            (policy.target_key, model_name)
            for policy, model_name in zip(self.policies, stellar_models)
        ]
        for (owner, model_name), count in sorted(disk_counts.items()):
            owned_models.extend([(owner, model_name)] * count)
        return ModelCandidate.from_owned_models(owned_models)

    def _decompose(self, candidate):
        stellar_by_owner = {
            component.owner_key: component.model_name
            for component in candidate.components
            if component.role == 'stellar'
        }
        stellar_models = tuple(
            stellar_by_owner[policy.target_key] for policy in self.policies
        )
        disk_counts = Counter(
            (component.owner_key, component.model_name)
            for component in candidate.components if component.role == 'disk'
        )
        return stellar_models, disk_counts
