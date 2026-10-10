"""Operator operations over retained fit candidates."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .candidate_library import CandidateLibrary, FitSelection
from .coordinator import FitCoordinator
from .dataset import FitDataset
from .joint_plotting import write_joint_results_page
from .model_policy import (
    ModelFamilyRegistry, TargetModelPolicy, resolve_package_policies,
)
from .model_selection import ModelCandidate, ModelComponentSpec
from .result_index import ResultIndex


@dataclass(frozen=True)
class AvailableFit:
    reference: str
    candidate_id: str
    display_name: str
    evidence: float | None
    completed_at: str
    is_default: bool
    display_position: int | None


class FitManager:
    """List and select retained candidates without running a sampler."""

    def __init__(self, package_directory):
        self.package_directory = Path(
            package_directory
        ).expanduser().resolve()
        self.library = CandidateLibrary(self.package_directory)

    def available(self):
        candidates = self.library.candidates()
        if not candidates:
            return ()
        try:
            selection = self.library.read_selection()
        except FileNotFoundError:
            selection = None
        positions = {
            candidate_id: position
            for position, candidate_id in enumerate(
                selection.displayed if selection is not None else ()
            )
        }
        default = selection.default if selection is not None else None
        return tuple(
            AvailableFit(
                reference=candidate.path.name,
                candidate_id=candidate.candidate_id,
                display_name=candidate.metadata['display_name'],
                evidence=candidate.posterior.evidence,
                completed_at=candidate.metadata['completed_at'],
                is_default=candidate.candidate_id == default,
                display_position=positions.get(candidate.candidate_id),
            )
            for candidate in sorted(
                candidates,
                key=lambda value: (
                    value.posterior.evidence is not None,
                    value.posterior.evidence
                    if value.posterior.evidence is not None else float('-inf'),
                ),
                reverse=True,
            )
        )

    def select(self, *, default=None, displayed=None):
        """Update selection using exact or unambiguous short references."""

        available = self.available()
        if not available:
            raise ValueError('package has no retained candidates')
        try:
            current = self.library.read_selection()
        except FileNotFoundError:
            current = FitSelection(
                default=available[0].candidate_id,
                displayed=tuple(fit.candidate_id for fit in available),
            )
        default_id = (
            current.default if default is None else self.resolve(default)
        )
        displayed_ids = (
            tuple(current.displayed) if displayed is None else
            tuple(self.resolve(reference) for reference in displayed)
        )
        displayed_ids = tuple(
            candidate_id for candidate_id in displayed_ids
            if candidate_id != default_id
        )
        selection = FitSelection(
            default=default_id,
            displayed=(default_id,) + displayed_ids,
            labels=current.labels,
        )
        self.library.write_selection(selection)
        index = ResultIndex.find_for_package(self.package_directory)
        if index is not None:
            index.update_package(self.package_directory)
        return selection

    def resolve(self, reference):
        """Resolve a full ID, storage name, or unambiguous prefix."""

        reference = str(reference)
        candidates = self.library.candidates()
        exact = [
            candidate for candidate in candidates
            if reference in {
                candidate.candidate_id, candidate.path.name,
                candidate.metadata['display_name'],
            }
        ]
        if len(exact) == 1:
            return exact[0].candidate_id
        prefixes = [
            candidate for candidate in candidates
            if candidate.candidate_id.startswith(reference) or
            candidate.path.name.startswith(reference)
        ]
        if len(prefixes) == 1:
            return prefixes[0].candidate_id
        if not prefixes:
            raise ValueError('unknown fit reference: {}'.format(reference))
        raise ValueError('ambiguous fit reference: {}'.format(reference))

    def render(self, output=None, fits=None, dataset=None, operator=False,
               policy_store=None, operator_api=None, operator_token=None,
               operator_actor=None, operator_actor_locked=False):
        candidate_ids = None
        if fits is not None:
            candidate_ids = tuple(self.resolve(reference) for reference in fits)
        return write_joint_results_page(
            self.package_directory, output=output,
            candidate_ids=candidate_ids, dataset=dataset, operator=operator,
            policy_store=policy_store, operator_api=operator_api,
            operator_token=operator_token, operator_actor=operator_actor,
            operator_actor_locked=operator_actor_locked,
        )

    def fit_exact(self, stellar, disks=(), *, force=False, nospec=False,
                  sampler=None, dataset=None):
        """Fit one explicitly owned model combination and retain it."""

        dataset = dataset or FitDataset.from_path(
            self.package_directory, nospec=nospec
        )
        candidate, policies = exact_candidate(
            dataset, stellar, disks,
            package_directory=self.package_directory,
        )
        coordinator = FitCoordinator(
            self.package_directory, sampler=sampler
        )
        return coordinator.fit_candidate(
            dataset, candidate, policy=policies, force=force
        )

    def refit(self, reference, *, nospec=False, sampler=None, dataset=None):
        """Force a fresh sampler run for one retained candidate."""

        candidate_id = self.resolve(reference)
        stored = self.library.find(candidate_id)
        candidate = ModelCandidate(components=tuple(
            ModelComponentSpec(**component)
            for component in stored.metadata['components']
        ))
        dataset = dataset or FitDataset.from_path(
            self.package_directory, nospec=nospec
        )
        policy = stored.metadata.get('policy', {}).get('targets', ())
        coordinator = FitCoordinator(
            self.package_directory, sampler=sampler
        )
        return coordinator.fit_candidate(
            dataset, candidate, policy=policy, force=True
        )


def exact_candidate(dataset, stellar, disks=(), registry=None,
                    package_directory=None):
    """Build and validate one exact hierarchy-owned candidate."""

    registry = registry or ModelFamilyRegistry.from_config()
    targets = {target.key: target for target in dataset.targets}
    aliases = {}
    for target in dataset.targets:
        for alias in (target.key, target.main_id, target.sdbid):
            if alias is not None:
                aliases.setdefault(_normalized(alias), set()).add(target.key)

    def owner_key(reference):
        matches = aliases.get(_normalized(reference), set())
        if not matches:
            raise ValueError('unknown target reference: {}'.format(reference))
        if len(matches) != 1:
            raise ValueError('ambiguous target reference: {}'.format(reference))
        return next(iter(matches))

    stellar_by_owner = {}
    for owner, model_name in stellar:
        key = owner_key(owner)
        registry.require(model_name, 'stellar')
        if key in stellar_by_owner:
            raise ValueError(
                'target {} has more than one stellar model'.format(owner)
            )
        stellar_by_owner[key] = model_name
    missing = set(targets) - set(stellar_by_owner)
    if missing:
        raise ValueError('stellar model missing for: {}'.format(
            ', '.join(targets[key].main_id or key for key in sorted(missing))
        ))

    disk_models = []
    for owner, model_name in disks:
        key = owner_key(owner)
        registry.require(model_name, 'disk')
        disk_models.append((key, model_name))
    disk_models.sort()
    owned_models = [
        (key, stellar_by_owner[key]) for key in sorted(targets)
    ] + disk_models
    candidate = ModelCandidate.from_owned_models(owned_models)

    defaults = {
        policy.target_key: policy
        for policy in resolve_package_policies(
            dataset, package_directory or Path('__no_package_policy__')
        )
    }
    policies = []
    for key in sorted(targets):
        owned_disks = [
            model_name for owner, model_name in disk_models if owner == key
        ]
        policies.append(TargetModelPolicy(
            target_key=key,
            stellar_models=(stellar_by_owner[key],),
            disk_models=tuple(sorted(set(owned_disks))),
            max_disks=len(owned_disks),
            search_mode='exact',
            reddening=defaults[key].reddening,
            origin='explicit',
        ))
    return candidate, tuple(policies)


def _normalized(value):
    return ' '.join(str(value).split()).casefold()
