"""Operator operations over retained fit candidates."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .candidate_library import CandidateLibrary, FitSelection
from .joint_plotting import write_joint_results_page
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

    def render(self, output=None):
        return write_joint_results_page(
            self.package_directory, output=output
        )
