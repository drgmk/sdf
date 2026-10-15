"""Thin orchestration for compiling, sampling, and retaining fits."""

from __future__ import annotations

from dataclasses import asdict, replace
from pathlib import Path
import subprocess

import numpy as np

from . import config as cfg
from . import utils
from .analysis import ANALYSIS_SCHEMA, FitAnalysis
from .candidate_library import (
    CandidateLibrary,
    FitSelection,
    fit_fingerprint,
    input_fingerprint,
)
from .candidate_search import CandidateSearch
from .dataset import FitDataset, canonical_fit_source
from .diagnostic_plots import ensure_candidate_diagnostic_plots
from .fit_priority import dataset_owner_order, recommend
from .likelihood import FitProblem
from .model_policy import resolve_package_policies
from .result_index import ResultIndex
from .samplers import MultiNestSampler, PosteriorRun


class FitCoordinator:
    """Run candidates while keeping scientific and storage layers separate."""

    def __init__(self, package_directory, sampler=None, provenance=None,
                 result_index=None):
        self.package_directory = Path(package_directory).expanduser().resolve()
        self.library = CandidateLibrary(self.package_directory)
        self.sampler = sampler or MultiNestSampler()
        self.result_index = (
            None if result_index is False else
            ResultIndex.find_for_package(self.package_directory)
            if result_index is None else result_index
        )
        self.provenance = _default_provenance()
        self.provenance.update(provenance or {})

    def fit_candidate(self, dataset, candidate, policy=(), force=False):
        """Fit and atomically retain one candidate, reusing it when current."""

        fit_state = _fit_state(candidate, self.sampler)
        existing = self.library.find(candidate.canonical_id)
        existing_is_current = (
            existing is not None and
            existing.is_current(dataset, fit_state=fit_state)
        )
        analysis_is_current = (
            existing_is_current and
            existing.metadata.get('analysis', {}).get('schema') ==
            ANALYSIS_SCHEMA
        )
        if analysis_is_current and not force:
            ensure_candidate_diagnostic_plots(existing)
            return existing

        problem = FitProblem.compile(dataset, candidate)
        if existing_is_current and not force:
            posterior = existing.posterior
        else:
            work_name = '{}-{}-{}'.format(
                candidate.storage_name.rsplit('--', 1)[-1],
                input_fingerprint(dataset)[:8],
                fit_fingerprint(fit_state)[:8],
            )
            work_directory = self.library.fits_directory / '.work' / work_name
            work_directory.mkdir(parents=True, exist_ok=True)
            # MultiNest's Fortran filenames are limited to about 100
            # characters. Run from the candidate work directory so even
            # deeply nested SDB package paths use the short `./run_` base.
            with utils.pushd(work_directory):
                posterior = self.sampler.run(problem, './run_')
            if posterior.native_output is not None:
                posterior = replace(
                    posterior,
                    native_output=str(
                        (work_directory / posterior.native_output).resolve().
                        relative_to(self.package_directory)
                    ),
                )
        posterior = _cap_posterior(posterior, cfg.fitting['n_samples_max'])
        analysis = FitAnalysis.run(problem, posterior)
        policy_payload = [
            asdict(item) if hasattr(item, '__dataclass_fields__') else item
            for item in policy
        ]
        result = self.library.commit(
            dataset,
            candidate,
            posterior,
            analysis=analysis.summary,
            analysis_arrays=analysis.arrays,
            policy={'targets': policy_payload},
            provenance=self.provenance,
            fit_state=fit_state,
        )
        ensure_candidate_diagnostic_plots(result, force=True)
        self._update_result_index()
        return result

    def evidence_guided(self, dataset, search, threshold=None, force=False):
        """Fit complete disk-count layers while evidence supports complexity."""

        threshold = (
            cfg.fitting['ev_threshold'] if threshold is None else threshold
        )
        fitted = {}
        layer = search.roots()
        roots = {candidate.canonical_id for candidate in layer}
        simpler_best = float('-inf')
        disk_count = 0

        while layer:
            layer_results = []
            for candidate in layer:
                result = self.fit_candidate(
                    dataset, candidate, policy=search.policies, force=force
                )
                fitted[candidate.canonical_id] = result
                layer_results.append(result)

            layer_best = max((_evidence(result) for result in layer_results))
            if disk_count and (
                    not np.isfinite(layer_best) or
                    layer_best <= simpler_best + threshold):
                break
            simpler_best = max(simpler_best, layer_best)
            layer = search.next_layer(layer)
            disk_count += 1

        results = tuple(
            fitted[key] for key in sorted(fitted)
        )
        self._write_evidence_selection(
            results, dataset, search.policies, threshold=threshold
        )
        return results, tuple(sorted(roots))

    def fit_all(self, dataset, candidates, policy=(), force=False):
        """Fit an explicit/exhaustive candidate collection exactly once."""

        results = tuple(
            self.fit_candidate(
                dataset, candidate, policy=policy, force=force
            )
            for candidate in candidates
        )
        self._write_evidence_selection(results, dataset, policy)
        return results

    def _write_evidence_selection(self, results, dataset, policy,
                                  threshold=None):
        threshold = (
            cfg.fitting['ev_threshold'] if threshold is None else threshold
        )
        ranked = tuple(
            result.candidate_id for result in sorted(
                results, key=_evidence, reverse=True
            )
        )
        if ranked:
            recommendation = recommend(
                results, policies=policy,
                owner_order=dataset_owner_order(dataset),
                threshold=threshold,
            )
            try:
                current = self.library.read_selection()
            except FileNotFoundError:
                current = None
            current_default = (
                None if current is None else
                self.library.find(current.default)
            )
            preserve_override = (
                current is not None and
                current.default_source == 'operator' and
                current_default is not None and
                current_default.is_current(dataset)
            )
            default = (
                current.default if preserve_override else
                recommendation.candidate_id
            )
            displayed = (default,) + tuple(
                candidate_id for candidate_id in ranked
                if candidate_id != default
            )
            self.library.write_selection(FitSelection(
                default=default,
                displayed=displayed,
                labels=(current.labels if current is not None else None),
                recommended=recommendation.candidate_id,
                default_source='operator' if preserve_override else 'automatic',
            ))
            self._update_result_index()
            return self.library.read_selection()
        return None

    def _update_result_index(self):
        if self.result_index is not None:
            self.result_index.update_package(self.package_directory)


def fit_package(source, *, sampler=None, force=False, nospec=False,
                provenance=None, policy_store=None):
    """Fit a discovered package using its resolved portable policy."""

    source = canonical_fit_source(source)
    dataset = FitDataset.from_path(source, nospec=nospec)
    package_directory = source if source.is_dir() else source.parent
    policies = resolve_package_policies(
        dataset, package_directory, policy_store=policy_store,
    )
    modes = {policy.search_mode for policy in policies}
    if len(modes) != 1:
        raise ValueError(
            'all target policies in one package must use the same search mode'
        )
    mode = modes.pop()
    search = CandidateSearch(policies)
    coordinator = FitCoordinator(
        package_directory, sampler=sampler, provenance=provenance
    )
    if mode == 'evidence-guided':
        results, _ = coordinator.evidence_guided(
            dataset, search, force=force
        )
        return results
    if mode == 'exhaustive':
        return coordinator.fit_all(
            dataset, search.exhaustive(), policy=policies, force=force
        )
    raise ValueError(
        '{} package fitting needs explicit candidate input, which is not '
        'implemented yet'.format(mode)
    )


def reprioritize_package(source, *, threshold=None, nospec=False,
                         policy_store=None, result_index=None):
    """Refresh selection metadata using current retained candidates only.

    This deliberately never compiles a fit or invokes a sampler. Candidates
    made from older package inputs are ignored, while a current explicit
    operator choice remains the default.
    """

    source = canonical_fit_source(source)
    dataset = FitDataset.from_path(source, nospec=nospec)
    package_directory = source if source.is_dir() else source.parent
    policies = resolve_package_policies(
        dataset, package_directory, policy_store=policy_store,
    )
    coordinator = FitCoordinator(
        package_directory, result_index=result_index,
    )
    current = tuple(
        candidate for candidate in coordinator.library.candidates()
        if candidate.is_current(dataset)
    )
    if not current:
        return None
    return coordinator._write_evidence_selection(
        current, dataset, policies, threshold=threshold,
    )


def _cap_posterior(posterior, maximum):
    if len(posterior.samples) <= maximum:
        return posterior
    selection = np.linspace(
        0, len(posterior.samples) - 1, maximum, dtype=int
    )
    return PosteriorRun(
        samples=posterior.samples[selection],
        parameter_names=posterior.parameter_names,
        sampler=posterior.sampler,
        log_likelihoods=(
            posterior.log_likelihoods[selection]
            if posterior.log_likelihoods is not None else None
        ),
        weights=(
            posterior.weights[selection]
            if posterior.weights is not None else None
        ),
        evidence=posterior.evidence,
        evidence_error=posterior.evidence_error,
        options=posterior.options,
        diagnostics=posterior.diagnostics,
        native_output=posterior.native_output,
    )


def _evidence(stored_candidate):
    evidence = stored_candidate.posterior.evidence
    if evidence is None:
        raise ValueError('evidence-guided search requires sampler evidence')
    return evidence


def _fit_state(candidate, sampler):
    model_files = {}
    for model_name in sorted(set(candidate.model_names)):
        root = Path(cfg.models['loc'][model_name])
        model_files[model_name] = [
            {
                'path': str(path.relative_to(root)),
                'size': path.stat().st_size,
                'mtime_ns': path.stat().st_mtime_ns,
            }
            for path in sorted(root.rglob('*')) if path.is_file()
        ]
    return {
        'candidate_id': candidate.canonical_id,
        'model_files': model_files,
        'model_normalization_range': cfg.fitting['model_om_range'],
        'spectrum_normalization_range': cfg.fitting['spectra_norm_range'],
        'spectral_model_preparation': 'gaussian-lsf-v1',
        'sampler': type(sampler).__name__,
        'sampler_options': {
            'n_live_points': cfg.fitting['n_live'],
            'n_iter_before_update': cfg.fitting['n_update'],
        },
    }


def _default_provenance():
    repository = Path(__file__).resolve().parents[1]
    provenance = {}
    try:
        completed = subprocess.run(
            ['git', 'rev-parse', 'HEAD'], cwd=repository,
            check=True, capture_output=True, text=True
        )
        provenance['sdf_revision'] = completed.stdout.strip()
        dirty = subprocess.run(
            ['git', 'status', '--porcelain', '--untracked-files=no'],
            cwd=repository, check=True, capture_output=True, text=True
        )
        provenance['sdf_dirty'] = bool(dirty.stdout.strip())
    except (OSError, subprocess.CalledProcessError):
        pass
    return provenance
