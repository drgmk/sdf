"""Thin orchestration for compiling, sampling, and retaining fits."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import subprocess

import numpy as np

from . import config as cfg
from .candidate_library import CandidateLibrary, FitSelection
from .likelihood import FitProblem
from .samplers import MultiNestSampler, PosteriorRun


class FitCoordinator:
    """Run candidates while keeping scientific and storage layers separate."""

    def __init__(self, package_directory, sampler=None, provenance=None):
        self.package_directory = Path(package_directory).expanduser().resolve()
        self.library = CandidateLibrary(self.package_directory)
        self.sampler = sampler or MultiNestSampler()
        self.provenance = _default_provenance()
        self.provenance.update(provenance or {})

    def fit_candidate(self, dataset, candidate, policy=(), force=False):
        """Fit and atomically retain one candidate, reusing it when current."""

        fit_state = _fit_state(candidate, self.sampler)
        existing = self.library.find(candidate.canonical_id)
        if (existing is not None and
                existing.is_current(dataset, fit_state=fit_state) and
                not force):
            return existing

        problem = FitProblem.compile(dataset, candidate)
        work_name = candidate.storage_name.rsplit('--', 1)[-1]
        work_directory = self.library.fits_directory / '.work' / work_name
        work_directory.mkdir(parents=True, exist_ok=True)
        output_base = str(work_directory / 'run_')
        posterior = self.sampler.run(problem, output_base)
        posterior = _cap_posterior(posterior, cfg.fitting['n_samples_max'])
        analysis = _basic_analysis(problem, posterior)
        policy_payload = [
            asdict(item) if hasattr(item, '__dataclass_fields__') else item
            for item in policy
        ]
        return self.library.commit(
            dataset,
            candidate,
            posterior,
            analysis=analysis,
            policy={'targets': policy_payload},
            provenance=self.provenance,
            fit_state=fit_state,
        )

    def evidence_guided(self, dataset, search, threshold=None, force=False):
        """Fit roots and expand every evidence-supported child branch."""

        threshold = (
            cfg.fitting['ev_threshold'] if threshold is None else threshold
        )
        fitted = {}
        expanded = set()
        queue = list(search.roots())
        roots = {candidate.canonical_id for candidate in queue}

        while queue:
            parent = queue.pop(0)
            if parent.canonical_id in expanded:
                continue
            parent_result = fitted.get(parent.canonical_id)
            if parent_result is None:
                parent_result = self.fit_candidate(
                    dataset, parent, policy=search.policies, force=force
                )
                fitted[parent.canonical_id] = parent_result
            expanded.add(parent.canonical_id)

            for child in search.children(parent):
                child_result = fitted.get(child.canonical_id)
                if child_result is None:
                    child_result = self.fit_candidate(
                        dataset, child, policy=search.policies, force=force
                    )
                    fitted[child.canonical_id] = child_result
                if (
                    _evidence(child_result) >
                    _evidence(parent_result) + threshold and
                    child.canonical_id not in expanded
                ):
                    queue.append(child)

        results = tuple(
            fitted[key] for key in sorted(fitted)
        )
        ranked = tuple(
            result.candidate_id for result in sorted(
                results, key=_evidence, reverse=True
            )
        )
        if ranked:
            self.library.write_selection(FitSelection(
                default=ranked[0], displayed=ranked
            ))
        return results, tuple(sorted(roots))


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


def _basic_analysis(problem, posterior):
    parameters = np.median(posterior.samples, axis=0)
    residuals, wavelengths, filters = problem.residual(parameters)
    return {
        'best_parameters': dict(zip(problem.parameter_names, parameters)),
        'chi_squared': float(np.sum(np.square(residuals))),
        'n_residuals': int(len(residuals)),
        'residuals': residuals,
        'wavelength_micron': wavelengths,
        'filters': filters,
    }


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
