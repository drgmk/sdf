from dataclasses import replace
import os

import numpy as np

from .context import sdf
from .test_likelihood import _dataset
from .test_dataset import _write_rawphot
import sdf.analysis
import sdf.candidate_search
import sdf.coordinator
import sdf.model_policy
import sdf.model_selection
import sdf.samplers


class FakeSampler:
    def __init__(self, evidence_by_disk_count=None):
        self.evidence_by_disk_count = evidence_by_disk_count or {}
        self.calls = []

    def run(self, problem, output_base):
        disk_count = sum(
            component.role == 'disk'
            for component in problem.candidate.components
        )
        self.calls.append((
            problem.candidate.canonical_id, output_base, os.getcwd()
        ))
        midpoint = problem.transform_prior(np.repeat(0.5, problem.ndim))
        return sdf.samplers.PosteriorRun(
            samples=np.vstack((midpoint, midpoint)),
            parameter_names=problem.parameter_names,
            sampler='fake',
            log_likelihoods=np.array([-1.0, -1.0]),
            evidence=self.evidence_by_disk_count.get(disk_count, 0.0),
        )


def _fit_dataset(tmp_path):
    cached_model = sdf.model.PhotModel.read_model('kurucz-0.0')
    filter_name = next(
        name for name, bases in zip(
            cached_model.filters, cached_model.colour_bases
        ) if len(bases) == 0
    )
    source = tmp_path / 'input-rawphot.txt'
    source.write_text('test input\n', encoding='utf-8')
    dataset = _dataset(filter_name)
    observation = replace(
        dataset.observations[0], source=str(source.resolve())
    )
    return replace(dataset, observations=(observation,))


def test_coordinator_commits_then_reuses_current_candidate(tmp_path):
    dataset = _fit_dataset(tmp_path)
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('kurucz-0.0',)
    )
    sampler = FakeSampler()
    coordinator = sdf.coordinator.FitCoordinator(
        tmp_path, sampler=sampler, provenance={'sdf_revision': 'test'}
    )

    first = coordinator.fit_candidate(dataset, candidate)
    second = coordinator.fit_candidate(dataset, candidate)

    assert len(sampler.calls) == 1
    assert first.candidate_id == second.candidate_id
    assert first.metadata['analysis']['n_residuals'] == 1
    assert first.metadata['provenance']['sdf_revision'] == 'test'
    assert sampler.calls[0][1] == './run_'
    assert '/fits/.work/' in sampler.calls[0][2]


def test_changed_inputs_use_a_fresh_resumable_work_directory(tmp_path):
    dataset = _fit_dataset(tmp_path)
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('kurucz-0.0',)
    )
    sampler = FakeSampler()
    coordinator = sdf.coordinator.FitCoordinator(tmp_path, sampler=sampler)
    coordinator.fit_candidate(dataset, candidate)
    with open(dataset.observations[0].source, 'a', encoding='utf-8') as handle:
        handle.write('changed input\n')

    coordinator.fit_candidate(dataset, candidate)

    assert len(sampler.calls) == 2
    assert sampler.calls[0][2] != sampler.calls[1][2]
    assert all('/fits/.work/' in call[2] for call in sampler.calls)


def test_coordinator_reanalyses_old_candidate_without_resampling(tmp_path):
    dataset = _fit_dataset(tmp_path)
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('kurucz-0.0',)
    )
    sampler = FakeSampler()
    coordinator = sdf.coordinator.FitCoordinator(tmp_path, sampler=sampler)
    first = coordinator.fit_candidate(dataset, candidate)
    result_path = first.path / 'result.json'
    payload = result_path.read_text(encoding='utf-8').replace(
        '"schema": {}'.format(sdf.analysis.ANALYSIS_SCHEMA),
        '"schema": 1',
    )
    result_path.write_text(payload, encoding='utf-8')

    refreshed = coordinator.fit_candidate(dataset, candidate)

    assert len(sampler.calls) == 1
    assert refreshed.metadata['analysis']['schema'] == \
        sdf.analysis.ANALYSIS_SCHEMA
    assert 'component_spectrum_median' in refreshed.analysis_arrays


def test_evidence_guided_search_expands_all_supported_siblings(tmp_path):
    dataset = _fit_dataset(tmp_path)
    policy = sdf.model_policy.TargetModelPolicy(
        target_key=dataset.targets[0].key,
        stellar_models=('kurucz-0.0',),
        disk_models=('modbb_disk_r',),
        max_disks=2,
    )
    search = sdf.candidate_search.CandidateSearch((policy,))
    sampler = FakeSampler({0: 0.0, 1: 6.0, 2: 7.0})
    coordinator = sdf.coordinator.FitCoordinator(tmp_path, sampler=sampler)

    fitted, roots = coordinator.evidence_guided(
        dataset, search, threshold=5.0
    )

    # The supported one-disk level causes the complete two-disk level to be
    # fitted. Its small improvement leaves the simpler one-disk fit preferred.
    assert len(fitted) == 3
    assert len(roots) == 1
    assert len(sampler.calls) == 3
    assert sorted(result.posterior.evidence for result in fitted) == [0.0, 6.0, 7.0]
    selection = coordinator.library.read_selection()
    assert selection.default == next(
        result.candidate_id for result in fitted
        if result.posterior.evidence == 6.0
    )
    assert selection.recommended == selection.default
    assert selection.default_source == 'automatic'
    assert len(selection.displayed) == 3


def test_evidence_guided_search_requires_more_than_configured_threshold(
        tmp_path):
    dataset = _fit_dataset(tmp_path)
    policy = sdf.model_policy.TargetModelPolicy(
        target_key=dataset.targets[0].key,
        stellar_models=('kurucz-0.0',),
        disk_models=('modbb_disk_r',),
        max_disks=2,
    )
    search = sdf.candidate_search.CandidateSearch((policy,))
    coordinator = sdf.coordinator.FitCoordinator(
        tmp_path, sampler=FakeSampler({0: 0.0, 1: 5.0, 2: 100.0})
    )

    fitted, _ = coordinator.evidence_guided(
        dataset, search, threshold=5.0
    )

    assert len(fitted) == 2
    assert all(
        sum(component['role'] == 'disk'
            for component in result.metadata['components']) < 2
        for result in fitted
    )
    selection = coordinator.library.read_selection()
    assert sum(
        component['role'] == 'disk'
        for component in coordinator.library.find(
            selection.default
        ).metadata['components']
    ) == 0


def test_evidence_selection_preserves_operator_default(tmp_path):
    dataset = _fit_dataset(tmp_path)
    policy = sdf.model_policy.TargetModelPolicy(
        target_key=dataset.targets[0].key,
        stellar_models=('kurucz-0.0',),
        disk_models=('modbb_disk_r',),
        max_disks=2,
    )
    search = sdf.candidate_search.CandidateSearch((policy,))
    coordinator = sdf.coordinator.FitCoordinator(
        tmp_path, sampler=FakeSampler({0: 0.0, 1: 6.0, 2: 7.0})
    )
    fitted, _ = coordinator.evidence_guided(dataset, search, threshold=5.0)
    explicit = coordinator.fit_candidate(
        dataset,
        sdf.model_selection.ModelCandidate.for_single_target(
            dataset.targets[0].key, ('phoenix_sol',)
        ),
    )
    coordinator.library.write_selection(sdf.candidate_library.FitSelection(
        default=explicit.candidate_id,
        displayed=(explicit.candidate_id,) + tuple(
            result.candidate_id for result in fitted
        ),
        recommended=coordinator.library.read_selection().recommended,
        default_source='operator',
    ))

    coordinator.evidence_guided(dataset, search, threshold=5.0)

    selection = coordinator.library.read_selection()
    assert selection.default == explicit.candidate_id
    assert selection.default != selection.recommended
    assert selection.default_source == 'operator'


def test_fit_package_loads_yaml_policy_and_runs_exhaustively(tmp_path):
    _write_rawphot(
        tmp_path / 'target-rawphot.txt',
        sdbid='target', main_id='HD 123'
    )
    (tmp_path / 'joint-fit.yml').write_text(
        'version: 1\nobservations: {}\n', encoding='utf-8'
    )
    (tmp_path / 'sdf-fit.yml').write_text(
        'targets:\n'
        '  HD 123:\n'
        '    stellar_models: [kurucz-0.0]\n'
        '    max_disks: 0\n'
        '    search_mode: exhaustive\n',
        encoding='utf-8',
    )
    sampler = FakeSampler()

    results = sdf.coordinator.fit_package(tmp_path, sampler=sampler)

    assert len(results) == 1
    assert len(sampler.calls) == 1
    assert results[0].metadata['policy']['targets'][0]['origin'] == 'package'
