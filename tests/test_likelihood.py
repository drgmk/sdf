import numpy as np

from .context import sdf
import sdf.dataset
import sdf.likelihood
import sdf.model_selection


def _dataset(filter_name='2MJ'):
    target = sdf.dataset.TargetIdentity(
        key='simbad:Test', main_id='Test', sdbid='sdb-test', keywords={}
    )
    photometry = sdf.photometry.Photometry(
        filters=np.array([filter_name]),
        fnujy=np.array([1.0]),
        e_fnujy=np.array([0.1]),
        upperlim=np.array([False]),
        ignore=np.array([False]),
        bibcode=np.array(['test']),
    )
    observation = sdf.dataset.Observation(
        observation_id='sdb-test',
        source='/tmp/test-rawphot.txt',
        target_key=target.key,
        contributor_keys=(target.key,),
        photometry=photometry,
    )
    return sdf.dataset.FitDataset(
        targets=(target,), observations=(observation,)
    )


def test_candidate_records_model_ownership_and_ordinals():
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        'simbad:Test',
        ('phoenix_sol', 'modbb_disk_r', 'modbb_disk_r'),
    )

    assert candidate.model_names == (
        'phoenix_sol', 'modbb_disk_r', 'modbb_disk_r'
    )
    assert [component.role for component in candidate.components] == [
        'stellar', 'disk', 'disk'
    ]
    assert [component.ordinal for component in candidate.components] == [
        1, 1, 2
    ]
    assert candidate.candidate_id == (
        'phoenix_sol+modbb_disk_r+modbb_disk_r'
    )


def test_fit_problem_wraps_existing_single_target_likelihood():
    # Other legacy tests mutate the cached model by selecting filters.  Use a
    # filter that is still present so this test characterizes FitProblem rather
    # than that independent cache behavior.
    filter_name = sdf.model.PhotModel.read_model('kurucz-0.0').filters[0]
    dataset = _dataset(filter_name)
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('kurucz-0.0',)
    )
    problem = sdf.likelihood.FitProblem.compile(dataset, candidate)
    parameters = problem.transform_prior(np.repeat(0.5, problem.ndim))

    assert problem.ndim == len(problem.parameter_ranges)
    assert problem.log_likelihood(parameters) == sdf.fitting.lnlike(
        parameters, problem.observations, problem.models
    )
