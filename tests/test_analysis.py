import numpy as np
import pytest

from .context import sdf
import sdf.analysis
import sdf.dataset
import sdf.likelihood
import sdf.model_selection
import sdf.samplers


def test_joint_analysis_uses_each_stellar_owners_parallax():
    cached_model = sdf.model.PhotModel.read_model('kurucz-0.0')
    filter_name = next(
        name for name, bases in zip(
            cached_model.filters, cached_model.colour_bases
        ) if len(bases) == 0
    )
    targets = tuple(
        sdf.dataset.TargetIdentity(
            key='simbad:' + name,
            main_id=name,
            sdbid=name,
            keywords={'plx_value': parallax, 'plx_err': 0.0},
        )
        for name, parallax in (('A', 10.0), ('B', 20.0))
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
        observation_id='AB', source='/tmp/AB-rawphot.txt',
        target_key='simbad:AB',
        contributor_keys=('simbad:A', 'simbad:B'),
        photometry=photometry,
    )
    dataset = sdf.dataset.FitDataset(
        targets=targets, observations=(observation,)
    )
    candidate = sdf.model_selection.ModelCandidate.from_owned_models((
        ('simbad:A', 'kurucz-0.0'),
        ('simbad:B', 'kurucz-0.0'),
    ))
    problem = sdf.likelihood.FitProblem.compile(dataset, candidate)
    midpoint = problem.transform_prior(np.repeat(0.5, problem.ndim))
    posterior = sdf.samplers.PosteriorRun(
        samples=np.vstack((midpoint, midpoint)),
        parameter_names=problem.parameter_names,
        sampler='test',
    )

    analysis = sdf.analysis.FitAnalysis.run(problem, posterior)

    target_results = analysis.summary['targets']
    radius_a = target_results['simbad:A']['stellar_components'][0][
        'radius_rsun'
    ]['median']
    radius_b = target_results['simbad:B']['stellar_components'][0][
        'radius_rsun'
    ]['median']
    assert radius_a == pytest.approx(2.0 * radius_b)
    assert target_results['simbad:A']['stellar_components'][0][
        'model_name'
    ] == 'kurucz-0.0'
    assert analysis.summary['data_scopes'][0] == {
        'observation_id': 'AB',
        'contributors': ['simbad:A', 'simbad:B'],
        'kind': 'Photometry',
        'element_start': 0,
        'element_stop': 1,
    }
    assert analysis.arrays['component_prediction_median'].shape == (2, 1)


def test_disk_analysis_uses_only_its_owner_stellar_luminosity():
    cached_model = sdf.model.PhotModel.read_model('kurucz-0.0')
    filter_name = next(
        name for name, bases in zip(
            cached_model.filters, cached_model.colour_bases
        ) if len(bases) == 0
    )
    targets = tuple(
        sdf.dataset.TargetIdentity(
            key='simbad:' + name, main_id=name, sdbid=name,
            keywords={'plx_value': 10.0, 'plx_err': 0.0},
        )
        for name in ('A', 'B')
    )
    photometry = sdf.photometry.Photometry(
        filters=np.array([filter_name]), fnujy=np.array([1.0]),
        e_fnujy=np.array([0.1]), upperlim=np.array([False]),
        ignore=np.array([False]), bibcode=np.array(['test']),
    )
    dataset = sdf.dataset.FitDataset(
        targets=targets,
        observations=(sdf.dataset.Observation(
            observation_id='AB', source='/tmp/AB-rawphot.txt',
            target_key='simbad:AB',
            contributor_keys=('simbad:A', 'simbad:B'),
            photometry=photometry,
        ),),
    )
    candidate = sdf.model_selection.ModelCandidate.from_owned_models((
        ('simbad:A', 'kurucz-0.0'),
        ('simbad:B', 'kurucz-0.0'),
        ('simbad:A', 'bb_disk_r'),
    ))
    problem = sdf.likelihood.FitProblem.compile(dataset, candidate)
    midpoint = problem.transform_prior(np.repeat(0.5, problem.ndim))
    posterior = sdf.samplers.PosteriorRun(
        samples=np.vstack((midpoint, midpoint)),
        parameter_names=problem.parameter_names, sampler='test',
    )

    analysis = sdf.analysis.FitAnalysis.run(problem, posterior)

    disk_a = analysis.summary['targets']['simbad:A']['disk_components']
    disk_b = analysis.summary['targets']['simbad:B']['disk_components']
    assert len(disk_a) == 1
    assert disk_b == []
    assert disk_a[0]['model_name'] == 'bb_disk_r'
    assert disk_a[0]['fractional_luminosity']['median'] > 0
    assert disk_a[0]['blackbody_radius_au']['median'] > 0


def test_analysis_rejects_posterior_for_another_problem():
    target = sdf.dataset.TargetIdentity(
        key='simbad:A', main_id='A', sdbid='A', keywords={}
    )
    cached_model = sdf.model.PhotModel.read_model('kurucz-0.0')
    filter_name = next(
        name for name, bases in zip(
            cached_model.filters, cached_model.colour_bases
        ) if len(bases) == 0
    )
    photometry = sdf.photometry.Photometry(
        filters=np.array([filter_name]), fnujy=np.array([1.0]),
        e_fnujy=np.array([0.1]), upperlim=np.array([False]),
        ignore=np.array([False]), bibcode=np.array(['test']),
    )
    dataset = sdf.dataset.FitDataset(
        targets=(target,), observations=(sdf.dataset.Observation(
            observation_id='A', source='/tmp/A-rawphot.txt',
            target_key=target.key, contributor_keys=(target.key,),
            photometry=photometry,
        ),),
    )
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        target.key, ('kurucz-0.0',)
    )
    problem = sdf.likelihood.FitProblem.compile(dataset, candidate)
    samples = np.ones((2, problem.ndim))
    posterior = sdf.samplers.PosteriorRun(
        samples=samples,
        parameter_names=tuple('wrong{}'.format(i) for i in range(problem.ndim)),
        sampler='test',
    )

    with pytest.raises(ValueError, match='parameter names'):
        sdf.analysis.FitAnalysis.run(problem, posterior)
