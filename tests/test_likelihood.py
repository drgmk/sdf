import numpy as np
import pytest

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


def test_joint_projection_adds_only_visible_target_models():
    cached_model = sdf.model.PhotModel.read_model('kurucz-0.0')
    filter_name = next(
        name for name, bases in zip(
            cached_model.filters, cached_model.colour_bases
        ) if len(bases) == 0
    )
    targets = tuple(
        sdf.dataset.TargetIdentity(
            key='simbad:' + name, main_id=name, sdbid=name, keywords={}
        )
        for name in ('A', 'B')
    )

    def observation(observation_id, contributors):
        photometry = sdf.photometry.Photometry(
            filters=np.array([filter_name]),
            fnujy=np.array([1.0]),
            e_fnujy=np.array([0.1]),
            upperlim=np.array([False]),
            ignore=np.array([False]),
            bibcode=np.array(['test']),
        )
        return sdf.dataset.Observation(
            observation_id=observation_id,
            source='/tmp/{}-rawphot.txt'.format(observation_id),
            target_key='simbad:' + observation_id,
            contributor_keys=tuple(
                'simbad:' + contributor for contributor in contributors
            ),
            photometry=photometry,
        )

    dataset = sdf.dataset.FitDataset(
        targets=targets,
        observations=(
            observation('A', ('A',)),
            observation('B', ('B',)),
            observation('AB', ('A', 'B')),
        ),
    )
    candidate = sdf.model_selection.ModelCandidate.from_owned_models((
        ('simbad:A', 'kurucz-0.0'),
        ('simbad:B', 'kurucz-0.0'),
    ))

    problem = sdf.likelihood.FitProblem.compile(dataset, candidate)
    parameters = problem.transform_prior(np.repeat(0.5, problem.ndim))
    _, _, _, _, _, obs_nel, _, _, _ = sdf.fitting.concat_obs(
        problem.observations
    )
    model_flux, component_flux = sdf.model.model_fluxes(
        problem.models, parameters, obs_nel,
        projection=problem.projection,
    )

    assert problem.projection.tolist() == [
        [True, False], [False, True], [True, True]
    ]
    assert [component.ordinal for component in candidate.components] == [1, 1]
    assert problem.parameter_names[0].startswith(
        'simbad:A:stellar:kurucz-0.0:1.'
    )
    assert problem.parameter_names[problem.ndim // 2].startswith(
        'simbad:B:stellar:kurucz-0.0:1.'
    )
    assert component_flux[0, 1] == 0.0
    assert component_flux[1, 0] == 0.0
    assert model_flux[2] == pytest.approx(model_flux[0] + model_flux[1])


def test_joint_projection_handles_colour_observations_without_nan_components():
    colour_name = next(
        name for name, bases in zip(
            sdf.model.PhotModel.read_model('kurucz-0.0').filters,
            sdf.model.PhotModel.read_model('kurucz-0.0').colour_bases,
        ) if len(bases) > 0
    )
    targets = tuple(
        sdf.dataset.TargetIdentity(
            key='simbad:' + name, main_id=name, sdbid=name, keywords={}
        )
        for name in ('A', 'B')
    )
    photometry = sdf.photometry.Photometry(
        filters=np.array([colour_name]), fnujy=np.array([0.1]),
        e_fnujy=np.array([0.01]), upperlim=np.array([False]),
        ignore=np.array([False]), bibcode=np.array(['test']),
    )
    dataset = sdf.dataset.FitDataset(
        targets=targets,
        observations=(sdf.dataset.Observation(
            observation_id='A', source='/tmp/A-rawphot.txt',
            target_key='simbad:A', contributor_keys=('simbad:A',),
            photometry=photometry,
        ),),
    )
    candidate = sdf.model_selection.ModelCandidate.from_owned_models((
        ('simbad:A', 'kurucz-0.0'),
        ('simbad:B', 'kurucz-0.0'),
    ))

    problem = sdf.likelihood.FitProblem.compile(dataset, candidate)
    parameters = problem.transform_prior(np.repeat(0.5, problem.ndim))
    _, _, _, _, _, obs_nel, _, _, _ = sdf.fitting.concat_obs(
        problem.observations
    )
    total, components = sdf.model.model_fluxes(
        problem.models, parameters, obs_nel,
        projection=problem.projection,
    )

    assert np.all(np.isfinite(total))
    assert np.all(np.isfinite(components))
    assert components[1, 0] == 0.0


def test_joint_upper_limit_applies_to_projected_sum():
    cached_model = sdf.model.PhotModel.read_model('kurucz-0.0')
    filter_name = next(
        name for name, bases in zip(
            cached_model.filters, cached_model.colour_bases
        ) if len(bases) == 0
    )
    targets = tuple(
        sdf.dataset.TargetIdentity(
            key='simbad:' + name, main_id=name, sdbid=name, keywords={}
        )
        for name in ('A', 'B')
    )
    photometry = sdf.photometry.Photometry(
        filters=np.array([filter_name]), fnujy=np.array([1e-300]),
        e_fnujy=np.array([1e-301]), upperlim=np.array([True]),
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
        ('simbad:A', 'kurucz-0.0'), ('simbad:B', 'kurucz-0.0'),
    ))
    problem = sdf.likelihood.FitProblem.compile(dataset, candidate)
    parameters = problem.transform_prior(np.repeat(0.5, problem.ndim))

    residual, _, _ = problem.residual(parameters)

    assert residual[0] < 0


def test_joint_spectrum_projection_and_parameter_labels():
    cached_model = sdf.model.PhotModel.read_model('kurucz-0.0')
    filter_name = next(
        name for name, bases in zip(
            cached_model.filters, cached_model.colour_bases
        ) if len(bases) == 0
    )
    target = sdf.dataset.TargetIdentity(
        key='simbad:A', main_id='A', sdbid='A', keywords={}
    )
    photometry = sdf.photometry.Photometry(
        filters=np.array([filter_name]), fnujy=np.array([1.0]),
        e_fnujy=np.array([0.1]), upperlim=np.array([False]),
        ignore=np.array([False]), bibcode=np.array(['test']),
    )
    spectrum = sdf.spectrum.ObsSpectrum(
        wavelength=np.array([1.0, 1.1]), fnujy=np.array([1.0, 1.0]),
        e_fnujy=np.array([0.1, 0.1]), bibcode='test', instrument='test',
    )
    dataset = sdf.dataset.FitDataset(
        targets=(target,), observations=(sdf.dataset.Observation(
            observation_id='A', source='/tmp/A-rawphot.txt',
            target_key=target.key, contributor_keys=(target.key,),
            photometry=photometry, spectra=(spectrum,),
        ),),
    )
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        target.key, ('kurucz-0.0',)
    )

    problem = sdf.likelihood.FitProblem.compile(dataset, candidate)
    parameters = problem.transform_prior(np.repeat(0.5, problem.ndim))
    residual, _, _ = problem.residual(parameters)

    assert problem.projection.tolist() == [[True], [True]]
    assert problem.parameter_names[-1] == 'observation:A.spectrum:1.norm'
    assert len(residual) == 3
    assert np.all(np.isfinite(residual))
