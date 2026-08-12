import pytest

from .context import sdf
import sdf.dataset
import sdf.model_policy


def _target(name, spectral_type):
    return sdf.dataset.TargetIdentity(
        key='simbad:' + name,
        main_id=name,
        sdbid='sdb-' + name,
        keywords={'sp_type': spectral_type},
    )


@pytest.mark.parametrize(
    'spectral_type,stellar_models,disk_models',
    [
        ('F3V', ('phoenix_sol',), ('modbb_disk_r',)),
        ('DA2', ('koester_wd',), ('bb_disk_r',)),
        ('L2', ('phoenix_cool',), ('modbb_disk_r',)),
        ('M7V', ('phoenix_cool', 'phoenix_sol'), ('modbb_disk_r',)),
        ('dM8', ('phoenix_cool', 'phoenix_sol'), ('modbb_disk_r',)),
    ],
)
def test_automatic_policy_matches_legacy_spectral_inference(
        spectral_type, stellar_models, disk_models):
    policy = sdf.model_policy.automatic_policy(
        _target('Test', spectral_type)
    )

    assert policy.stellar_models == stellar_models
    assert policy.disk_models == disk_models
    assert policy.max_disks == 2
    assert policy.origin == 'automatic'


def test_package_policy_uses_simbad_name_and_overrides_only_named_target(
        tmp_path):
    targets = (_target('HD  1 A', 'F3V'), _target('HD 1 B', 'DA2'))
    dataset = sdf.dataset.FitDataset(targets=targets, observations=())
    (tmp_path / 'sdf-fit.yml').write_text(
        'targets:\n'
        '  HD 1 B:\n'
        '    stellar_models: [phoenix_cool]\n'
        '    disk_models: [modbb_disk_r]\n'
        '    max_disks: 0\n'
        '    search_mode: exact\n',
        encoding='utf-8',
    )

    policies = sdf.model_policy.resolve_package_policies(dataset, tmp_path)
    by_key = {policy.target_key: policy for policy in policies}

    assert by_key['simbad:HD  1 A'].origin == 'automatic'
    assert by_key['simbad:HD 1 B'].stellar_models == ('phoenix_cool',)
    assert by_key['simbad:HD 1 B'].max_disks == 0
    assert by_key['simbad:HD 1 B'].search_mode == 'exact'
    assert by_key['simbad:HD 1 B'].origin == 'package'


def test_package_policy_rejects_wrong_model_role(tmp_path):
    dataset = sdf.dataset.FitDataset(
        targets=(_target('HD 1', 'F3V'),), observations=()
    )
    (tmp_path / 'sdf-fit.yml').write_text(
        'targets:\n  HD 1:\n    stellar_models: [modbb_disk_r]\n',
        encoding='utf-8',
    )

    with pytest.raises(sdf.utils.SdfError, match='disk model, not stellar'):
        sdf.model_policy.resolve_package_policies(dataset, tmp_path)
