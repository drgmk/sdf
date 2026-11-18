from .context import sdf
import sdf.candidate_search
import sdf.model_policy


def _policy(key, stellar_models=('phoenix_sol',), max_disks=1):
    return sdf.model_policy.TargetModelPolicy(
        target_key=key,
        stellar_models=stellar_models,
        disk_models=('modbb_disk_r',),
        max_disks=max_disks,
    )


def test_joint_roots_cover_stellar_alternatives_per_target():
    search = sdf.candidate_search.CandidateSearch((
        _policy('simbad:B', ('koester_wd',)),
        _policy('simbad:A', ('phoenix_sol', 'phoenix_cool')),
    ))

    roots = search.roots()

    assert len(roots) == 2
    assert [
        [(component.owner_key, component.model_name)
         for component in candidate.components]
        for candidate in roots
    ] == [
        [('simbad:A', 'phoenix_sol'), ('simbad:B', 'koester_wd')],
        [('simbad:A', 'phoenix_cool'), ('simbad:B', 'koester_wd')],
    ]


def test_children_include_each_possible_disk_owner_without_order_bias():
    search = sdf.candidate_search.CandidateSearch((
        _policy('simbad:A'), _policy('simbad:B')
    ))
    root = search.roots()[0]

    first_level = search.children(root)
    second_level_from_a = search.children(next(
        candidate for candidate in first_level
        if any(component.owner_key == 'simbad:A' and component.role == 'disk'
               for component in candidate.components)
    ))
    second_level_from_b = search.children(next(
        candidate for candidate in first_level
        if any(component.owner_key == 'simbad:B' and component.role == 'disk'
               for component in candidate.components)
    ))

    common = {
        candidate.canonical_id for candidate in second_level_from_a
    } & {
        candidate.canonical_id for candidate in second_level_from_b
    }
    assert len(first_level) == 2
    assert len(common) == 1


def test_exhaustive_search_deduplicates_joint_candidate_graph():
    search = sdf.candidate_search.CandidateSearch((
        _policy('simbad:A'), _policy('simbad:B')
    ))

    candidates = search.exhaustive()

    # no disks, disk A, disk B, disks A+B
    assert len(candidates) == 4
    assert len({candidate.canonical_id for candidate in candidates}) == 4


def test_next_layer_completes_all_equivalent_disk_allocations():
    search = sdf.candidate_search.CandidateSearch((
        _policy('simbad:A', max_disks=2),
        _policy('simbad:B', max_disks=2),
    ))

    one_disk = search.next_layer(search.roots())
    two_disks = search.next_layer(one_disk)

    assert len(one_disk) == 2
    assert {
        tuple(sorted(
            component.owner_key for component in candidate.components
            if component.role == 'disk'
        ))
        for candidate in two_disks
    } == {
        ('simbad:A', 'simbad:A'),
        ('simbad:A', 'simbad:B'),
        ('simbad:B', 'simbad:B'),
    }
