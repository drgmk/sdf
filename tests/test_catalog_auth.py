import os

import pytest

import sdf.catalog_auth


def test_htgroup_parser_supports_comments_blank_lines_and_duplicates():
    groups = sdf.catalog_auth.parse_htgroup('''
        # Current SDB access groups
        admin: grant alice
        reasons: alice bob
        reasons: grant
    ''')

    assert groups == {
        'admin': frozenset(('grant', 'alice')),
        'reasons': frozenset(('alice', 'bob', 'grant')),
    }


def test_htgroup_parser_rejects_malformed_lines():
    with pytest.raises(ValueError, match='line 2'):
        sdf.catalog_auth.parse_htgroup('admin: grant\nnot-a-group-line\n')


def test_proxy_authorizer_reads_and_reloads_group_membership(tmp_path):
    path = tmp_path / '.htgroup'
    path.write_text('admin: grant\nreasons: alice\n', encoding='utf-8')
    auth = sdf.catalog_auth.CatalogAuthorizer(
        user_header='X-SDF-User', group_file=path,
    )

    grant = auth.principal({'X-SDF-User': 'grant'})
    alice = auth.principal({'X-SDF-User': 'alice'})
    anonymous = auth.principal({})

    assert grant.groups == frozenset(('admin',))
    assert auth.can_access_sample(grant, 'reasons')
    assert auth.can_access_sample(alice, 'reasons')
    assert not anonymous.authenticated

    path.write_text('admin: grant\nreasons: alice bob\n', encoding='utf-8')
    os.utime(path, ns=(path.stat().st_atime_ns, path.stat().st_mtime_ns + 1))
    assert auth.principal({'X-SDF-User': 'bob'}).groups == frozenset(
        ('reasons',)
    )


def test_local_authorizer_preserves_loopback_operator_workflow():
    auth = sdf.catalog_auth.CatalogAuthorizer(local_actor='Grant')
    principal = auth.principal({})

    assert principal.username == 'Grant'
    assert auth.require_admin(principal) is principal
    assert auth.can_access_sample(principal, 'any-sample')


def test_proxy_authorizer_requires_both_header_and_group_file(tmp_path):
    with pytest.raises(ValueError, match='group file'):
        sdf.catalog_auth.CatalogAuthorizer(user_header='X-SDF-User')
    with pytest.raises(ValueError, match='proxy user header'):
        sdf.catalog_auth.CatalogAuthorizer(group_file=tmp_path / '.htgroup')
