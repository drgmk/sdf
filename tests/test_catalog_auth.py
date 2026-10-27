import sqlite3

import pytest

import sdf.catalog_auth


def _access_database(tmp_path):
    path = tmp_path / 'sdb.sqlite'
    with sqlite3.connect(path) as connection:
        connection.executescript('''
            CREATE TABLE users (
                username TEXT PRIMARY KEY,
                is_admin INTEGER NOT NULL,
                disabled INTEGER NOT NULL
            );
            CREATE TABLE samples (
                id INTEGER PRIMARY KEY,
                name TEXT UNIQUE NOT NULL,
                owner_username TEXT,
                visibility TEXT NOT NULL
            );
            CREATE TABLE sample_grants (
                sample_id INTEGER NOT NULL,
                username TEXT NOT NULL,
                PRIMARY KEY (sample_id, username)
            );
            INSERT INTO users VALUES ('grant', 1, 0);
            INSERT INTO users VALUES ('alice', 0, 0);
            INSERT INTO users VALUES ('bob', 0, 0);
            INSERT INTO users VALUES ('disabled', 0, 1);
            INSERT INTO samples VALUES (1, 'public', NULL, 'public');
            INSERT INTO samples VALUES (2, 'private', 'alice', 'private');
            INSERT INTO samples VALUES (3, 'shared', NULL, 'private');
            INSERT INTO sample_grants VALUES (3, 'bob');
        ''')
    return path


def test_proxy_authorizer_reads_database_roles_and_sample_policy(tmp_path):
    path = _access_database(tmp_path)
    auth = sdf.catalog_auth.CatalogAuthorizer(
        user_header='X-SDF-User', access_database=path,
    )

    grant = auth.principal({'X-SDF-User': 'grant'})
    alice = auth.principal({'X-SDF-User': 'alice'})
    bob = auth.principal({'X-SDF-User': 'bob'})
    anonymous = auth.principal({})

    assert grant.is_admin
    assert auth.can_access_sample(grant, 'private')
    assert auth.can_access_sample(alice, 'private')
    assert auth.can_access_sample(bob, 'shared')
    assert not auth.can_access_sample(bob, 'private')
    assert not anonymous.authenticated

    with pytest.raises(sdf.catalog_auth.AuthorizationDenied, match='disabled'):
        auth.principal({'X-SDF-User': 'disabled'})
    with pytest.raises(sdf.catalog_auth.AuthorizationDenied, match='registered'):
        auth.principal({'X-SDF-User': 'unknown'})


def test_public_samples_and_signed_session_are_fail_closed(tmp_path):
    database = _access_database(tmp_path)
    auth = sdf.catalog_auth.CatalogAuthorizer(
        user_header='X-SDF-User', access_database=database,
        session_secret=b'test-secret',
    )

    anonymous = auth.principal({})
    assert auth.can_access_sample(anonymous, 'public')
    assert not auth.can_access_sample(anonymous, 'private')

    token = auth.issue_session('alice')
    alice = auth.principal({'Cookie': 'other=x; sdf_session=' + token})
    assert alice.username == 'alice'
    assert auth.can_access_sample(alice, 'public')
    assert auth.can_access_sample(alice, 'private')

    tampered = token[:-1] + ('a' if token[-1] != 'a' else 'b')
    assert not auth.principal({'Cookie': 'sdf_session=' + tampered}).authenticated
    expired = auth.issue_session('alice', now=0)
    assert not auth.principal({'Cookie': 'sdf_session=' + expired}).authenticated

    with sqlite3.connect(database) as connection:
        connection.execute(
            "INSERT INTO samples VALUES (4, 'new-public', NULL, 'public')"
        )
    assert auth.is_public_sample('new-public')


def test_local_authorizer_preserves_loopback_operator_workflow():
    auth = sdf.catalog_auth.CatalogAuthorizer(local_actor='Grant')
    principal = auth.principal({})

    assert principal.username == 'Grant'
    assert auth.require_admin(principal) is principal
    assert auth.can_access_sample(principal, 'any-sample')


def test_proxy_authorizer_requires_both_header_and_database(tmp_path):
    with pytest.raises(ValueError, match='authorization database'):
        sdf.catalog_auth.CatalogAuthorizer(user_header='X-SDF-User')
    with pytest.raises(ValueError, match='proxy user header'):
        sdf.catalog_auth.CatalogAuthorizer(
            access_database=_access_database(tmp_path)
        )
