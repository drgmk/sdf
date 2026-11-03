# Catalog authentication and authorization

This is a server-operator reference. Website behavior for ordinary accounts is
described in [Using the `sdf` website](web-guide.md), and administrator fitting
controls in [Fit administration](admin-guide.md). The recommended combined
proxy configuration is maintained in the separate
[`sdb-server`](https://github.com/drgmk/sdb-server) deployment repository.

`sdf-fits serve` can sit behind an HTTPS reverse proxy while retaining an
existing Apache `.htpasswd` file. Apache verifies the password only at the
login handoff; `sdf` never reads or stores passwords. Users, administrator roles,
sample ownership, visibility, and explicit sample grants live in the
authoritative `sdb` SQLite database.

A successful login creates a signed, expiring session shared by the `sdf` catalog
and `sdb` review applications. `.htgroup`, per-sample `.htaccess`, and a separate
public-sample allowlist are not used.

The route boundary is:

| Route | Access |
| --- | --- |
| `/`, `/search`, `/target/...`, `/system/...` | public |
| `/samples` | public samples, plus owned/shared samples when signed in |
| `/samples/NAME` | public, owner, explicit grantee, or administrator |
| `/samples/new`, `/sample-imports/...` | any signed-in user; each import is private by default |
| `/sample-edits` | owning user for a private sample, or an administrator |
| `/auth/login` | Apache-authenticated handoff that creates the signed session |
| `/imports`, `/imports/ID` | any signed-in user; jobs are visible to their requester and administrators |
| `/admin/...`, `/review/...` | `sdb` users with `is_admin` enabled |

Register accounts and manage roles with `sdb`. Usernames must match `.htpasswd`:

```shell
sdb --database /srv/sdb/sdb.sqlite user add grant@example.com --admin
sdb --database /srv/sdb/sdb.sqlite user add collaborator@example.com
sdb --database /srv/sdb/sdb.sqlite user set-admin collaborator@example.com
sdb --database /srv/sdb/sdb.sqlite user remove-admin collaborator@example.com
sdb --database /srv/sdb/sdb.sqlite user disable collaborator@example.com
```

Sample commands can set access without changing membership:

```shell
sdb --database /srv/sdb/sdb.sqlite sample set-access reasons \
  --owner grant@example.com --visibility private
sdb --database /srv/sdb/sdb.sqlite sample share reasons collaborator@example.com
sdb --database /srv/sdb/sdb.sqlite sample unshare reasons collaborator@example.com
sdb --database /srv/sdb/sdb.sqlite sample delete reasons --yes \
  --actor grant@example.com --reason 'Remove obsolete sample'
```

Deleting a sample removes the grouping, membership history, and grants but not
the targets or their data. `sdb` retains a compact deletion audit record.

Start both applications with the same stable secret (at least 32 bytes):

```shell
sdf-fits serve /srv/sdb/export \
  --workers 4 \
  --auth-user-header X-SDF-User \
  --auth-database /srv/sdb/sdb.sqlite \
  --auth-session-secret-file /srv/sdb/session-secret \
  --review-base-url /review

sdb --database /srv/sdb/sdb.sqlite review serve --all \
  --host 127.0.0.1 --port 8782 --base-path /review \
  --auth-user-header X-SDF-User \
  --auth-session-secret-file /srv/sdb/session-secret \
  --auth-login-url /auth/login
```

Without proxy-auth options, `sdf` remains in local development mode and treats
the CLI actor as an administrator.

Apache needs to authenticate only `/auth/login`; application code authorizes
all other routes from the signed session and current database state:

```apache
# Required modules: auth_basic, authn_file, headers, proxy, and proxy_http.
RequestHeader unset X-SDF-User early

<Location "/auth/login">
    AuthType Basic
    AuthName "sdb"
    AuthUserFile /srv/sdb/.htpasswd
    Require valid-user
    RequestHeader set X-SDF-User "expr=%{REMOTE_USER}"
</Location>

ProxyPass        /review/ http://127.0.0.1:8782/
ProxyPassReverse /review/ http://127.0.0.1:8782/
ProxyPass        / http://127.0.0.1:8765/
ProxyPassReverse / http://127.0.0.1:8765/
```

With HTTPS, add `--auth-cookie-secure`. Keep `.htpasswd` and the session secret
outside the document root. `sdf` refuses non-loopback bind addresses, so the
trusted identity header is accepted only from the local reverse proxy.

The authenticated username is also the audit actor. POST bodies cannot
override it, and page tokens provide cross-site request-forgery protection.

Signed-in users can import a missing SIMBAD target from search or create a
private sample from up to 1000 text/CSV/TSV identifiers. Sample acquisition is
preceded by a local check that distinguishes existing targets, new identifiers
that require SIMBAD, and ambiguous aliases. Existing targets are reused without
remote name resolution. Acquisition is durable and uses the same single `sdb`
writer as individual imports; successful targets are exported and fitted even
when some identifiers fail. Only administrators can make a newly imported
sample public. Public individual SED pages remain public independently of
sample access.

An owner may add or remove members of their private sample from its Properties
page. Explicit grantees have read access only. Administrators may edit all
samples and are the only users permitted to change a public sample. Every
membership change is recorded by `sdb`; the asynchronous job then republishes
the sample snapshot before the web view changes.

Public system routes never serve `operator.html`. Live configuration pages are
rendered ephemerally, avoiding persistence of the page token or authenticated
username in public package directories.
