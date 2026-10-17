# Catalog authentication and authorization

`sdf-fits serve` can sit behind an HTTPS reverse proxy while preserving the
legacy SDB `.htpasswd` and `.htgroup` files. Apache verifies passwords; SDF
never reads or stores them. Apache passes the authenticated username in a
header, and SDF reads `.htgroup` to authorize sample and fit-configuration
views. A successful sample sign-in creates a signed, expiring SDF session so
the public and user-specific samples can share one page and one URL tree.

The route boundary is:

| Route | Access |
| --- | --- |
| `/`, `/search`, `/target/...`, `/system/...` | public |
| `/samples` | public samples, plus permitted private samples when signed in |
| `/samples/NAME` | public allowlist, group `NAME`, or administrator group |
| `/auth/login` | Apache-authenticated handoff back to `/samples` |
| `/imports`, `/imports/ID` | any signed-in user; jobs are visible to their requester and administrators |
| `/admin`, `/admin/jobs`, `/admin/system/...` | administrator group |
| `/admin/system/.../api/...` | administrator group plus the page CSRF token |

Group files use the existing Apache syntax:

```text
admin: grant alice
reasons: grant bob
arks: alice
```

Group names must match exported sample names. The configured administrator
group can access every sample. Changes to `.htgroup` are picked up without
restarting SDF.

Public access is fail-closed and explicit. Put one exact exported sample name
per line in a separate file; blank lines and lines beginning with `#` are
ignored:

```text
# Available without signing in
reasons
DEBRIS
```

Changes to this file are also picked up without restarting SDF. A missing or
empty public-samples file does not make private samples public.

Start the loopback-only application with explicit proxy authentication:

```shell
sdf-fits serve /srv/sdb/export \
  --workers 4 \
  --auth-user-header X-SDF-User \
  --auth-group-file /srv/sdb/auth/.htgroup \
  --public-samples-file /srv/sdb/auth/public-samples \
  --auth-admin-group admin
```

When the SDB review application is mounted on the same origin, add
`--review-base-url /review`. Administrator SED pages then link directly to the
corresponding SDB review target. SDB can link back with
`sdb review serve --sdf-base-url ''`, where the empty URL denotes the shared
origin.

`--auth-user-header` and `--auth-group-file` must be supplied together. With
neither option, SDF remains in local development mode and treats the CLI actor
as an administrator.

An Apache 2.4 virtual host can authenticate only the protected route trees and
proxy everything to SDF:

```apache
# Required modules: auth_basic, authn_file, authz_groupfile, headers,
# proxy, proxy_http and ssl.

# Never allow a browser to assert the trusted application identity.
RequestHeader unset X-SDF-User early

<Location "/auth/login">
    AuthType Basic
    AuthName "SDB samples"
    AuthUserFile /srv/sdb/auth/.htpasswd
    Require valid-user
    RequestHeader set X-SDF-User "expr=%{REMOTE_USER}"
</Location>

<Location "/admin">
    AuthType Basic
    AuthName "SDB administration"
    AuthUserFile /srv/sdb/auth/.htpasswd
    AuthGroupFile /srv/sdb/auth/.htgroup
    Require group admin
    RequestHeader set X-SDF-User "expr=%{REMOTE_USER}"
</Location>

<Location "/review">
    AuthType Basic
    AuthName "SDB administration"
    AuthUserFile /srv/sdb/auth/.htpasswd
    AuthGroupFile /srv/sdb/auth/.htgroup
    Require group admin
    RequestHeader set X-SDF-User "expr=%{REMOTE_USER}"
</Location>

ProxyPass        /review/ http://127.0.0.1:8782/
ProxyPassReverse /review/ http://127.0.0.1:8782/

ProxyPass        / http://127.0.0.1:8765/
ProxyPassReverse / http://127.0.0.1:8765/
```

With HTTPS, add `--auth-cookie-secure` so browsers never send the sample login
session over a plaintext connection. `/samples` itself remains public: an
anonymous visitor sees the explicit public allowlist, while a signed-in user
sees the union of public samples and their current `.htgroup` memberships.
Private direct links send anonymous visitors through `/auth/login` and return
them to the original sample page.

Keep the password and group files outside the document root. SDF continues to
refuse non-loopback bind addresses, so the trusted identity header is only
accepted from a local reverse proxy. The proxy should be the only network path
to the service.

The authenticated username is the operator actor. Proxy-authenticated pages
therefore do not offer an editable actor field, and POST bodies cannot override
the proxy identity. The page token remains useful for cross-site request
forgery protection but is not an authentication mechanism.

When SDB import is configured, an unsuccessful public search offers anonymous
visitors a sign-in link. Any signed-in account can then request an exact SIMBAD
target; the username is retained on the durable import record. This is broader
than fit configuration, which remains administrator-only, but avoids allowing
anonymous writes to the authoritative SDB database.

Public system routes never serve `operator.html`. Live operator pages are
rendered ephemerally, avoiding persistence of the page token or authenticated
username inside a public package directory.
