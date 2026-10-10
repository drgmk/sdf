# Catalog authentication and authorization

`sdf-fits serve` can sit behind an HTTPS reverse proxy while preserving the
legacy SDB `.htpasswd` and `.htgroup` files. Apache verifies passwords; SDF
never reads or stores them. Apache passes the authenticated username in a
header, and SDF reads `.htgroup` to authorize sample and fit-configuration
views.

The route boundary is:

| Route | Access |
| --- | --- |
| `/`, `/search`, `/target/...`, `/system/...` | public |
| `/samples` | authenticated, filtered to permitted samples |
| `/samples/NAME` | group `NAME` or the administrator group |
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

Start the loopback-only application with explicit proxy authentication:

```shell
sdf-fits serve /srv/sdb/export \
  --workers 4 \
  --auth-user-header X-SDF-User \
  --auth-group-file /srv/sdb/auth/.htgroup \
  --auth-admin-group admin
```

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

<Location "/samples">
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

ProxyPass        / http://127.0.0.1:8765/
ProxyPassReverse / http://127.0.0.1:8765/
```

Keep the password and group files outside the document root. SDF continues to
refuse non-loopback bind addresses, so the trusted identity header is only
accepted from a local reverse proxy. The proxy should be the only network path
to the service.

The authenticated username is the operator actor. Proxy-authenticated pages
therefore do not offer an editable actor field, and POST bodies cannot override
the proxy identity. The page token remains useful for cross-site request
forgery protection but is not an authentication mechanism.

Public system routes never serve `operator.html`. Live operator pages are
rendered ephemerally, avoiding persistence of the page token or authenticated
username inside a public package directory.
