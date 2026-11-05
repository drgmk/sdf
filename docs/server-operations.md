# Server and fitting operations

This guide covers running the catalog-wide website and fit service. The
recommended combined deployment is maintained in the separate
[`sdb-server`](https://github.com/drgmk/sdb-server) repository; it starts
Apache, `sdf-fits serve`, and `sdb review serve` with one configuration and
shared authentication. This document describes the `sdf` side of that stack.

## The service boundary

One `sdf-fits serve` process operates over a complete `sdb` export root:

```shell
sdf-fits serve /srv/sdb/export --workers 2
```

It provides:

- the public search, SED, and sample website;
- the administrator dashboard and SED fit controls;
- an administrator activity page combining catalogue inventory, audited
  changes, queue health, and optional aggregate access-log statistics;
- a rebuildable result index;
- durable fit policy, selection, and job state;
- one or more isolated fit workers; and
- a watcher for new or changed `sdb` exports.

The service binds only to a loopback address. Publish it through the trusted
reverse proxy rather than binding it directly to an external interface.

For local development without proxy authentication:

```shell
sdf-fits serve /path/to/export --actor Grant --port 8765
```

Local mode treats the named actor as an administrator. Production should use
the authenticated configuration in
[Authentication and authorization](catalog-auth.md).

For local page-use statistics, give the service the frontend access log:

```shell
sdf-fits serve /srv/sdb/export \
  --access-log /srv/sdb-server/logs/apache-access.log \
  ...
```

The activity page parses both Apache common logs and the extended format used
by the reference `sdb-server` deployment. Client addresses are discarded while
the report is built. The extended format adds user-agent and request-duration
fields, allowing known crawlers to be excluded and median request time to be
reported.

## Export-root state

By default `sdf` creates two SQLite files beside the exported packages:

| File | Purpose | Recovery |
| --- | --- | --- |
| `sdf-results.sqlite` | Fast projection of package-local retained results | Rebuild from `fits/` artifacts |
| `sdf-operator.sqlite` | Policies, selections, fit/import jobs, and audit state | Back up; this is durable state |

The result index is disposable. The administrative database is not: include it
in normal server backups together with the authoritative `sdb` database and
private authentication material.

Package-local numerical results remain under each package's `fits/` directory.
Large sampler-native output may exist temporarily in `fits/.work/` but is not
required for normal result display after a candidate has been committed.

## Workers and export watcher

The default service supervises one worker. Increase concurrency only for
independent systems:

```shell
sdf-fits serve /srv/sdb/export --workers 4
```

Each worker is a separate process and each active MultiNest fit can consume
substantial memory. Choose the worker count from available RAM as well as CPU
cores. On the current 10-core, 16 GB development host, two workers are the
conservative interactive setting and four is a reasonable starting point when
the machine is dedicated to the service. Monitor representative simultaneous
fits before increasing further.

Atomic queue claims allow workers to share one administrative database. A
package-level lock prevents overlapping fits of the same system. If a worker
dies, only its claimed job is recovered and its slot is restarted.

Only one worker-owning service may use an operator database at a time. A
second `sdf-fits serve` process fails at startup instead of recovering or
duplicating jobs owned by the first. Additional read-serving processes must
use `--no-worker`.

The export watcher indexes new and materially changed packages and queues them
according to their state. It does not queue the entire pre-existing unfit
backlog unless requested:

```shell
sdf-fits serve /srv/sdb/export --autofit-existing
```

Useful controls are:

- `--workers N`: number of concurrent isolated fit workers;
- `--no-worker`: serve and queue without consuming jobs;
- `--no-watch`: disable export-tree polling;
- `--autofit-existing`: queue the existing unfit backlog on startup;
- `--poll-seconds N`: worker/watcher polling interval; and
- `--rebuild-index`: rebuild the result index before serving.

Interrupted running jobs are returned to the queue at startup. `SIGTERM` and
an interactive interrupt both request an orderly worker shutdown.

## `sdb` integration and imports

Target and sample import require the authoritative `sdb` paths:

```shell
sdf-fits serve /srv/sdb/export \
  --sdb-executable /opt/sdb/bin/sdb \
  --sdb-database /srv/sdb/sdb.sqlite \
  --sdb-reference-database /srv/sdb/sdb-reference.sqlite \
  --sdb-cache-database /srv/sdb/sdb-cache.sqlite \
  --sdb-config /srv/sdb/sdb.toml
```

Imports are disabled unless `--sdb-database` is supplied. Acquisition uses a
single `sdb` writer, exports complete packages under a shared lock, updates the
`sdf` result index, and queues needed first fits. The shared lock prevents a fit
from seeing a package halfway through export.

Pass `--review-base-url /review` when the `sdb` review application is mounted on
the same site. Administrative SED pages then link directly to the corresponding
review target.

## Authentication

A production `sdf` service accepts identity only from a trusted loopback reverse
proxy and checks users/roles in the authoritative `sdb` database:

```shell
sdf-fits serve /srv/sdb/export \
  --auth-user-header X-SDF-User \
  --auth-database /srv/sdb/sdb.sqlite \
  --auth-session-secret-file /srv/sdb/session-secret \
  --auth-cookie-secure
```

If `--auth-database` is omitted it defaults to `--sdb-database`. Both `sdf` and
`sdb` review must use the same stable session secret. See
[Authentication and authorization](catalog-auth.md) for the complete route and
proxy model.

## Queue diagnostics and headless workers

The web Jobs page is normally sufficient. `sdf-jobs` is useful for diagnostics
or when the catalog service runs with `--no-worker`:

```shell
sdf-jobs --state /srv/sdb/export/sdf-operator.sqlite list
sdf-jobs --state /srv/sdb/export/sdf-operator.sqlite list --status failed
sdf-jobs --state /srv/sdb/export/sdf-operator.sqlite run --once
sdf-jobs --state /srv/sdb/export/sdf-operator.sqlite run
```

Do not start an independent consumer pool against the same state database
without coordinating it with the service's supervised workers.

## Result-index maintenance

The service updates the result index as packages change. It can always be
rebuilt or checked explicitly:

```shell
sdf-index rebuild /srv/sdb/export
sdf-index verify /srv/sdb/export
```

Rebuilding does not alter retained fit artifacts or the administrative
database.

Reapply current evidence-priority rules without sampling:

```shell
sdf-fits reprioritize /srv/sdb/export --all
```

Refresh older posterior NPZ metadata without sampling:

```shell
sdf-fits refresh-posteriors /srv/sdb/export --all
```

## Package-level commands

Most routine work belongs in the web interface. These commands remain useful
for inspection, automation, and recovery:

```shell
sdf-fits list PACKAGE
sdf-fits list PACKAGE --json
sdf-fits render PACKAGE
sdf-fits select PACKAGE --default FIT --display FIT OTHER_FIT
sdf-fits refit PACKAGE FIT
```

`refit` and `fit --force` always use a new MultiNest work directory. They are
therefore genuine from-scratch runs rather than resumes of retained native
files. Ordinary automatic jobs continue to reuse a current validated
candidate.

Fit one exact owned component combination synchronously:

```shell
sdf-fits fit PACKAGE \
  --stellar 'HD 113766A=phoenix_sol' \
  --stellar 'HD 113766B=phoenix_sol' \
  --disk 'HD 113766A=modbb_disk_r'
```

Every physical contributor needs exactly one `--stellar` assignment. Repeat
`--disk` for multiple disk components. An identical current candidate is
reused unless `--force` is supplied.

## Legacy commands

`sdf-fit` remains the compatibility entry point. A package directory is routed
through the joint-fitting core; an ordinary rawphot file still uses the legacy
single-file `Result` workflow. New catalog-wide operation should use
`sdf-fits serve`.

`sdf-sample --index ...` generates hierarchy-aware products from the current
result index. The live service renders equivalent sample pages directly. The
form without `--index`, legacy result-database writing, and `sdf-cleanup` belong
to the earlier MySQL/static-site workflow and should not be used as the model
for new operations.

Run `COMMAND --help` for the authoritative option list.
