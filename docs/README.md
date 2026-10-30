# `sdf` documentation

`sdf`'s documentation is organized by what a reader is trying to do. Public web
features and administrative fitting controls deliberately have separate
guides: most visitors should not need to understand model search, fit queues,
or the `sdb` export layout.

## Website users

- [Using the `sdf` website](web-guide.md) explains search, SED pages, fit tabs,
  sample tables and plots, downloads, sign-in, and user-created samples.
- [Using `sdf` fit outputs](fit-outputs.md) explains the downloadable JSON and
  NPZ products with Python examples.

## Fit administrators

- [Fit administration](admin-guide.md) covers the administrative dashboard,
  data review, model configuration, fit jobs, and choosing the public fit.

## Server operators

- [Server and fitting operations](server-operations.md) covers
  `sdf-fits serve`, worker and watcher behavior, state files, `sdb` integration,
  and maintenance commands.
- [Authentication and authorization](catalog-auth.md) describes the trusted
  reverse-proxy boundary, accounts, roles, sessions, and sample access.
- The separate
  [`sdb-server`](https://github.com/drgmk/sdb-server) repository contains the
  combined Apache/`sdf`/`sdb` deployment scripts.

## Fit and data-model reference

- [Joint fitting reference](joint-fitting.md) defines targets, observation
  scopes, package manifests, spectra, model policy, evidence selection, and
  retained candidates.
- [Python API](api.rst) and [command-line tools](user.rst) are included in the
  Sphinx documentation.

## Historical design material

- [Joint-fitting architecture plan](history/joint-fitting-plan.md) records the
  design process. It is not an operating guide; current behavior is documented
  in the references above.

Build the combined Markdown and API documentation with:

```shell
python -m pip install '.[docs]'
make -C docs html
```
