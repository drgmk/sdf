# `sdf`

all your seds are belong to us

`sdf` fits stellar and circumstellar-disk models to photometry and spectra. It
supports single targets and hierarchical systems whose resolved and composite
measurements must be fitted together.

## Install

Install the working tree into an environment that provides MultiNest:

```shell
python -m pip install .
```

The historical `sdf.yml` environment remains in the repository as a reference,
but development installations use normal `pip` packaging.

## Configure

Configuration is optional when the bundled models are sufficient. For a full
model library, copy `sdf.example.toml` to `sdf.toml` and set its
`[models].root`. SDF also discovers `~/.config/sdf/config.toml`, or accepts an
explicit path through `SDF_CONFIG`. Relative paths are resolved from the TOML
file. See [the configuration guide](docs/config.rst) for the small set of
process-wide fitting controls; per-system model choices remain with the fit
package or web UI.

To build the documentation, install the optional documentation dependencies:

```shell
python -m pip install '.[docs]'
make -C docs html
```

## Documentation

Start with [the documentation guide](docs/README.md). The principal entry
points are:

- [Using the `sdf` website](docs/web-guide.md)
- [Fit administration](docs/admin-guide.md)
- [Server and fitting operations](docs/server-operations.md)
- [Joint fitting reference](docs/joint-fitting.md)
- [Using JSON and posterior outputs](docs/fit-outputs.md)

`sdf-fits serve` is the current catalog-wide website and fitting service.
Run `sdf-fits --help` for its package-management commands and
`sdf-fits serve --help` for service options.

For a standalone IPAC photometry table, `sdf-fit TABLE` creates and fits a
self-contained result package without starting the service.
