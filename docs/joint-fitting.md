# Joint-fitting development workflow

Joint fitting is currently being developed on `codex/joint-fitting`. Ordinary
single-file fitting remains on the established `Result` path. A directory
containing `joint-fit.yml`, or any rawphot file inside that directory, is
recognized as one package and routed through the new fitting core.

Install the branch into the `sdf` environment after changing source:

```shell
conda run -n sdf python -m pip install .
```

Fit an exported package with:

```shell
sdf-fit --file /path/to/package
```

Passing several rawphot members of the same joint package still runs the
package once. It uses a package-level `.sdf_lock`; unrelated rawphot files
retain their per-file locks and existing output behavior.

## Inputs

The topology file remains deliberately small:

```yaml
version: 1
observations:
  system-sdbid:
    - component-a-sdbid
    - component-b-sdbid
```

Omitted observations self-map. The SDB IDs only connect files inside the
package. SDF resolves model ownership and saved policy to canonical SIMBAD
names when loading it.

An optional `sdf-fit.yml` applies portable model policy by SIMBAD name:

```yaml
targets:
  HD 48915 A:
    stellar_models: [phoenix_sol]
    disk_models: [modbb_disk_r]
    max_disks: 2
    search_mode: evidence-guided
  HD 48915 B:
    stellar_models: [koester_wd]
    disk_models: [bb_disk_r]
    max_disks: 1
    search_mode: evidence-guided
```

Without this file, stellar grids are inferred independently from each physical
target's spectral type. Disk support is never inferred from stellar class;
`max_disks` is a search limit, not a requirement to fit disks.

The whole package currently needs one search mode. `evidence-guided` and
`exhaustive` are runnable. `exact` and `posterior-only` are represented in
policy but still need an explicit candidate CLI/API input.

## Outputs

Every completed candidate is committed before the search continues:

```text
package/
  fits/
    selection.json
    <readable-model-name>--<identity-hash>/
      result.json
      posterior-<content-hash>.npz
    .work/
      ... sampler-native resumable output ...
```

`result.json` contains identities, component ownership, evidence, material fit
state, code provenance, policy, posterior summaries, observation scopes, and
target-owned stellar/disk results. The compressed array contains capped
posterior samples plus observation and component predictions. Re-fitting the
same exact candidate replaces its artifact; it does not create dated result
history. `selection.json` independently controls the default and displayed
candidates, so changing presentation does not refit anything.

Changing rawphot/topology makes retained candidates non-current. Changes to
models or material sampler settings also make them non-current. Merely moving
the package, changing display selection, or changing SDF code does not; the
code revision and dirty flag are provenance.

## Current boundary

The branch can discover, compile, sample, retain, and reuse joint candidates.
It supports resolved/composite photometry, colours, upper limits, and spectra.
Analysis uses each stellar or disk component's owner and each target's own
parallax. Every joint-package run also writes a standalone Bokeh `index.html`
from the compact retained candidates. The page follows `fits/selection.json`,
shows the default fit first, and includes scoped photometry, component/system
models, and residuals. It does not need the native MultiNest output.

The HD 113766 package is the first end-to-end joint rehearsal. Result-database
writing and indexed sample generation are intentionally not enabled yet.
Requesting database output for a joint package raises a clear error rather
than silently writing a single-target-shaped legacy result. The next milestone
is a rebuildable result index, followed by fit-selection commands and the
operator workflow. A normal-star/white-dwarf package remains the second
hierarchy-aware validation case.

The first result-index implementation can build a disposable SQLite query view
of all retained joint candidates below an export root:

```shell
sdf-index rebuild /path/to/export
sdf-index verify /path/to/export
```

The default database is `/path/to/export/sdf-results.sqlite`; use `--database`
to place it elsewhere. The index records packages, canonical target identities,
candidates, component ownership, evidence, input freshness, and display/default
selection. It contains no irreplaceable state and can always be rebuilt from
the package-local `fits` directories. Automatic index updates, fit-selection
commands, and `sdf-sample` integration are the next implementation slice.

Run the regression suite with:

```shell
PYTHONDONTWRITEBYTECODE=1 conda run -n sdf \
  pytest -q -p no:cacheprovider
```
