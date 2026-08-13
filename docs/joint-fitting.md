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

The HD 113766 package is the first end-to-end joint rehearsal. HD 283750 plus
its wide DC8 companion HD 283750B is the mixed-family rehearsal: two fully
resolved observation scopes fitted explicitly with `phoenix_sol` and
`koester_wd`. The latter verifies mixed main-sequence/white-dwarf compilation,
sampling, retained ownership, analysis, and plotting. Because the 124-arcsec
pair has no shared measurements, its likelihood correctly factorizes; it is a
model-family and ownership test rather than a composite-photometry test.

The HD 283750 rehearsal also demonstrates that current atmosphere-grid `logg`
and flux normalization priors are independent. Its SED and radii are sensible,
but the fitted gravity medians imply unphysical masses for both stars. Do not
interpret broadband-only `logg` scientifically until SDF gains family-specific
stellar mass/radius priors or deliberately treats gravity as a nuisance
parameter.

Result-database writing and indexed sample generation are intentionally not
enabled yet.
Requesting database output for a joint package raises a clear error rather
than silently writing a single-target-shaped legacy result. The next milestone
is a rebuildable result index, followed by fit-selection commands and the
operator workflow.

The first result-index implementation can build a disposable SQLite query view
of all retained joint candidates below an export root:

```shell
sdf-index rebuild /path/to/export
sdf-index verify /path/to/export
```

The default database is `/path/to/export/sdf-results.sqlite`; use `--database`
to place it elsewhere. The index records packages, canonical target identities,
candidates, component ownership, evidence, input freshness, and display/default
selection. Fitted and derived scalar quantities are indexed per physical target
and owned model component. The `default_targets()` and `default_parameters()`
query APIs accept sample-sized lists of SIMBAD names, SDB IDs, or canonical
keys; these are the data boundary for migrating `sdf-sample`.

The index contains no irreplaceable state and can always be rebuilt from the
package-local `fits` directories. Once that default index exists, fitting a
package anywhere below its recorded root refreshes that package's rows
transactionally. SDF does not create an index implicitly. Reworking the sample
tables and plots to represent multiple physical targets and owned components is
the remaining `sdf-sample` integration step.

Retained fits can be inspected and presented without rerunning MultiNest:

```shell
sdf-fits list /path/to/package
sdf-fits select /path/to/package \
  --default FIT-REFERENCE --display FIT-REFERENCE [FIT-REFERENCE ...]
sdf-fits render /path/to/package
sdf-fits render /path/to/package --fit FIT-REFERENCE \
  --output /tmp/fit-preview.html
```

`list` prints the candidate storage names used as concise fit references,
evidence, and current selection state; `--json` provides the same data for a
UI. A full candidate ID or an unambiguous storage-name prefix can also be used
as a reference. Selecting a new default automatically puts it first in the
displayed tabs, regenerates `index.html`, and refreshes an existing result
index. It does not run a fit.

An explicit model combination can be fitted without changing the public
selection. Repeat `--disk` to request more than one disk component:

```shell
sdf-fits fit /path/to/package \
  --stellar 'HD 113766A=phoenix_sol' \
  --stellar 'HD 113766B=phoenix_sol' \
  --disk 'HD 113766A=modbb_disk_r' \
  --disk 'HD 113766A=modbb_disk_r'
```

Every physical contributor must have exactly one `--stellar` assignment.
Targets can be identified by canonical SIMBAD name, SDB ID, or the stored
identity key. The command reuses an identical current candidate unless
`--force` is supplied. To force a fresh MultiNest run for an already retained
candidate, use `sdf-fits refit PACKAGE FIT-REFERENCE`. Both fitting commands
remain synchronous; the operator service will put the same operations behind
the fit-job queue.

Run the regression suite with:

```shell
PYTHONDONTWRITEBYTECODE=1 conda run -n sdf \
  pytest -q -p no:cacheprovider
```
