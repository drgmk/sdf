# Joint-fitting development workflow

Joint fitting is currently being developed on `codex/joint-fitting`. Ordinary
single-file fitting remains on the established `Result` path. A directory is
recognized as one package when it contains `joint-fit.yml`, `spectra.yml`, or
rawphot files carrying the explicit `sdf_package='.'` keyword. Any rawphot
member resolves to that package and is routed through the new fitting core.

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

SDB package exports put this portable marker in every rawphot member:

```text
\sdf_package='.'
```

It distinguishes a package containing several self-mapped rawphot files from
an ordinary directory of unrelated inputs. A spectral package also declares
`spectra_manifest='spectra.yml'`; SDF requires the referenced manifest to be
present and requires all package members to agree. These pointers contain no
installation-specific paths.

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

An optional `spectra.yml` adds spectra without embedding installation-specific
paths in rawphot metadata:

```yaml
version: 1
spectra:
  - id: cassis-irs:14982144:preferred
    file: spectra/cassis_yaaar_spcfw_14982144t.fits
    format: cassis-irs
    contributors:
      - component-a-sdbid
      - component-b-sdbid
```

Files must be below the package directory. An optional `sha256` verifies the
asset before fitting. `contributors` are independent of the rawphot scopes,
so an unresolved spectrum can see A and B while resolved photometry remains
separate. All segments inherit the product contributors by default. Rare
module-specific cases can override ownership or eligibility:

```yaml
    segment_overrides:
      SL1:
        contributors: [component-a-sdbid]
      LL1:
        include: false
```

The CASSIS reader forms independently normalized SL/LL or SH/LH segment
observations. Generic `csv` spectra use wavelength (micron), flux density (Jy),
and uncertainty (Jy) as their first three columns. Several products may be
listed and all included products are fitted; there is no implicit first-file
selection. `--no-spectra` ignores both this manifest and legacy rawphot
spectrum keywords.

For segments with resolving-power metadata, model grids are convolved once at
fit compilation and then sampled at the observed wavelengths. The prepared
model is cached across candidates. A grid already coarser than the observation
is interpolated rather than spuriously sharpened and is marked under-resolved
in the preparation metadata.

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
`exhaustive` are runnable, and `sdf-fits fit` provides an explicit exact-model
CLI/API path. `posterior-only` remains represented in policy but does not yet
have a dedicated command.

### Evidence-guided priority

The configured evidence threshold remains 5. Evidence-guided searches fit
complete disk-count levels rather than following individual ownership branches.
For a binary, for example, every one-disk placement is evaluated before the
one-disk level is judged; if that level is supported, every two-disk placement
is then evaluated. Search stops after the best fit in a complete level fails to
improve on the best simpler level by more than the threshold. This deliberately
assumes that a jointly useful pair of added components will normally have a
useful one-component precursor.

After fitting, SDF finds the maximum evidence and treats candidates within 5
of it as statistically competitive (with evidence uncertainty handled
conservatively near the boundary). It recommends the competitive candidate
with the fewest disk components. Equivalent placements then use deterministic
preferences: concentrate disks around fewer physical contributors, prefer the
earlier hierarchy contributor (A before B), follow configured model-family
order, then prefer fewer numerical parameters. These rules never override an
evidence difference greater than the threshold.

`selection.json` records both the automatic `recommended` fit and the public
`default`. Choosing a default in the operator UI marks it as an operator
override, which subsequent automatic fits preserve. The UI continues to show
the recommendation separately and offers **Use recommended** to return control
to the automatic procedure.

Priority rules can be reapplied to retained fits without sampling again:

```shell
sdf-fits reprioritize /path/to/package
sdf-fits reprioritize /path/to/export --all
```

The bulk form uses `sdf-operator.sqlite` when present, rewrites only package
selection metadata, and rebuilds the disposable result index once at the end.
Candidates made from superseded package inputs are ignored and current
operator-selected defaults are preserved.

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
models, residuals, target inputs, fit products, and external resources. It does
not need the native MultiNest output. Physical contributors and observed
scopes are separate: resolved scopes share the contributor colour, while
composite scopes use a neutral system colour and plot the summed model.
The public layout keeps the SED above the fold in a legacy-sized central
column. A compact left rail projects composite observations as a hierarchy,
without displaying SDB identifiers, and the fit tabs use `★`/`○` shorthand
for stellar and disk components (for example `A ★○ · B ★`). Per-owner summary
blocks repeat the plot colours and carry the physical results and product
links. The plot keeps only a small in-panel legend for data, model photometry,
and upper limits; component and total-model lines stay out of it. The SED and
residual panels share one toolbar and x range, following the established SDF
axis layout.

The operator view is an extension of this same document rather than a second
page design. It adds a right-hand drawer while retaining the public content and
styling:

```shell
sdf-fits render /path/to/package --operator
```

This writes a read-only `operator.html` development view by default. For normal
operation, run one loopback-only service over the whole SDB export root:

```shell
sdf-fits serve /path/to/export --actor Grant
```

For a backlog, fit independent systems concurrently with a bounded worker
pool:

```shell
sdf-fits serve /path/to/export --actor Grant --workers 4 --autofit-existing
```

The catalog home page includes every currently exported fitting package,
including unfit systems. It provides target/SDB-ID search, fit-state filters,
sample views sourced from the latest complete sample-view snapshot in the SDB
export manifests (independent of whether that run selected one target, a
sample, or all targets), a global
job list, and links to the shared public/operator SED page for each fitted
system. An unfit system instead gets a low-friction initial-fit action with a
pre-filled reason.

In deployed use this operational catalog now lives below `/admin`, while the
root page and individual `/system/...` SEDs remain public. Authenticated sample
views live below `/samples` and retain the legacy per-sample group semantics.
Apache can reuse the existing `.htpasswd` and `.htgroup` files; see
[`catalog-auth.md`](catalog-auth.md) for the route model, CLI options, and a
reverse-proxy example. Local loopback use without those options retains the
low-friction operator workflow.

By default the same process supervises one isolated fit worker and watches the
export tree for new or changed rawphot/topology inputs. `--workers N` starts N
isolated worker processes, each fitting one system at a time. Atomic queue
claims allow those workers to share the operator database, while package-level
exclusion prevents an older fit and a newly requested refit of the same system
from overlapping. Start conservatively because each MultiNest process has its
own model arrays and memory footprint; the useful value is normally bounded by
both CPU cores and available RAM.

New and stale packages are indexed and queued automatically; the existing
unfit backlog is left alone unless `--autofit-existing` is supplied. Useful
development controls are `--workers`, `--no-worker`, `--no-watch`,
`--rebuild-index`, and `--poll-seconds`. Interrupted `running` jobs return to
the queue at service startup. If an individual pool worker fails, only that
worker's job is recovered before its slot is restarted. The result index
defaults to `/path/to/export/sdf-results.sqlite` and durable policy/job state
to `/path/to/export/sdf-operator.sqlite`. `SIGTERM` follows the same orderly
shutdown path as an interactive interrupt, allowing supervised fit workers to
stop before the service exits.

The drawer exposes model choices for each physical contributor, but actor,
reason, and the queue action apply once to the system. Changed contributor
configs are saved transactionally and produce one system fit job. A custom
config is shared wherever that physical target is used unless the export
already needs a system-specific exception. The retained default fit can seed
an otherwise-missing config when the service first opens; this preserves an
explicit white-dwarf or cool-star choice even when the exported rawphot
metadata lacks a spectral type. Choosing the default records that choice so a
later service restart does not reconstruct the removed custom config.

The service also lists every retained candidate. `Preview` does not change the
public selection; `shown` controls the public tabs, and `Make default` changes
the default fit immediately without running MultiNest.

The standalone job command remains useful for inspection or headless worker
operation when the catalog service is started with `--no-worker`:

```shell
sdf-jobs --state /path/to/sdf-operator.sqlite list --status queued
sdf-jobs --state /path/to/sdf-operator.sqlite run
```

Use `run --once` for a single queued job. The operator database is intentionally
separate from the rebuildable result index, so rebuilding result queries cannot
erase policy or job state.

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

Legacy result-database writing remains intentionally disabled for joint
packages. Requesting it raises a clear error rather than silently writing a
single-target-shaped result. The replacement result index, fit-selection
commands, persistent policy, catalog-wide loopback service, supervised fit-job
queue, export watcher, and hierarchy-aware sample products are available. SDF
does not yet write photometry decisions or other edits back to SDB.

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
and owned model component. The `default_targets()`, `default_parameters()`, and
`default_packages()` query APIs accept sample-sized lists of SIMBAD names, SDB
IDs, canonical keys, or composite observation IDs. Observation rows retain
their ordered physical contributors, so a sample entry for an unresolved `AB`
measurement can find the package and expose both A and B without treating
either as a direct sample member.

The index contains no irreplaceable state and can always be rebuilt from the
package-local `fits` directories. Once that default index exists, fitting a
package anywhere below its recorded root refreshes that package's rows
transactionally. SDF does not create an index implicitly; rebuild it when SDF
reports that its disposable schema is out of date. Hierarchy-aware sample
products are available through an explicit indexed path:

```shell
sdf-sample --index /path/to/export/sdf-results.sqlite \
  --samples reasons --tables --plots \
  --output-root /path/to/www/samples \
  --sed-base-url https://example.test/seds
```

Sample membership still comes from SDB. Fitted values come from the default
candidate in the disposable SDF index. The generated table has one row per
physical contributor, marks direct sample members separately from system
context, and keeps disks under their physical owner. Its VOTable follows the
same projection. The long-form photometry CSV has one row per observed datum
and fitted component, with explicit observation contributors, component owner,
flags, observed flux, total-model posterior interval, and component-model
posterior interval. This avoids variable component columns and remains usable
for single stars, resolved components, and composite measurements.

The plots provide a hierarchy-aware HR diagram, owned-disk fractional
luminosity versus temperature, and disk flux/size planning plots in the legacy
MIRI, PACS 70, and 880-micron bands. The latter convolve the retained median
component spectra with SDF filters on demand, so arbitrary-band products do not
inflate candidate artifacts. Disk/star contrast uses stellar spectra belonging
to the disk's physical owner. Plot points link back to their retained SED page.

Without `--index`, `sdf-sample` retains the old MySQL-shaped implementation.
The indexed path now replaces the summary table, VOTable, photometry CSV, HR
and disk plots, and flux/size plots without reconstructing old `Result` objects
or querying legacy result tables.

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
