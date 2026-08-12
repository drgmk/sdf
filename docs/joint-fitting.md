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
parallax.

Joint HTML/plot generation and result-database writing are intentionally not
enabled yet. They require the remaining hierarchy-aware presentation/index
work; requesting those legacy outputs for a joint package raises a clear
error rather than silently producing single-target results.

Run the regression suite with:

```shell
PYTHONDONTWRITEBYTECODE=1 conda run -n sdf \
  pytest -q -p no:cacheprovider
```
