# Joint fitting reference

SDF can fit all resolved and composite observations of a hierarchical system
in one likelihood. The central distinction is between a physical contributor
that owns a model and an observation scope that sees one or more contributors.

This document defines the package and fitting behavior. Website use is covered
by [Using the SDF website](web-guide.md), fit decisions by
[Fit administration](admin-guide.md), and process management by
[Server and fitting operations](server-operations.md).

## Physical contributors and observation scopes

A physical contributor owns stellar and optional disk model components. An
observation scope describes what an instrument measured. In an A/B binary with
resolved A and B photometry plus unresolved AB photometry:

| Observation scope | Model contributors |
| --- | --- |
| A | A |
| B | B |
| AB | A + B |

The AB observation does not own a third star. Its predicted flux is the sum of
the components owned by A and B. This separation is also used for spectra and
for catalogues whose beam blends several physical objects.

Every fitted parameter name is qualified by its owning component. Derived
stellar and disk quantities use that owner's astrometry and stellar
luminosity.

## Discovering a fitting package

A directory is an explicit fitting package when it contains `joint-fit.yml`,
`spectra.yml`, or rawphot members carrying:

```text
\sdf_package='.'
```

The marker lets a directory contain several self-mapped rawphot files without
requiring a non-trivial topology file. Any rawphot member resolves to the
package and the package is fitted only once.

An ordinary single rawphot file without a package marker remains supported by
the legacy `Result` path. Package inputs use the joint-fitting core and compact
candidate artifacts.

## Rawphot and topology

Each physical or composite exported target has its own rawphot file. The
optional `joint-fit.yml` records only non-default observation mappings:

```yaml
version: 1
observations:
  system-sdbid:
    - component-a-sdbid
    - component-b-sdbid
```

An observation omitted from the file maps to itself. The identifiers connect
files within the package; they are not persistent fit-policy identities. SDF
resolves physical ownership and saved configuration to canonical SIMBAD
identities where possible.

The package loader validates that every contributor exists, that mappings are
unambiguous, and that package members agree about their sidecar manifests.

## Spectra

Rawphot members can point to a portable package-relative `spectra.yml`:

```text
\spectra_manifest='spectra.yml'
```

The manifest declares products and their physical contributors:

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

Files must remain below the package directory. An optional `sha256` verifies
an asset before fitting. Several products may be listed and all included
products are fitted.

CASSIS IRS products are divided into their independently normalized SL/LL or
SH/LH segments. Every segment inherits the product contributor list unless a
rare segment-specific override is supplied:

```yaml
    segment_overrides:
      SL1:
        contributors: [component-a-sdbid]
      LL1:
        include: false
```

Generic `csv` spectra use wavelength in microns, flux density in Jy, and
uncertainty in Jy as their first three columns.

For products with resolving-power metadata, numerical model grids are
convolved once at fit compilation and sampled at the observed wavelengths. The
prepared model is reused across candidate fits. A model already coarser than
the observation is interpolated and marked under-resolved rather than
spuriously sharpened.

`--no-spectra` excludes both manifest spectra and legacy rawphot spectrum
keywords for diagnostic fits.

## External resources

An optional `resources.yml` attaches provider-neutral links to the relevant
observation scope without making those resources numerical fit inputs:

```yaml
version: 1
resources:
  - provider: tldr_pacs
    label: Herschel/PACS images
    url: https://example.org/herschel/sources/source-key.html
    observation_ids: ['1342223848', '1342223850']
    observation_sdbids: [system-sdbid]
    target_sdbids: [component-a-sdbid, component-b-sdbid]
```

A composite product is displayed once on its parent observation rather than
repeated on every contributor. Resource links do not affect fit freshness.

## Model configuration

An optional package-local `sdf-fit.yml` supplies portable configuration by
canonical SIMBAD name:

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

Without an explicit configuration, SDF infers a stellar family independently
for every physical contributor from its metadata. Disk availability is not
inferred from stellar class: white dwarfs and other unusual stars can still
host disks. `max_disks` is a search limit, not a required disk count.

Configuration precedence is:

1. an exact CLI or API request;
2. package-local `sdf-fit.yml`;
3. saved administrative configuration;
4. compatible legacy single-target configuration; and
5. automatic inference from exported target metadata.

The resolved configuration and its origin are copied into every result.
Administrative configuration is keyed by external physical identity rather
than generated SDB IDs, with package-specific state used only when necessary.

## Candidates and search

A candidate is one exact owned model combination: one stellar component per
physical contributor plus zero or more disks assigned to particular owners.
For a binary, a one-disk candidate around A is distinct from the otherwise
identical candidate with the disk around B.

The whole package uses one search mode:

- `evidence-guided` evaluates complete complexity levels and stops when the
  next level is unsupported;
- `exhaustive` evaluates every configuration permitted by the limits; and
- an exact `sdf-fits fit` request evaluates only the supplied owned model
  combination.

Evidence-guided search fits all placements at a disk-count level before
judging the level. For example, every one-disk placement in a binary is
evaluated before SDF decides whether to try two disks. The configured
log-evidence threshold is 5.

After fitting, candidates within 5 of the maximum evidence are treated as
competitive, with evidence uncertainty handled conservatively near the
boundary. SDF recommends the simplest competitive candidate. Equivalent
placements use deterministic tie-breakers: concentrate disks around fewer
owners, prefer earlier hierarchy contributors, follow configured model-family
order, and then prefer fewer numerical parameters.

The recommendation and the published default are separate. An administrator
can override the default without changing or deleting retained candidates.

## Retained candidates and selection

Each completed candidate is committed before the search proceeds:

```text
package/
  fits/
    selection.json
    <readable-model-name>--<identity-hash>/
      result.json
      posterior-<content-hash>.npz
      corner.png
      distributions.png
    .work/
      ... sampler-native resumable output ...
```

`result.json` contains identities, component ownership, evidence, fit state,
resolved configuration, summaries, observation scopes, and provenance.
`posterior-*.npz` contains capped posterior samples, labelled derived samples,
observation/model arrays, component predictions, and model spectra.

`selection.json` independently records:

- the automatic recommendation;
- the public default;
- which retained candidates are shown as tabs; and
- whether the default came from automatic selection or an administrator.

Changing the selection never invokes the sampler. Re-fitting the same exact
candidate replaces its artifact rather than creating dated result history.
Other retained candidates remain available for comparison.

See [Using SDF fit outputs](fit-outputs.md) for the numerical schemas and code
examples.

## Freshness

SDF keeps separate fingerprints for input data and material fit state.

A candidate becomes non-current when, for example:

- rawphot values or flags change;
- observation topology changes;
- a spectrum is added, removed, or changed; or
- numerical model grids, priors, or material sampler settings change.

Moving the package, changing shown/default selection, changing an alias-only
identity file, or changing presentation code does not by itself make the
numerical fit stale. The SDF revision remains recorded as provenance.

Old candidates can be inspected with an **input changed** warning, but a stale
candidate should not remain the public default after a replacement fit is
available.

## Result index and samples

Package files are the numerical source of truth. `sdf-results.sqlite` is a
disposable projection that makes catalog and sample queries fast. It contains
packages, identities, candidates, ownership, evidence, freshness, selection,
and scalar fitted/derived quantities. It can always be rebuilt from retained
artifacts.

Sample membership remains an SDB concept. SDF projects visible samples through
the current default candidate for each system. Physical contributors remain
separate rows, disks stay under their owner, and composite observation rows
retain their contributor list.

## Command-line reference

Routine fitting decisions are intended to be made through the administrative
website. Package-level commands remain available for automation and recovery:

```shell
sdf-fits list PACKAGE
sdf-fits render PACKAGE
sdf-fits select PACKAGE --default FIT --display FIT OTHER_FIT
sdf-fits fit PACKAGE --stellar 'TARGET=MODEL' [--disk 'TARGET=MODEL']
sdf-fits refit PACKAGE FIT
```

Catalog-wide serving, workers, indexing, and maintenance are documented in
[Server and fitting operations](server-operations.md). `COMMAND --help` is the
authoritative option reference.

## Current compatibility boundary

Package fitting supports resolved and composite photometry, colours, upper
limits, multiple physical contributors, mixed stellar model families, and
segmented spectra. It produces shared public/admin SED pages and hierarchy-
aware sample products from compact retained candidates.

Ordinary unpackaged rawphot files still use the established single-target
`Result` workflow. Legacy MySQL result writing is intentionally unavailable
for joint packages; the retained artifacts and rebuildable result index replace
that single-target-shaped output. SDF does not currently write photometry or
association decisions back to SDB.
