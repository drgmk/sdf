# Using the `sdf` website

This guide describes the public website and the additional tools available to
ordinary signed-in users. It does not cover fit configuration or server
operation; those are documented in the [fit administration](admin-guide.md)
and [server operations](server-operations.md) guides.

## Searching for a system

The search box appears in the common page header. Search by a target name,
component name, or known catalogue identifier. Search ignores case and is
lenient about spaces, so `HD 38`, `hd38`, and similar compact catalogue forms
find the same indexed aliases.

Search results represent fitted systems rather than isolated catalogue rows.
A multiple system therefore has one SED page containing its physical
components and any composite observations.

Individual system pages are public. Signing in is not required to inspect a
fit or download its products.

## Reading an SED page

The page header gives the system name and a compact description of the default
model. In the shorthand used throughout the site:

- `★` is a stellar model component;
- `○` is a disk model component; and
- labels such as `A`, `B`, or `AB` identify hierarchy components or composite
  observations.

For example, `A ★○ · B ★` means a star and disk owned by A plus a star owned by
B. Model ownership matters: a disk around A is physically different from the
same numerical disk model placed around B.

The system panel at the left shows the hierarchy, component names, spectral
types and distances where available, links to the exported input data, and
external services such as SIMBAD. Composite observation scopes can appear
alongside the physical contributors they contain.

### Fit tabs and results

The first tab is the current public default. Other tabs are retained
alternatives that an administrator has chosen to show. Switching tabs changes
the model and result summary but not the underlying observed data.

Each fit contains:

- a log-log SED and residual panel;
- physical results grouped by the component that owns them;
- log evidence, chi-squared, and the number of residual elements;
- downloadable metadata and posterior arrays; and
- fitted-parameter and derived-quantity distribution plots when available.

The compact fit-result links are:

- **metadata (JSON)** for model identity, summaries, evidence, and provenance;
- **posterior data (NPZ)** for labelled fitted and derived samples;
- **parameters (PNG)** for fitted-parameter distributions; and
- **derived (PNG)** for physical derived-quantity distributions.

See [Using `sdf` fit outputs](fit-outputs.md) for Python examples.

### Plot symbols

The in-plot legend distinguishes data, upper limits, excluded data, and model
photometry. Used measurements are solid points, ignored measurements are
greyed, upper limits use downward triangles, and model photometry uses open
rings. Stellar, disk, total-model, and observed-spectrum lines are labelled
separately when present.

The residual panel uses the same data-state symbols. Residuals are in units of
the reported measurement uncertainty. Points retained for context but excluded
from fitting are still shown.

The toolbar below the plots controls both panels. Panning or zooming one panel
updates the shared wavelength range.

## Samples

The **Samples** link lists every sample visible to the current visitor:

- public samples are visible without signing in;
- private samples are visible to their owner;
- shared samples are visible to their collaborators; and
- administrators can see all samples.

A sample page has the same navigation on each view:

- **Properties** is a filterable, sortable table of stellar and disk results;
- **H–R and disks** contains population-level diagnostic plots; and
- **Flux and size** contains observing-planning plots.

Clicking a target or plot point opens its system SED page. A line connecting
plot symbols indicates components from the same hierarchy. Context components
needed to understand a system can appear even when only one member belongs
directly to the sample; these rows are marked **system context**.

The right side of the sample navigation provides downloads:

- **VOTable** contains the sample result projection; and
- **Photometry** is a long-form CSV retaining observation contributors, model
  ownership, flags, observed flux, and model intervals.

An **Awaiting a current fit** section means that membership is known but the
system has no current default result yet.

Sample owners and collaborators can use **Edit sample** on the Properties page
to add identifiers or remove current members. Additions are checked locally
first, then existing targets are reused and genuinely new targets are imported
through SIMBAD. Removal changes only sample membership: the target, its fit,
and memberships in other samples are retained. The Sharing tab accepts the
exact username or account email of an existing user; collaborators can view
and edit the sample and can add or remove other collaborators. Administrators
can edit any sample; changing a public sample is restricted to administrators.

## Signing in

Signing in adds private/shared samples and enables imports. It does not change
the public content of individual SED pages. An account must exist both in the
site's password system and in the `sdb` user database.

Ordinary accounts can:

- see public, owned, and shared samples in one list;
- import a missing SIMBAD target from an unsuccessful search; and
- create a private sample from a target list.

Fit configuration and `sdb` data review require an administrator account.

## Importing a missing target

If a search finds no fitted system and imports are enabled, a signed-in user
can request the exact search text as a SIMBAD identifier. The resulting status
page follows the job through `sdb` acquisition, provider collection, export, and
the first `sdf` fit. When complete it links to the new SED page.

Imports may take several minutes because remote catalogues and MultiNest are
not part of the web request itself. Provider warnings can be reported even
when a usable target and fit were created.

## Creating a sample

Choose **Create sample** from the sample list. Supply a sample name and either:

- identifiers pasted one per line; or
- a text, CSV, or TSV file containing up to 1000 identifiers.

CSV/TSV files need a column named `target`, `simbad_name`, `identifier`,
`name`, or `sdbid`. Blank lines and lines beginning with `#` are ignored.

The first step checks the list locally and reports:

- targets already present in `sdb`;
- new identifiers that require SIMBAD; and
- ambiguous local aliases that must be resolved before import.

Existing targets are reused rather than re-imported. After a successful check,
submit the sample and follow its status page. New samples are private by
default and visible only to their creator and administrators. Only an
administrator can make a sample public; explicit sharing is managed by an
administrator through `sdb`.
