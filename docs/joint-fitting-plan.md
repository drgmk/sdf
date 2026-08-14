# Joint fitting and result architecture plan

## Status and purpose

This document is the working design for extending SDF from fitting one
raw-photometry file at a time to fitting a package of resolved and unresolved
observations together.  It also proposes the limited restructuring needed to
make that extension maintainable.

The immediate compatibility requirement is that an ordinary single rawphot
file continues to fit as it does now.  Historical `Result` pickles and the
current results-database schema are not compatibility constraints.  Pickles
may be treated as disposable caches, and the database may be replaced by a
query index over the current useful fit for each target or system.

The operational goal is not to preserve successive historical versions of a
fit.  Most targets will fit automatically; a minority will need corrected
photometry or an explicit model choice.  For the current target data, however,
it is useful to retain every completed model fit compactly: an operator can
compare it with later alternatives or add it back to the public SED tabs
without rerunning MultiNest.  SDF should distinguish this candidate library
from the smaller display selection that downstream users normally see.  Each
candidate should record the SDF code state, model grids, fit configuration,
and input state so that it can be understood and its analysis regenerated.

Implementation is underway on `codex/joint-fitting`. The branch now has the
dataset/candidate/problem boundaries, local MultiNest callbacks, joint-package
discovery and projections, compact candidate artifacts, SIMBAD-keyed portable
policy, deterministic candidate search, orchestration, hierarchy-aware
analysis, joint plotting, fit selection commands, and a rebuildable result
index. Single-file fitting still follows the legacy facade. The operator slice
now has a shared public/operator SED page, persistent policy controls, retained
fit selection, and a fit-job queue/worker boundary. See `joint-fitting.md` for
the runnable development interface and current boundary.

The SDB export contract remains deliberately small:

```yaml
version: 1
observations:
  HD 113766:
    - HD 113766A
    - HD 113766B
```

In actual exports these names currently appear as SDB IDs.  SDF should resolve
the files to target identities on load and should use canonical SIMBAD names
for persistent user configuration and published identity.  The joint-fit file
describes observation topology only; it should not accumulate fitting policy.

## Design principles

1. **Single-file fitting is the one-target form of the same data model.**  It
   should not remain a separate implementation beside a joint pipeline.
2. **Identity, observation scope, and model ownership are different things.**
   A composite observation can see several physical target models without
   itself owning a physical model.
3. **Automatic fitting should remain conservative and evidence-led.**  User
   policy constrains which models may be tried; it should not normally choose
   a more complicated model without evidence.
4. **Do not infer disk policy from stellar class.**  Any type of stellar model
   can host a disk.
5. **Persistent operator choices use external identities, not generated SDB
   IDs.**  The primary key is the canonical SIMBAD main identifier, with a
   SIMBAD object identifier retained as an additional anchor when available.
6. **Retained fits and displayed fits are different sets.**  Compact results
   from completed candidates are reusable working products.  An explicit
   selection chooses the default fit and the smaller set shown to downstream
   users.  The database provides fast queries plus job/policy state, while the
   files contain the numerical results.
7. **Sampler-specific details stay behind a small interface.**  MultiNest is
   the initial evidence-capable sampler; emcee and other samplers should be
   possible without remodelling the fitting problem.

## Concepts and objects

The names below are provisional, but the responsibilities should remain
separate even if modules or classes are combined during implementation.

### `TargetIdentity`

Represents one physical model contributor.

It contains:

- a stable configuration key such as `simbad:HD 113766A`;
- the exact canonical SIMBAD `main_id` used for display;
- a SIMBAD numeric object identifier when available;
- known aliases needed to reconnect renamed or re-exported targets;
- the current SDB ID as a source locator, not persistent identity;
- coordinates, spectral type, parallax, and their provenance;
- the rawphot file from which this target's resolved observations and metadata
  were loaded.

SIMBAD names require a normalized lookup representation for matching, while
the exact canonical spelling should be retained for display.  A target without
a SIMBAD identity needs a fallback key; anything that has photometry to fit
should exist in some catalogue from which we can generate a name, for example
many 2MASS objects are not in simbad.
Its generated SDB ID should not silently become a long-term policy key.

### `Observation`

Represents one observed SED scope, not necessarily one physical object.

It contains:

- its file and file/header identifier;
- its display identity, which may be a physical or composite SIMBAD object;
- photometry and optional spectra;
- the ordered set of physical contributor keys visible to the observation;
- observation-specific spectrum-normalization parameters;
- provenance needed to label plotted points and result rows.

For an A/B/AB package the observations are:

| Observation | Contributors seen by the likelihood |
|---|---|
| A | A |
| B | B |
| AB | A + B |

An observation omitted from `joint-fit.yml` is a self-observation.  The YAML
only records non-default contributor mappings.

### `FitDataset`

The validated, sampler-independent input to fitting.

It contains:

- all `TargetIdentity` objects;
- all `Observation` objects;
- the observation-to-target topology;
- concatenated observed arrays and their observation boundaries;
- per-target metadata used by derived analysis;
- a compact input-state marker covering the rawphot files, topology, spectra
  dependencies, and resolved fit policy.

Responsibilities:

- discover a package from a rawphot file or directory;
- read and validate `joint-fit.yml`;
- ensure every contributor is a physical target with an input file;
- ensure observation and contributor identifiers resolve unambiguously;
- flatten observations without losing their scope labels;
- make the one-file/no-YAML case a normal self-mapped dataset;
- provide an input-state marker for freshness checks.

It does not select model families or perform fitting.

### `ModelFamilyRegistry`

Describes the available numerical model grids and their capabilities.

For each model family it records:

- model name and kind (`stellar`, `disk`, or another future kind);
- parameter names and ranges;
- whether spectral and photometric grids exist;
- automatic-selection rules, such as supported spectral classes;
- whether it may be used as an alternative to another family;
- any constraints required to combine it with observations.

This replaces string-substring classification as the eventual source of truth.
The existing configuration-derived lists can populate the registry during the
transition.

### `TargetModelPolicy`

The resolved fitting policy for one physical target.

Example portable form:

```yaml
targets:
  HD 48915 A:
    stellar_models: [phoenix_sol]
    max_disks: 2
  HD 48915 B:
    stellar_models: [koester_wd]
    max_disks: 2
```

The target keys are SIMBAD main identifiers.  `max_disks` limits the search; it
does not require that number of disks.  A value of zero is an explicit user
constraint, not a default consequence of the chosen stellar model.

Policy fields should initially include:

- allowed stellar model family or families;
- allowed disk model family or families;
- maximum disk components;
- automatic, evidence-guided, exhaustive, or exact-candidate search mode;
- optional reddening policy;
- origin (`automatic`, central policy, package file, CLI, or API);
- actor, reason, and modification time for operator changes.

Policy precedence should be:

1. exact API or command-line request;
2. package-local `sdf-fit.yml`, for portable/reproducible work;
3. central operator policy store;
4. legacy single-target model configuration while it is supported;
5. spectral-type-based automatic selection.

The fully resolved policy and its origin are always copied into the result.

### `ModelComponentSpec`

Describes one physical model component without loading its numerical grid.

It contains:

- a stable component ID within the candidate;
- owner `TargetIdentity` key;
- role (`stellar` or `disk`);
- exact model family name;
- ordinal for repeated components;
- parameter labels qualified by owner and component;
- the observations in which it contributes, derived from topology.

Ownership prevents two identical stellar grids fitted to A and B from being
mistaken for exchangeable components.  It also tells derived analysis which
parallax and stellar luminosity apply to a disk.

### `ModelCandidate`

A complete model hypothesis whose evidence can be evaluated.

It contains an ordered tuple of `ModelComponentSpec` objects and a stable,
human-readable candidate ID.  It can be compiled against a `FitDataset`, but
does not contain posterior results.

Examples for a binary include:

- one stellar component for A and one for B;
- the same plus a disk owned by A;
- the same plus a disk owned by B;
- two disks around A, two around B, or one around each.

Stellar alternatives are also candidate dimensions.  A Sirius-like system can
therefore use a normal stellar model for A and `koester_wd` for B, while a cool
companion can use `phoenix_cool` independently of the primary.

### `CandidateSearch`

Generates and evaluates a graph of `ModelCandidate` objects.

The current binary tree is evidence-led but assumes one stellar contributor
and a linear sequence of disk additions.  Joint fitting requires a graph:

1. Generate minimal roots with one allowed stellar model per contributor.
2. Fit each unique root when stellar-family alternatives exist.
3. From every evidence-supported root, generate all unique candidates formed
   by adding one permitted disk component to one contributor.
4. Fit those children and continue only branches whose evidence exceeds their
   parent by the configured threshold.
5. Stop at target and package component limits.

Several children at the same depth must be evaluated before pruning.  A
greedy choice of only the best first disk would make the result depend on
generation order.  The search process needs fitted candidates and parent
evidence while it is running.  Every completed candidate should be committed
to the package's compact candidate library before the search continues.  This
both protects slow MultiNest work and makes all previously fitted alternatives
available to an operator.  Search-frontier state and reasons that an unfit
branch was not expanded do not need permanent storage.

Search modes:

- **evidence-guided**: the default graph expansion above;
- **exhaustive**: fit every candidate within explicit limits, useful for
  validation and small unusual systems;
- **exact**: fit only explicitly supplied complete candidates;
- **posterior-only**: run a sampler for an already selected candidate without
  requiring evidence.

Evidence-guided pruning can theoretically miss a model whose two-component
extension is supported even though neither one-component parent is.  This is
also a limitation of the current tree.  Exhaustive mode provides the check;
we can later consider a cautious look-ahead policy if real systems require it.

### `FitProblem`

A compiled, pure description of one dataset/candidate likelihood.

It contains:

- flattened observations and uncertainties;
- loaded/interpolated numerical models;
- parameter definitions, bounds, and prior transforms;
- observation boundaries;
- the observation-by-component projection matrix;
- spectrum-normalization parameters;
- methods for model prediction, residuals, chi-squared, and log likelihood.

The projection matrix is the core joint-fitting addition.  A component is
evaluated against all required filters/spectra but contributes only to
observation scopes that include its owner.  The default single-target
projection is all ones and must reproduce the current calculation.

`FitProblem` must have no plotting, filesystem, database, or sampler logic.  It
should be easy to test with synthetic A, B, and AB fluxes.

### `Sampler` protocol and implementations

The sampler interface is approximately:

```python
class Sampler(Protocol):
    capabilities: SamplerCapabilities

    def run(
        self,
        problem: FitProblem,
        output_dir: Path,
        options: SamplerOptions,
    ) -> PosteriorRun: ...
```

`SamplerCapabilities` states whether a sampler supplies evidence, weighted
samples, resumable runs, and convergence diagnostics.  Candidate selection
requires an evidence-capable sampler; posterior-only refinement does not.

Initial implementations:

- `MultiNestSampler`: evidence and posterior samples, used for automatic model
  search;
- `EmceeSampler`: possible later posterior refinement for an exact selected
  candidate, but not itself a replacement for evidence-based selection.

The MultiNest adapter should use local closures that capture `FitProblem` and
its prior.  PyMultiNest retains these callbacks for the duration of `run`, so
the current module globals are unnecessary.  The wrapper's `context` value is
not passed through to the user likelihood.  Process-based parallelism remains
the safe default because the underlying library should not be assumed to be
thread-safe.

### `PosteriorRun`

Sampler-neutral output from fitting one candidate:

- samples and parameter labels;
- weights or equal-weight status;
- log likelihoods;
- evidence and its uncertainty when available;
- sampler name, version, options, and diagnostics;
- run state and timestamps;
- paths to sampler-native checkpoint/output files.

This object contains no stellar/disk derived quantities.

### `FitAnalysis`

Converts a `PosteriorRun` plus its `FitProblem` into scientific result values.

It computes:

- posterior summaries of fitted parameters;
- predictions and residuals for each observation scope;
- per-component spectra and photometry;
- per-target stellar quantities using that target's parallax;
- per-target disk quantities using the owning target's stellar luminosity;
- combined system spectra and photometry;
- uncertainty distributions needed by plotting and sample tables.

The existing analysis code should be migrated here in stages.  It currently
assumes one parallax and often totals all stellar components before deriving
disk quantities, which is not sufficient for hierarchy-aware results.

### `CandidateResult`, `CandidateLibrary`, and `FitSelection`

`CandidateResult` combines one candidate, its posterior, and its analyzed
products.  Once fitting succeeds it is saved as a compact reusable artifact.
It contains:

- package and input identity;
- the selected model components and their ownership;
- fitted and derived results;
- observations, predictions, and residuals;
- warnings, degeneracies, and unresolved ownership;
- the resolved model policy;
- fit time plus the SDF revision, model-grid names/revisions, sampler, and
  material configuration state;
- enough posterior samples to rerun `FitAnalysis` and regenerate plots without
  rerunning the sampler.

`CandidateLibrary` discovers and indexes every completed candidate artifact in
one target/package.  It determines whether a candidate is current for the
present data and fitting inputs, can find a previously computed exact model,
and supplies candidates to the operator UI.  It does not decide which fits are
scientifically preferred.  Evidence rankings are only meaningful among
candidates with compatible data state, priors, and evidence calculation; the
library should not rank stale or materially incompatible candidates together.

`FitSelection` is a small mutable record containing:

- the default candidate;
- the ordered candidates shown as SED tabs;
- optional short operator labels or notes.

A package can therefore retain several or many computed candidates but expose
only a useful subset.  Changing the default or displayed tabs is immediate and
does not run MultiNest.  The automatic search can propose a selection, while
an operator can revise it.  Sample extraction normally uses the default fit.

Near-degenerate models that differ only in ownership of an unresolved disk
should be reported as such rather than arbitrarily presented as a secure host
assignment.

### `ResultStore`

Owns filesystem serialization and freshness checks.  It has no scientific
analysis logic.

Responsibilities:

- commit each completed candidate atomically;
- store a compact, normally equal-weight and capped posterior sample;
- write the display/default selection atomically;
- retain sampler-native output only while useful for resuming the current run;
- compare the current data and fit state with each stored candidate;
- replace an exact candidate when that same model is deliberately refit,
  without creating a chain of historical versions;
- provide loaders used by plotting without requiring Python pickles.

Historical `Result` pickles are not required to load.  A converter may be
written later if a specific old result is valuable.  New pickles, if retained
at all, are disposable acceleration caches rather than candidate artifacts.
There is still no chain of historical versions: refitting the same candidate
for the same/current inputs replaces its compact artifact.

### `FitCoordinator`

Orchestrates a run without absorbing the responsibilities above:

1. load a `FitDataset`;
2. resolve central/package/CLI policy;
3. ask `CandidateSearch` for candidates;
4. build a `FitProblem` and invoke a compatible sampler for each candidate;
5. analyze completed posteriors;
6. commit every completed candidate to the candidate library;
7. update or preserve the display/default selection;
8. update the result index;
9. generate or enqueue plots and pages.

The current `fit_results()` API can become a compatibility wrapper around the
coordinator for a single rawphot file.

## Model-selection UI and policy persistence

### Operator interaction

The operator SED page should show one row per physical contributor with:

- canonical SIMBAD name and aliases;
- current spectral type and automatically inferred stellar model;
- a stellar-model selector (`Auto` plus available model families);
- allowed disk family and maximum disk components, defaulting to `Auto`;
- the origin of the active policy;
- reset-to-automatic and refit actions.

System-level advanced controls can expose search mode, evidence threshold, and
exact candidates.  Common changes should not require editing YAML or locating
a package directory.

The operator page should also expose an **Available fits** table populated by
the candidate library.  For each completed candidate it should show the owned
model components, evidence relative to the best current candidate, fit time,
and whether it is current for the present inputs.  Actions should include:

- preview this candidate without changing the public page;
- add it to or remove it from the SED tabs;
- make it the default fit;
- regenerate its analysis/plots from the saved posterior;
- rerun MultiNest when the candidate is stale or explicitly requested.

If an operator requests a model combination that is already current in the
library, the UI should reveal or select it immediately.  Only a genuinely new
or stale combination should enqueue a slow fit.  The public SED page shows the
ordered `FitSelection` tabs, with the default first; it does not need to expose
the full candidate library.

Changing a policy should show a concise preview, record actor/reason, save the
central policy, and enqueue a refit.  It should not block an HTTP request while
MultiNest runs.  Fit job state belongs in the central state database. Interaction
should be low-friction; auto filled actor/reason fields, minimal clicks.

### Public versus operator pages

Current SED pages are static public products and remain read-only. The joint
result renderer uses one template and one set of component-aware projections
for both audiences. Operator mode adds a right-hand controls drawer to that
same page; it is not a separately designed application. The static operator
page is a useful rendering/development shell, but persistence and job launch
must be supplied by an authenticated or local backend endpoint.

For multiple systems, the page distinguishes physical contributors from
observed scopes. Each physical contributor has its own identity, input, and
external links. A resolved observation uses that contributor's plot colour;
an unresolved/composite observation gets a neutral system colour and is
labelled by the contributors it sees. Candidate plots show thin owner-coloured
component spectra and thicker scope-level total spectra, with photometry and
residuals grouped by observed scope. This keeps the plot scientifically honest
without inventing per-component photometry for unresolved measurements.

### Central policy store

Operator policy is centralized rather than written into one of many target
folders. A small `sdf-operator.sqlite` store owns policy and the fit-job queue.
It is separate from the rebuildable `sdf-results.sqlite`, since operator state
must not disappear during an index rebuild. Policy is normally target-global
and keyed by canonical SIMBAD identity; a package-scoped exception is used
when the same object genuinely needs different treatment in a particular
analysis.

SDF should snapshot the resolved policy into every fit result.  A portable
package-local `sdf-fit.yml` remains useful for fits made outside the main SDB
installation and for reproducible transfer to another machine.

## Result products

### Directory layout

The exact names can be adjusted, but each exported package remains
self-contained:

```text
<package>/
  <target-a>-rawphot.txt
  <target-b>-rawphot.txt
  <composite>-rawphot.txt
  joint-fit.yml
  index.html
  fits/
    index.json
    selection.json
    A-phoenix_sol+B-koester_wd/
      result.json
      posterior.npz
      corner.png
      distributions.png
    A-phoenix_sol+B-koester_wd+A-modbb_disk_r/
      result.json
      posterior.npz
      corner.png
      distributions.png
    .work/
      ... interrupted/current sampler runs ...
```

Candidate directories should have readable names based on their owned model
components.  A single target can retain familiar names such as
`phoenix_sol+modbb_disk_r`.  A joint candidate qualifies models with short
package-local component labels, for example
`A-phoenix_sol+B-koester_wd+A-modbb_disk_r`.  Repeated components gain an
ordinal.  The full SIMBAD owner identities and model specification live in
`result.json`; directory names are labels, not identity keys.  Very long names
may be shortened with a deterministic suffix.

`fits/index.json` is the inexpensive catalogue of available candidate
directories and their evidence/freshness.  It can be regenerated by scanning
the candidate directories.  `fits/selection.json` records the default and
ordered SED tabs.  The central result index mirrors both for queries and the
operator UI.

MultiNest writes into `fits/.work/<candidate>/`.  On successful completion SDF
first commits the compact candidate artifact, then may remove the large native
sampler files.  Interrupted output can remain resumable.  Keeping completed
native output may be a diagnostic option, but is not required to regenerate
the analysis or plots.

### Candidate `result.json` and `posterior.npz`

Each candidate's compact result records:

- candidate model components and physical ownership;
- evidence and fitted-parameter summaries;
- package identity plus separate data-state and fit-state markers;
- physical targets, SIMBAD identities, aliases, and current SDB locators;
- observations and contributor mappings;
- fitted parameters plus per-target and system derived values;
- observation-level predictions, residuals, exclusions, and provenance;
- references to the posterior and plots;
- warnings and unresolved degeneracies;
- fit time, SDF revision, sampler, model-grid identity, and material
  configuration state.

`posterior.npz` retains the equal-weight posterior used by analysis, normally
capped at the existing `n_samples_max`, along with parameter labels and any
weights/log likelihoods that are still needed.  Existing examples suggest
that this is compact: equal-weight text files of roughly 287--461 kB compressed
to approximately 62--109 kB as NPZ.  Dense per-sample predictions for every
filter need not be retained because `FitAnalysis` can regenerate them from the
posterior and model grids.  A regenerable analysis cache can be added only if
that calculation proves inconveniently slow.

There is no requirement to retain old result formats or build a general schema
migration system.  If the stored structure later changes incompatibly, it is
normally simpler to refit or regenerate the result.  A small format marker can
be added when an actual second format exists.  Large posterior samples and
spectra should not be expanded into JSON; they belong in the accompanying
array product.

The data-state marker answers whether the rawphot, spectra, or joint topology
has changed.  The fit-state marker covers the numerical model grids, priors,
and sampler settings that affect this candidate.  Merely changing which tabs
are displayed, or changing search policy without changing this exact
candidate, does not make its posterior stale.  A stale candidate can remain
visible to operators until rerun, but should not be used as the public default.
Refitting the same model combination replaces that candidate directory rather
than creating dated copies.

The recorded SDF revision is primarily provenance.  A new commit should not
automatically make every fit stale.  If a code change materially alters the
likelihood or parameter interpretation, the operator can explicitly require a
refit (and a lightweight fit-compatibility marker can be introduced if this
becomes common).  Changes confined to analysis or plotting should regenerate
those products from `posterior.npz` without MultiNest.

### Plot loading

Plotting should load selected candidate results and numerical model grids, not
unpickle a large live `Result` object.  A small view/adapter object may provide
convenient arrays to existing plotting functions while they are migrated.  The
public HTML generator reads `selection.json` and creates one tab per selected
candidate, preserving the spirit of the current SDF page.  Operator previews
can load any candidate in the library without changing that selection.

## Query index and sample generation

Scanning ten thousand package directories and all their candidate subfolders
for every sample page is avoidable.  Use a central SQLite result index as a
materialized view of candidate libraries and display selections.

Likely tables include:

- `fit_packages`: package identity, input state, default fit, status,
  and update time;
- `fit_targets`: canonical SIMBAD identity, aliases, current SDB ID, and
  package membership;
- `fit_candidates`: every retained candidate's model signature, evidence,
  freshness, location, and fit time;
- `fit_selection`: default and ordered displayed candidates;
- `fit_components`: target ownership, model family, role, and ordinal;
- `fit_parameters`: selected scalar fitted/derived values and units;
- `fit_jobs`: queued/running/completed/failed refits;
- `fit_policy`: centralized operator policy and audit fields.

The index should be updated transactionally whenever a candidate or selection
is committed.  Refitting an exact candidate replaces that candidate's rows;
other model combinations remain available.  The index should also be
rebuildable from the candidate directories and selections, for example with:

```text
sdf-index rebuild <export-root>
sdf-index verify <export-root>
```

`sdf-sample` normally queries the default candidate for target/system
quantities and follows its result path for richer material.  It may explicitly
request another retained candidate when needed.  Sample membership can remain
owned by SDB;
the export manifest or identity mapping connects sample members to result
packages.  If the index is absent, SDF can build it once rather than repeatedly
trawling the tree.

This preserves the performance reason for the old results database without
requiring dated result history.  Rebuilding means indexing the retained model
candidates and current display selection, not reconstructing fitting attempts
that failed or were overwritten for the same candidate.

## Proposed source structure

An incremental target layout is:

```text
sdf/
  dataset.py
  model_policy.py
  model_selection.py
  likelihood.py
  samplers/
    __init__.py
    base.py
    multinest.py
    emcee.py             # later
  analysis.py
  result_io.py
  candidate_library.py
  result_index.py
  coordinator.py
  result.py              # compatibility facade during migration
```

This is not a requirement to create every module immediately.  Modules should
appear when code and tests establish a real boundary.

## Implementation sequence

### Phase 0: branch, baseline, and direct defects

- Work on `codex/joint-fitting`; the installed non-editable SDF remains
  unchanged while the source branch is developed.
- Preserve the existing `.gitignore` change and generated example files.
- Record current single-file CLI and result behaviour with characterization
  tests.
- Fix the current Python 3.12/NumPy PHOENIX-reader test failure.
- Fix defects that obstruct multiple observations: `concat_obs()` appending
  filters from the wrong array, the NumPy `nused` identity comparison, the
  zero-flux `sum` typo, and exact rawphot suffix removal.

### Phase 1: pure single-target fitting core

- Add single-target `TargetIdentity`, `Observation`, and `FitDataset`.
- Add `ModelComponentSpec`, `ModelCandidate`, and `FitProblem`.
- Move prior/likelihood state out of module globals.
- Add `MultiNestSampler` using closures.
- Run existing single-file fitting through this core and compare likelihoods,
  evidence, parameters, paths, and plots with the baseline.

### Phase 2: compact candidate library

- Define per-candidate `result.json` and compact posterior arrays.
- Make pickles optional/disposable and add product-based plot loading.
- Introduce `ResultStore`, `CandidateLibrary`, and separate data/fit freshness
  state.
- Add `FitSelection` for the default and ordered displayed tabs.
- Commit every completed candidate before removing its native MultiNest files.
- Retain a `result.Result` facade only where it keeps current callers simple.

### Phase 3: policy and model search

- Build the model-family registry and spectral-type inference.
- Add portable `sdf-fit.yml` parsing.
- Add a central policy-store interface keyed by SIMBAD identity.
- Implement evidence-guided candidate graphs plus exhaustive and exact modes.
- Reuse already-current candidates during search and retain every newly
  completed candidate in the library.

### Phase 4: joint datasets and likelihoods

- Discover and validate `joint-fit.yml` packages.
- Deduplicate package members collected by `--file`, `--dir`, and sample runs.
- Add package-level locks.
- Compile observation/component projections.
- Test synthetic A, B, and AB photometry, colours, exclusions, upper limits,
  and spectra.
- Verify the no-YAML projection remains numerically identical to the
  single-target baseline.

### Phase 5: hierarchy-aware analysis and presentation

- Derive stellar quantities with each target's own parallax.
- Associate disk quantities with their owning target and flag unresolved host
  degeneracy.
- Produce observation-, target-, and system-level predictions.
- Update plots, HTML, JSON, and thumbnails to label scopes and model owners.
- Run an end-to-end HD 113766 rehearsal and at least one normal-star/white-
  dwarf system.

### Phase 6: operator workflow and indexing

- Add the SQLite policy/job/result index.
- Add the Available fits UI, display/default controls, and asynchronous refit
  requests.
- Keep public SED pages read-only.
- Change `sdf-sample` to query the default candidate through the result index.
- Provide index rebuild and verification commands.

The result index, rebuild/verification commands, fit library, selection
commands, shared SED page, persistent policy state, functional drawer controls,
and asynchronous job boundary are implemented. Operator policy/job state uses
a separate `sdf-operator.sqlite`: it is durable operational state, whereas
`sdf-results.sqlite` remains a disposable materialized view.

The first hierarchy-aware `sdf-sample --index` path is also implemented. The
index now records composite observation scopes and their ordered contributors.
Sample IDs may therefore resolve either physical targets or unresolved
observations. The new HTML table and VOTable show direct members plus physical
system context, stellar results, and explicitly owned disks. The new Bokeh
plots provide an HR diagram and disk fractional luminosity versus temperature.
The legacy path remains the fallback when `--index` is omitted.

Completing phase 6 now means deciding whether the legacy raw-photometry CSV is
still useful, and adding the compact per-band model predictions needed to
replace flux/size plots without rebuilding old `Result` objects. Operational
worker hardening can follow if concurrent use requires it.

### Phase 7: optional samplers and performance

- Add `EmceeSampler` for exact-candidate posterior refinement if useful.
- Measure candidate-search and repeated model-interpolation costs.
- Cache compiled/interpolated models at safe dataset/candidate boundaries.
- Consider process-level parallel candidate fitting where MultiNest and file
  locking permit it.

## Validation gates

The implementation should not advance past each gate without:

1. single-file likelihood and result regression tests;
2. deterministic candidate and parameter naming;
3. candidates become stale when data or material fit inputs change, while a
   display-only or irrelevant policy change does not force a refit;
4. no fitting outcome dependent on rawphot discovery order;
5. no model ownership inferred from component ordering;
6. explicit units and identity provenance in candidate results;
7. the SDF revision, model-grid identity, sampler, and material configuration
   are recorded in each retained candidate;
8. ability to rebuild the query index from candidate directories and
   `FitSelection`;
9. a real joint-system rehearsal whose resolved and composite residuals can be
   inspected in the SED UI.

## Decisions still to validate with real fits

- The default evidence threshold and whether it should differ for alternative
  stellar families versus adding disk components.
- Whether evidence-guided search needs limited two-step look-ahead to catch
  jointly supported additions whose individual parents are weak.
- How much posterior and model-prediction detail is needed by plotting and
  downstream use, and therefore belongs in `result.json` versus the array
  product.
- Whether the existing cap of 1,000 equal-weight posterior samples is adequate
  for every reanalysis/plotting use, and whether any weights or likelihood
  values must also be retained.
- Whether completed native MultiNest files should be removed by default or
  retained for selected fits only; the compact candidate must be sufficient
  for ordinary regeneration either way.
- Whether the central policy/job/index store belongs operationally in SDF, in
  the main SDB SQLite database, or behind a small interface with either
  implementation.  The domain model should not depend on this deployment
  choice.
- How SIMBAD aliases and numeric object identifiers are exported so policies
  survive a changed canonical main ID as well as a future SDB ID generation.
