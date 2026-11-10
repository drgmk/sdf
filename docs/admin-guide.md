# Fit administration

This guide is for users with the `sdb` administrator role. It covers decisions
made through the web interface: reviewing inputs, configuring model searches,
queueing fits, and selecting the result shown publicly. Deployment and process
management are covered separately in
[Server and fitting operations](server-operations.md).

## Administrator pages

After signing in, the common header adds:

- **Admin**, the fitting dashboard;
- **Data review**, the `sdb` association and photometry-review application; and
- **Jobs**, the persistent target-import and fit queues; and
- **Activity**, catalogue inventory, website use, recent actions, and job
  health.

The normal and administrative SED pages share the same layout and results. The
administrative page adds a **Fit controls** drawer and a **Review data** link.
This keeps what an administrator evaluates visually identical to what a public
visitor sees.

## Site activity

The **Activity** link opens `/admin/activity`. It combines three existing
sources without creating a second domain-event store:

- the disposable result index for fitted-system inventory;
- the `sdb` and `sdf` audit/job databases for imports, sample edits, review
  actions, exports, fits, and their authenticated actors; and
- an optional Apache access log for aggregate page views, downloads, popular
  SEDs and samples, searches, errors, and request timing.

Choose 7, 30, or 90 days, or all recorded time. Access-log processing ignores
known crawlers, assets, redirects, failed page requests, and API polling. It
does not display or copy client addresses. Cloudflare Web Analytics remains
the better source for visits, referrers, geography, device information, and
browser-side performance.

Pass `--access-log PATH` to `sdf-fits serve` to enable the local website-use
sections. Use `--activity-timezone ZONE`, for example
`--activity-timezone Pacific/Auckland`, to group daily views and display event
times in that IANA timezone. The page states the active zone and abbreviation.
Inventory, changes, and job health remain available without an access log.

## Dashboard states

The dashboard lists every package currently exported by `sdb`, not just packages
that already have a result. It can be searched and filtered by state:

- **current**: the retained default fit matches the exported inputs;
- **unfit**: no completed candidate is available;
- **stale**: photometry, spectra, hierarchy, or another material input changed
  after fitting; and
- **needs selection**: retained results exist but no valid public selection is
  available.

The dashboard is paginated so it remains usable for catalogs containing many
thousands of systems. Opening an unfit package presents a single action to
queue its initial automatic fit. The authenticated account supplies the audit
actor; the pre-filled reason can be edited.

## Review data before changing models

Use **Review data** when the problem is an association, hierarchy, inclusion,
or provider-data decision. Those changes belong in `sdb`. `sdf` fit controls alter
model configuration and result presentation; they do not edit photometry or
write review decisions back to `sdb`.

After `sdb` exports changed inputs, the `sdf` watcher indexes the package and
queues new or stale exports according to the server configuration. The old fit
may remain visible with an **input changed** warning until the replacement
finishes.

## Model configuration

Open **Fit controls** on an administrative SED page. Configuration is shown
separately for every physical contributor because components in one hierarchy
can need different stellar models or disk limits.

The normal state is **default config**. It displays the choices `sdf` would use
from exported metadata and automatic rules without storing a user override.
Choose **Specify custom config** only when the default is unsuitable.

A custom contributor configuration contains:

- **Stellar model**: the allowed stellar atmosphere family;
- **Disk model**: the disk family available to the search, or not used;
- **Maximum disks**: an upper bound, not a request to fit that many disks; and
- **Search**: normally evidence guided, with exhaustive search available for
  exceptional cases.

Configuration is normally keyed by the physical SIMBAD identity so it follows
that target across exports and future `sdb` identifiers. The page indicates when
a package-specific configuration is being used.

Changing fields in the drawer does not fit immediately. Enter or accept the
reason and press **Queue fit** once for the complete system. All contributor
changes are saved together and one background job evaluates the joint package.
An unchanged candidate whose inputs and numerical configuration are current is
validated and reused; queueing is not a forced rerun. Use `sdf-fits refit` for
the uncommon diagnostic case that needs a completely fresh sampler run.

Use **Use default config** to remove the custom override and return to the
automatic choices. This is different from queueing an unchanged automatic fit:
it changes which configuration will be used for future fits.

## Evidence-guided model search

The default search evaluates complete complexity levels. For example, every
allowed one-disk placement in a binary is evaluated before `sdf` decides whether
to proceed to two disks. The configured log-evidence threshold is 5.

After fitting, candidates within 5 of the maximum evidence are treated as
competitive, with evidence uncertainty handled conservatively near the
boundary. `sdf` recommends the simplest competitive candidate. When equally
complex placements are competitive, deterministic preferences favor:

1. disks concentrated around fewer physical contributors;
2. the earlier hierarchy contributor, such as A before B;
3. configured model-family order; and
4. fewer numerical parameters.

These rules never prefer a simpler model across an evidence difference greater
than the threshold. They provide a consistent automatic recommendation; an
administrator can still choose a different retained candidate.

## Available fits

The lower part of the drawer lists every retained candidate for the current
inputs. Each entry shows its component shorthand, evidence, and state.

- **Preview** opens that candidate without changing the public page.
- **shown** controls whether it appears as a public fit tab.
- **Make default** makes it the first public tab and the result used by sample
  products.
- **Use recommended** returns the default to `sdf`'s current evidence-based
  recommendation.

The default candidate must always be shown. Previewing is harmless: it neither
changes selection nor queues MultiNest. Changing shown/default selection is
also immediate and does not refit.

**Recommended** is `sdf`'s current automatic choice. **Default** is the choice
actually published. They differ when an administrator has intentionally
overridden the recommendation.

## Jobs and failures

The **Jobs** page separates target imports from fits. Jobs persist across
server restarts and show their requester, reason, timestamps, state, and error
message. Fit states are queued, running, completed, or failed.

If a fit fails, first decide whether the failure is:

- an input-data problem, handled through `sdb` review and re-export;
- an inappropriate model configuration, corrected in Fit controls; or
- an operational failure, diagnosed from the server/worker logs.

Queueing a corrected system creates a new job. A package-level lock prevents
two workers fitting the same system simultaneously.

## Target and sample imports

Any signed-in user can request a missing target or create a private sample.
Administrators can also create public samples and inspect all import jobs.
Private-sample owners and collaborators can add or remove members and manage
collaborators from the sample Properties page. Administrators can edit any
sample, including public and legacy ownerless samples. Membership edits retain
the target and fitting products and publish a fresh sample export automatically;
sharing changes take effect immediately without an export.
Account roles, sample ownership, sharing, visibility, and deletion are managed
with the `sdb` user/sample commands described in
[Authentication and authorization](catalog-auth.md).
