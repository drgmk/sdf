# Plotting a multiple system

This is an extension of [the single-system plotting
recipe](fit-output-plots.md). Use it for a joint fit containing two or more
physical contributors, such as A, B, and unresolved AB observations.

No different file format is involved. The extra information is carried by
three sets of labels:

- each model component has an `owner_key`;
- each observation scope has one or more `contributors`; and
- the JSON `analysis.targets` maps owner keys to display names and physical
  summaries.

An owner is a physical target. A component is a stellar or disk model belonging
to that target. A scope describes one imported photometry file or observed
spectrum and says which targets contribute to it.

## Inspect the structure first

Load `result.json` and `sed-data.npz` as in the main plotting guide, then print
a compact description:

```python
targets = result["analysis"]["targets"]

for component in result["components"]:
    owner = component["owner_key"]
    print(
        targets[owner]["main_id"],
        component["role"],
        component["model_name"],
        component["component_id"],
    )

for scope in sed_metadata["data_scopes"]:
    names = [targets[key]["main_id"] for key in scope["contributors"]]
    print(scope["observation_id"], "sees", " + ".join(names))
```

For an unresolved AB observation, `contributors` contains both A and B. Do not
infer that from an observation filename or an `AB` suffix.

## Colour physical contributors consistently

The component arrays are parallel. Build an owner-to-colour mapping once and
use it for smooth spectra, data, posterior draws, and annotations:

```python
owner_keys = list(result["analysis"]["targets"])
palette = plt.get_cmap("tab10")
owner_colour = {key: palette(index) for index, key in enumerate(owner_keys)}

component_ids = sed["component_ids"].astype(str)
component_owners = sed["component_owner_keys"].astype(str)
component_roles = sed["component_roles"].astype(str)
component_spectra = sed["component_spectrum_median"]
model_wave = sed["spectrum_wavelength_micron"]

fig, ax = plt.subplots(figsize=(9, 5.5))
for component_id, owner, role, spectrum in zip(
    component_ids, component_owners, component_roles, component_spectra
):
    label = f"{targets[owner]['main_id']} {role}"
    ax.loglog(
        model_wave, spectrum,
        color=owner_colour[owner],
        ls="--" if role == "stellar" else ":",
        alpha=0.75, label=label,
    )
```

## Plot the model seen by each observation

A composite observation must be compared with the sum of only the components
owned by its contributors. A resolved A observation must not include B. The
following adds one smooth total curve for every distinct contributor set:

```python
seen_contributor_sets = set()

for scope in sed_metadata["data_scopes"]:
    contributors = tuple(scope["contributors"])
    contributor_set = frozenset(contributors)
    if contributor_set in seen_contributor_sets:
        continue
    seen_contributor_sets.add(contributor_set)

    visible = np.isin(component_owners, list(contributor_set))
    total = component_spectra[visible].sum(axis=0)
    label = " + ".join(targets[key]["main_id"] for key in contributors)
    colour = (
        owner_colour[contributors[0]] if len(contributors) == 1 else "black"
    )
    ax.loglog(model_wave, total, color=colour, lw=2.0,
              alpha=0.8, label=f"total seen by {label}")
```

This is the principal difference from the single-system example, where
`component_spectra.sum(axis=0)` is sufficient.

## Add resolved and composite data

Use the scope slice exactly as for a single system. Colour a resolved scope by
its one contributor and a composite scope in black. This abbreviated example
shows detections; carry over ignored points, upper limits, spectra, model
photometry, and residuals from the complete single-system recipe if needed.

```python
wave = sed["wavelength_micron"]
flux = sed["observed_flux"]
error = sed["observed_uncertainty"]
ignored = sed["observed_ignore"].astype(bool)
upper_limit = sed["observed_upper_limit"].astype(bool)

for scope in sed_metadata["data_scopes"]:
    section = slice(scope["element_start"], scope["element_stop"])
    contributors = scope["contributors"]
    colour = owner_colour[contributors[0]] if len(contributors) == 1 else "black"
    label = " + ".join(targets[key]["main_id"] for key in contributors)

    x = wave[section]
    y = flux[section]
    dy = error[section]
    omit = ignored[section]
    limit = upper_limit[section]
    valid = (
        np.isfinite(x) & (x > 0) & np.isfinite(y) & (y > 0) &
        ~omit & ~limit
    )

    if scope["kind"] == "ObsSpectrum":
        ax.plot(x[valid], y[valid], color=colour, lw=1.0, label=label)
    else:
        ax.errorbar(x[valid], y[valid], yerr=dy[valid], fmt="o", ms=4,
                    color=colour, lw=0.7, label=label)

ax.set_xlabel("wavelength / micron")
ax.set_ylabel("flux density / Jy")
ax.legend(fontsize="small", ncols=2)
plt.show()
```

## Draw joint posterior realizations correctly

The posterior has one row per *joint* realization. When drawing a total model:

1. choose one row index;
2. evaluate every desired component from that same row;
3. include components whose owners occur in the observation's contributors;
4. sum those spectra; and only then
5. draw the combined line.

Never select a separate posterior row for A and B. Their parameters can be
correlated by composite photometry and spectra, so independently selected rows
do not represent samples from the fitted joint posterior.

The disk-evaluation code in [Overplot disk models drawn from the
posterior](fit-output-plots.md#overplot-disk-models-drawn-from-the-posterior)
already selects parameters by full component ID. Put it in a function, call it
for each component using a common row index, and sum according to
`scope["contributors"]`.

## Limits for multiple systems

The package can contain one `limits-data/component-<n>.npz` per physical
contributor. Identify each archive from its metadata:

```python
from pathlib import Path

for path in sorted((package_dir / "limits-data").glob("component-*.npz")):
    arrays = load_npz(path)
    metadata = json.loads(arrays["metadata_json"].item())
    owner = metadata["owner_key"]
    print(path, targets[owner]["main_id"])
```

Run the limits plotting recipe once per archive. Disk posterior points in that
archive already correspond to its `owner_key`; do not combine A and B limits
unless that is the scientific comparison you intend to make.

## Common mistakes

- Do not use row order to decide which component is A or B; use
  `component_owner_keys`.
- Do not treat a composite observation as a third physical target. It is a
  scope whose `contributors` contains multiple owner keys.
- Do not sum component percentile curves to construct a percentile of the
  total. Use `prediction_median` for observed samples, or evaluate and sum
  complete posterior rows for a new spectral grid.
- Do not assume `component-1.npz` belongs to A; inspect `owner_key`.
