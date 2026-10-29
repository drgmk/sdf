# Using SDF fit outputs

Each retained SDF fit has two downloadable data products:

- `result.json` is a compact, human-readable description and summary of the
  fit; and
- `posterior-<hash>.npz` contains the numerical posterior samples and arrays
  used to inspect or reconstruct the fit.

The JSON is usually the easiest place to start. Use the NPZ when you need the
sample distributions, correlations, component predictions, or spectra. The
NPZ is self-describing and can also be used by itself.

## Finding the files

Within an exported target package, every model configuration has its own
directory:

```text
<package>/
  fits/
    selection.json
    <model-name>--<candidate-hash>/
      result.json
      posterior-<posterior-hash>.npz
      corner.png
      distributions.png
```

The fit page links directly to the JSON, NPZ, and diagnostic PNG files. If you
are working with a downloaded candidate directory, load the JSON first and use
it to find the content-addressed posterior filename:

```python
from pathlib import Path
import json

candidate_dir = Path("phoenix_sol+modbb_disk_r--2ba2d639")
result = json.loads(
    (candidate_dir / "result.json").read_text(encoding="utf-8")
)
posterior_path = candidate_dir / result["posterior"]["file"]

print(result["display_name"])
print(posterior_path)
```

`fits/selection.json` says which retained candidates are displayed and which
one is the current default. It does not duplicate the fit results.

## Reading `result.json`

The most useful top-level entries are:

| Entry | Contents |
| --- | --- |
| `candidate_id` | Stable identity of the exact component/model combination |
| `display_name` | Short model name used in the UI |
| `completed_at` | Fit completion time |
| `components` | Model components, physical owners, roles, and model names |
| `posterior` | Posterior filename, evidence, sampler, and sample count |
| `analysis` | Parameter and physical-result summaries, fit statistics, and data scopes |
| `policy` | Model-search configuration used for the target or components |
| `provenance` | Available code and input provenance |

For example:

```python
p = result["posterior"]
print(f"sampler: {p['sampler']}")
print(f"samples retained: {p['n_samples']}")
print(f"ln Z = {p['evidence']:.2f} +/- {p['evidence_error']:.2f}")

analysis = result["analysis"]
print(f"chi-squared = {analysis['chi_squared']:.1f}")
print(f"residuals = {analysis['n_residuals']}")
```

More positive log-evidence is preferred when comparing candidates for the same
input data. Evidence should not be compared between unrelated systems or fits
that used different data.

### Components and parameter summaries

Component IDs are deliberately explicit because one system may contain more
than one star or disk. Do not infer ownership by splitting these strings; use
the `owner_key`, `role`, and `model_name` fields:

```python
for component in result["components"]:
    print(
        component["role"],
        component["model_name"],
        "belongs to",
        component["owner_key"],
    )
```

The fitted-parameter summaries are under `analysis.parameters`. Each summary
is the median and the distances from the median to the 16th and 84th
percentiles:

```python
for name, summary in analysis["parameters"].items():
    print(
        f"{name}: {summary['median']:.5g} "
        f"-{summary['minus']:.2g} +{summary['plus']:.2g}"
    )
```

The same summaries are grouped by model component in `analysis.components`.
Physical results are grouped by target in `analysis.targets`, including such
quantities as stellar radius and luminosity, disk temperature, fractional
luminosity, and blackbody radius when they can be calculated:

```python
for target_key, target in analysis["targets"].items():
    print(target["main_id"])
    for star in target.get("stellar_components", []):
        luminosity = star.get("luminosity_lsun")
        if luminosity:
            print("  Lstar/Lsun:", luminosity["median"])
    for disk in target.get("disk_components", []):
        fraction = disk.get("fractional_luminosity")
        if fraction:
            print("  Ldisk/Lstar:", fraction["median"])
```

Use the NPZ rather than these summaries if you need correlations or want to
calculate different credible intervals.

## Reading the posterior NPZ

The file contains ordinary NumPy arrays and no pickled Python objects, so it
can be loaded safely with `allow_pickle=False`:

```python
import json
import numpy as np

with np.load(posterior_path, allow_pickle=False) as archive:
    arrays = {name: archive[name] for name in archive.files}

print(arrays["artifact_type"].item())   # sdf.posterior
print(arrays["schema_version"].item())  # currently 1
metadata = json.loads(arrays["metadata_json"].item())
print(metadata["description"])
print(metadata["array_axes"])
```

Inspecting `archive.files`, the shapes, and dtypes is a useful first step for
any artifact:

```python
for name, values in arrays.items():
    print(f"{name:40s} {str(values.shape):16s} {values.dtype}")
```

Not every fit needs every optional array. Test for a name with
`if name in arrays` when writing code that will process many systems.

## Fitted and derived samples

`samples` contains the parameters sampled directly by the fitter. Its columns
are named, in order, by `parameter_names`:

```python
fitted_names = arrays["parameter_names"].astype(str)
fitted_samples = arrays["samples"]

fitted = {
    name: fitted_samples[:, column]
    for column, name in enumerate(fitted_names)
}

temperature = fitted[
    "simbad:HD 123:stellar:phoenix_sol:1.Teff"
]
```

The exact component IDs vary by system, so general code should normally select
names from the supplied metadata rather than hard-code the example above.

`derived_samples` contains quantities calculated from the fitted samples. It
also propagates target information such as parallax uncertainty. Its columns
are described by several parallel arrays:

```python
derived_samples = arrays["derived_samples"]
derived_names = arrays["derived_parameter_names"].astype(str)
derived_labels = arrays["derived_parameter_labels"].astype(str)
derived_units = arrays["derived_parameter_units"].astype(str)
derived_owners = arrays["derived_parameter_owner_keys"].astype(str)
derived_components = arrays["derived_parameter_component_ids"].astype(str)
derived_roles = arrays["derived_parameter_roles"].astype(str)

for column in range(len(derived_names)):
    q16, median, q84 = np.percentile(
        derived_samples[:, column], [16, 50, 84]
    )
    unit = f" {derived_units[column]}" if derived_units[column] else ""
    print(
        f"{derived_labels[column]} ({derived_owners[column]}): "
        f"{median:.4g} -{median-q16:.2g} +{q84-median:.2g}{unit}"
    )
```

The rows of `samples` and `derived_samples` correspond. Preserve that row
alignment when examining correlations. For example:

```python
def named_column(values, names, wanted):
    matches = np.flatnonzero(names == wanted)
    if len(matches) != 1:
        raise KeyError(f"expected one {wanted!r}, found {len(matches)}")
    return values[:, matches[0]]

radius_name = derived_names[derived_labels == "stellar radius"][0]
luminosity_name = derived_names[derived_labels == "stellar luminosity"][0]

radius = named_column(derived_samples, derived_names, radius_name)
luminosity = named_column(derived_samples, derived_names, luminosity_name)

import matplotlib.pyplot as plt
plt.scatter(radius, luminosity, s=5, alpha=0.25)
plt.xlabel("stellar radius / R_sun")
plt.ylabel("stellar luminosity / L_sun")
plt.show()
```

For multiple systems or components, select with metadata rather than labels
alone. For example, all disk-derived quantities belonging to one target are:

```python
target_key = "simbad:HD 123"
use = (derived_roles == "disk") & (derived_owners == target_key)

for name, label, unit, values in zip(
    derived_names[use],
    derived_labels[use],
    derived_units[use],
    derived_samples[:, use].T,
):
    print(name, label, unit, np.median(values))
```

An optional Pandas view is convenient for exploratory work:

```python
import pandas as pd

fitted_frame = pd.DataFrame(fitted_samples, columns=fitted_names)
derived_frame = pd.DataFrame(derived_samples, columns=derived_names)
```

Current MultiNest artifacts normally contain capped equal-weight samples, so
ordinary NumPy percentiles are appropriate. If a future or alternative
sampler writes a `weights` array, use it when calculating distributions. A
`log_likelihoods` array, when present, has one value per retained sample.

## Photometry, predictions, and residuals

The observation axis is described by arrays of equal length:

```python
wavelength = arrays["wavelength_micron"]
filters = arrays["filter"].astype(str)
observed = arrays["observed_flux"]
uncertainty = arrays["observed_uncertainty"]
upper_limit = arrays["observed_upper_limit"].astype(bool)
ignored = arrays["observed_ignore"].astype(bool)

model = arrays["prediction_median"]
model_p16 = arrays["prediction_p16"]
model_p84 = arrays["prediction_p84"]
residual_sigma = arrays["residual"]
```

`observed_flux` and `observed_uncertainty` are in Jy. For fitted spectra they
include the fitted segment normalization used by the likelihood;
`observed_flux_native` and `observed_uncertainty_native` preserve the imported
values. `residual` is in uncertainty units.

A minimal photometry plot is:

```python
finite = np.isfinite(wavelength) & (wavelength > 0)
detection = finite & ~upper_limit

plt.errorbar(
    wavelength[detection & ~ignored],
    observed[detection & ~ignored],
    yerr=uncertainty[detection & ~ignored],
    fmt=".",
    label="used data",
)
plt.scatter(
    wavelength[detection & ignored],
    observed[detection & ignored],
    color="0.7",
    label="ignored data",
)
plt.scatter(
    wavelength[finite & upper_limit],
    observed[finite & upper_limit],
    marker="v",
    label="upper limits",
)
plt.scatter(wavelength[finite], model[finite], facecolors="none", label="model")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("wavelength / micron")
plt.ylabel("flux density / Jy")
plt.legend()
plt.show()
```

The prediction intervals are calculated from a deterministic posterior subset,
currently no more than 256 rows. `analysis.prediction_samples` records how many
were used.

### Contributions from individual model components

The first axis of each `component_prediction_*` array follows the parallel
`component_*` metadata arrays:

```python
component_ids = arrays["component_ids"].astype(str)
component_owners = arrays["component_owner_keys"].astype(str)
component_roles = arrays["component_roles"].astype(str)
component_models = arrays["component_model_names"].astype(str)
component_flux = arrays["component_prediction_median"]

for component_id, owner, role, model_name, flux in zip(
    component_ids,
    component_owners,
    component_roles,
    component_models,
    component_flux,
):
    print(component_id, owner, role, model_name)
    plt.scatter(wavelength[finite], flux[finite], label=component_id)
```

These are contributions to each observed element. `prediction_median` is the
median total prediction; it is preferable to summing independently computed
component percentiles.

## Model spectra

Median component spectra are stored on a common wavelength grid:

```python
wave = arrays["spectrum_wavelength_micron"]
component_spectra = arrays["component_spectrum_median"]

for component_id, spectrum in zip(component_ids, component_spectra):
    plt.loglog(wave, spectrum, label=component_id)

plt.loglog(wave, component_spectra.sum(axis=0), color="black", label="total")
plt.xlabel("wavelength / micron")
plt.ylabel("flux density / Jy")
plt.legend()
plt.show()
```

This array contains the fitted model spectra, not the original observed
spectral products. Observed photometry and spectral samples share the
observation arrays described above.

To divide that axis into its original photometry or spectral products, use
`analysis.data_scopes` from `result.json`:

```python
for scope in result["analysis"]["data_scopes"]:
    section = slice(scope["element_start"], scope["element_stop"])
    print(scope["kind"], scope["observation_id"], observed[section].shape)
    if scope["kind"] == "Spectrum":
        print(scope.get("instrument"), scope.get("segment_id"))
```

## Comparing retained candidates

To make a small evidence table for all retained configurations in a package:

```python
rows = []
for path in sorted((candidate_dir.parent).glob("*/result.json")):
    candidate = json.loads(path.read_text(encoding="utf-8"))
    rows.append((
        candidate["posterior"]["evidence"],
        candidate["posterior"]["evidence_error"],
        candidate["display_name"],
        path.parent,
    ))

for evidence, error, name, path in sorted(
    rows, key=lambda row: row[0], reverse=True
):
    print(f"{evidence:9.2f} +/- {error:.2f}  {name:30s}  {path.name}")
```

The public default is not necessarily the numerically highest-evidence model:
SDF chooses the simplest model supported within its evidence threshold, and an
administrator may override that choice. Read `fits/selection.json` to
reproduce the display/default selection.

## Important details

- Fitted and derived columns must always be identified by their supplied name
  arrays; do not guess their order.
- Fitted parameters retain the parameterization used by the model. For
  example, `log_Temp` and `log_lam0` remain logarithmic. Derived temperature
  and luminosity arrays are in the physical units given alongside them.
- Posterior files retain a compact sample, normally capped by SDF's
  `n_samples_max`; they are not the sampler's complete native output.
- The JSON summaries use the 16th, 50th, and 84th percentiles and can be
  asymmetric.
- `observed_ignore` and `observed_upper_limit` are distinct. Ignored data are
  retained for inspection but did not constrain the fit.
- Component predictions and spectra are labelled by the parallel
  `component_*` arrays. Do not assume that row 0 is always the primary star.
- `metadata_json` makes the NPZ independently understandable, while
  `result.json` remains the most convenient source for summaries and fit
  configuration.

Artifacts written before the self-describing labels were introduced can be
updated without rerunning the fit:

```shell
sdf-fits refresh-posteriors /path/to/package
sdf-fits refresh-posteriors /path/to/export-root --all
```
