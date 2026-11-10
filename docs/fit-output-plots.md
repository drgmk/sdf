# Reproducing SED and limits plots

This guide turns the downloadable `sdf` arrays into Matplotlib figures. It
covers a single physical system; resolved multiple systems are treated in the
[multiple-system extension](fit-output-multiple-systems.md).

The plotting archives contain the reduced coordinates displayed by the site,
so reproducing a served plot does not require loading raw photometry, FITS
spectra, or an `sdf` model grid. Only the final section, which draws new model
spectra from posterior rows, requires an installed and configured `sdf`.

## Before starting

Download these links from the fit tab:

- **fit (JSON)**: identifies the fit and its other files;
- **posteriors (NPZ)**: needed only for new posterior calculations or model
  draws;
- **SED plot data (NPZ)**: observations, spectra, predictions, residuals, and
  median component spectra; and
- **limits plot data (NPZ)**: available from each limits tab.

The examples need NumPy and Matplotlib:

```shell
python -m pip install numpy matplotlib
```

If you have the exported package tree, load the selected candidate from its
`fits` directory. If you downloaded the three fit files individually, put
them together and set `candidate_dir` to that directory instead.

```python
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt

package_dir = Path("sdbid-v3-000552.54-414511.0")
candidate_dir = package_dir / "fits" / "phoenix_sol+modbb_disk_r--2ba2d639"
result = json.loads((candidate_dir / "result.json").read_text())

def load_npz(path):
    with np.load(path, allow_pickle=False) as archive:
        return {name: archive[name] for name in archive.files}

sed = load_npz(candidate_dir / result["sed_data"]["file"])
sed_metadata = json.loads(sed["metadata_json"].item())
```

## Reproduce the SED and residual panels

The served figure distinguishes photometry from observed spectra using
`data_scopes`. Each scope gives a half-open `[element_start, element_stop)`
slice into every observation-length array. This is important: a spectrum is
not a separate, hidden product in the plotting archive.

The following is a complete two-panel plot. It includes used and ignored
photometry, upper limits, model photometry, every observed spectral segment,
the predicted spectrum at those samples, component model spectra, and
residuals.

```python
wave = sed["wavelength_micron"]
flux = sed["observed_flux"]
error = sed["observed_uncertainty"]
prediction = sed["prediction_median"]
residual = sed["residual"]
ignored = sed["observed_ignore"].astype(bool)
upper_limit = sed["observed_upper_limit"].astype(bool)

fig, (ax, residual_ax) = plt.subplots(
    2, 1, figsize=(9, 6.6), sharex=True,
    gridspec_kw={"height_ratios": (3, 1), "hspace": 0.04},
)
colours = plt.get_cmap("tab10")

# Smooth median model spectra. The rows use the parallel component arrays.
model_wave = sed["spectrum_wavelength_micron"]
component_spectra = sed["component_spectrum_median"]
component_ids = sed["component_ids"].astype(str)
component_roles = sed["component_roles"].astype(str)
for index, (component_id, role, model_flux) in enumerate(zip(
    component_ids, component_roles, component_spectra
)):
    ax.plot(
        model_wave, model_flux,
        ls="--" if role == "stellar" else ":",
        lw=1.2, alpha=0.7, color=colours(index),
        label=f"{role}: {component_id}",
    )
ax.plot(
    model_wave, component_spectra.sum(axis=0),
    color="black", lw=2.0, alpha=0.8, label="total model",
)

for scope_index, scope in enumerate(sed_metadata["data_scopes"]):
    section = slice(scope["element_start"], scope["element_stop"])
    x = wave[section]
    y = flux[section]
    dy = error[section]
    model_y = prediction[section]
    r = residual[section]
    omit = ignored[section]
    limit = upper_limit[section]
    colour = colours(scope_index % 10)
    label = scope.get("instrument") or scope["observation_id"]
    finite = np.isfinite(x) & (x > 0) & np.isfinite(y) & (y > 0)

    if scope["kind"] == "ObsSpectrum":
        # Spectral normalization fitted by sdf is already included in y and dy.
        used = finite & ~omit
        ax.plot(x[used], y[used], color=colour, lw=1.0, label=label)
        ax.plot(x[used], model_y[used], color=colour, lw=1.7, ls="--")
        good_r = used & np.isfinite(r)
        residual_ax.plot(x[good_r], r[good_r], color=colour, lw=0.9)
        continue

    detection = finite & ~limit
    used = detection & ~omit
    excluded = detection & omit
    limits = finite & limit

    ax.errorbar(
        x[used], y[used], yerr=dy[used], fmt="o", ms=4,
        color=colour, lw=0.7, label=label,
    )
    ax.scatter(x[excluded], y[excluded], s=18, color="0.7", zorder=3)
    ax.scatter(x[limits], y[limits], s=28, marker="v", color=colour, zorder=3)

    model_ok = np.isfinite(x) & (x > 0) & np.isfinite(model_y) & (model_y > 0)
    ax.scatter(
        x[model_ok], model_y[model_ok], s=42,
        facecolors="none", edgecolors=colour, linewidths=1.0,
    )

    for mask, marker, point_colour in (
        (used, "o", colour),
        (excluded, "o", "0.7"),
        (limits & ~omit, "v", colour),
        (limits & omit, "v", "0.7"),
    ):
        valid = mask & np.isfinite(r)
        residual_ax.scatter(x[valid], r[valid], s=18, marker=marker,
                            color=point_colour)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylabel("flux density / Jy")
ax.legend(fontsize="small", ncols=2)
ax.grid(alpha=0.18)

residual_ax.axhline(0, color="0.5", lw=0.8)
residual_ax.axhline(-3, color="0.7", lw=0.8, ls="--")
residual_ax.axhline(3, color="0.7", lw=0.8, ls="--")
finite_residuals = residual[np.isfinite(residual)]
extent = max(4.0, np.max(np.abs(finite_residuals)) * 1.08)
residual_ax.set_ylim(-extent, extent)
residual_ax.set_xscale("log")
residual_ax.set_xlabel("wavelength / micron")
residual_ax.set_ylabel("residual / sigma")
residual_ax.grid(alpha=0.18)
plt.show()
```

The example intentionally uses the fitted `observed_flux`. For a spectral
segment this includes its fitted multiplicative normalization. Use
`observed_flux_native` and `observed_uncertainty_native` to inspect the
original imported values instead.

### What exactly matches the served plot?

- `prediction_median` supplies the open model points and sampled spectral
  prediction.
- `component_spectrum_median` supplies the smooth stellar and disk curves.
- `residual` is already in uncertainty units; do not recompute it as a simple
  quotient because upper limits and ignored data require special handling.
- `observed_ignore` and `observed_upper_limit` are independent flags.
- Prediction intervals are available as `prediction_p16` and
  `prediction_p84` if a shaded interval is useful, although the standard
  served SED does not currently draw it.

## Reproduce the limits plots

One limits archive supplies both the temperature and blackbody-radius views.
The blue region begins at `envelope`, the lowest detection-limit curve at each
temperature. Solid curves are measurement limits and dashed curves are
catalogue upper limits.

```python
limits = load_npz(package_dir / "limits-data" / "component-1.npz")

views = [("temperature_k", "disk temperature / K")]
if "radius_au" in limits:
    views.append(("radius_au", "blackbody radius / au"))
fig, axes = plt.subplots(1, len(views), figsize=(6 * len(views), 5),
                         sharey=True, squeeze=False)

for ax, (coordinate_name, xlabel) in zip(axes[0], views):
    x = limits[coordinate_name]
    colours = plt.get_cmap("tab20")
    ax.fill_between(
        x, limits["envelope"], 2e-1,
        color="#dbeafe", alpha=0.65, label="detectable",
    )

    for index, (band, basis) in enumerate(zip(
        limits["filters"].astype(str), limits["basis"].astype(str)
    )):
        ax.plot(
            x, limits["curves"][:, index],
            color=colours(index), lw=1.5,
            ls="--" if basis == "upper limit" else "-",
            label=band,
        )

    for curve, label, dash in zip(
        limits["reference_curves"],
        limits["reference_labels"].astype(str),
        limits["reference_dashes"].astype(str),
    ):
        line_style = {
            "solid": "-", "dashed": "--", "dotted": ":", "dotdash": "-."
        }.get(dash, "-")
        ax.plot(x, curve, color="0.55", lw=1.2,
                ls=line_style, label=label)

    disk_x_name = (
        "disk_temperature_k" if coordinate_name == "temperature_k"
        else "disk_radius_au"
    )
    for component_id, disk_x, fractional_luminosity in zip(
        limits["disk_component_ids"].astype(str),
        limits[disk_x_name],
        limits["disk_fractional_luminosity"],
    ):
        valid = (
            np.isfinite(disk_x) & (disk_x > 0) &
            np.isfinite(fractional_luminosity) &
            (fractional_luminosity > 0)
        )
        ax.scatter(
            disk_x[valid], fractional_luminosity[valid],
            s=7, alpha=0.12, label=f"fitted {component_id}",
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-7, 2e-1)
    ax.set_xlabel(xlabel)
    ax.grid(alpha=0.18)

axes[0, 0].set_ylabel("fractional luminosity")
axes[0, -1].legend(fontsize="x-small", ncols=2)
plt.tight_layout()
plt.show()
```

If `radius_au` is absent, the stellar luminosity needed for the radius
conversion was unavailable. Plot only the temperature view in that case. For
a multiple system there can be one limits archive per physical contributor;
read `owner_key` from each archive's `metadata_json` rather than assuming that
`component-1.npz` is the primary.

## Overplot disk models drawn from the posterior

`sed-data.npz` intentionally contains only the median model spectrum. To draw
many correlated posterior realizations, use complete rows from the posterior
NPZ and evaluate the corresponding `sdf` model grid. This preserves
correlations among temperature, turnover wavelength, emissivity slope, and
normalization.

This example chooses the first disk component from the metadata and draws 50
of its spectra. It works for any disk grid whose fitted parameters and model
file are available to the local `sdf` installation:

```python
from sdf.model import SpecModel

posterior = load_npz(candidate_dir / result["posterior"]["file"])
names = posterior["parameter_names"].astype(str)
samples = posterior["samples"]
columns = {name: index for index, name in enumerate(names)}

disk = next(
    component for component in result["components"]
    if component["role"] == "disk"
)
grid = SpecModel.read_model(disk["model_name"])

# SpecModel.fnujy expects its grid parameters followed by log10 normalization.
wanted_names = [
    f"{disk['component_id']}.{parameter}"
    for parameter in (*grid.parameters, "norm")
]
wanted_columns = [columns[name] for name in wanted_names]

rng = np.random.default_rng(12345)
weights = posterior.get("weights")
probability = None if weights is None else weights / weights.sum()
draws = rng.choice(len(samples), size=min(50, len(samples)),
                   replace=True, p=probability)

fig, ax = plt.subplots(figsize=(8, 5))
for row in draws:
    model_flux = grid.fnujy(samples[row, wanted_columns])
    ax.loglog(grid.wavelength, model_flux, color="tab:orange", alpha=0.10)

# Put the observed SED on top for context.
photometry = np.array([
    scope["kind"] == "Photometry"
    for scope in sed_metadata["data_scopes"]
    for _ in range(scope["element_stop"] - scope["element_start"])
])
shown = (
    photometry & ~sed["observed_ignore"].astype(bool) &
    np.isfinite(sed["wavelength_micron"]) &
    (sed["wavelength_micron"] > 0) &
    (sed["observed_flux"] > 0)
)
ax.scatter(sed["wavelength_micron"][shown], sed["observed_flux"][shown],
           s=16, color="black", zorder=3, label="data")
ax.set_xlabel("wavelength / micron")
ax.set_ylabel("flux density / Jy")
ax.legend()
plt.show()
```

Do not draw each disk parameter independently from its marginal distribution:
that creates combinations the posterior did not support. Select rows, as
above. In a two-disk fit, repeat the component selection for each disk; sum
spectra from the *same posterior row* when drawing a total realization.

## Related reference

- [Using `sdf` fit outputs](fit-outputs.md) describes every artifact and shows
  how to extract summaries and posterior columns.
- [Plotting a multiple system](fit-output-multiple-systems.md) explains owner
  keys, contributors, resolved observations, and composite observations.
