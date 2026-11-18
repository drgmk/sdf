# Fitting examples and science regressions

These deterministic examples exercise the same fitting core through a single
rawphot table and through resolved/composite binary packages. Their photometry
is generated directly from the `phoenix_sol` grid at `Teff = 6000 K` and
`log(g) = 4.5`, with fixed two-percent statistical and one-percent systematic
uncertainties. The selected bands have no measurement-level recalibration, so
the committed catalogue values map exactly back to the generated model fluxes.

- `single-star` is the reference one-component fit.
- `faint-binary` gives A the reference photometry, makes B 1000 times fainter,
  and supplies resolved A/B plus composite AB measurements.
- `equal-binary` duplicates the reference component and supplies resolved A/B
  plus composite AB measurements at twice the component flux.
- `spectrum-star` adds a portable 0.8–2.5 μm CSV spectrum to the reference
  photometry and exercises the fitted spectral-normalization parameter.

The committed IPAC and YAML files are ready to fit. Regenerate them after an
intentional model/filter change with:

```bash
python examples/fitting/generate.py
```

## One-file fit

Keep generated fit products out of the source tree:

```bash
sdf-fit examples/fitting/single-star/single-star-rawphot.txt \
  --output /tmp/sdf-single-example
```

With no model options, `sdf-fit` infers an automatic policy from the table
metadata and runs the normal evidence-guided search.
An exact one-model fit is also available:

```bash
sdf-fit examples/fitting/single-star/single-star-rawphot.txt \
  --output /tmp/sdf-single-exact --stellar phoenix_sol
```

## Binary fits

Copy a binary package before fitting so retained results and MultiNest work
files do not enter the repository:

```bash
cp -R examples/fitting/faint-binary /tmp/sdf-faint-binary
sdf-fits fit /tmp/sdf-faint-binary \
  --stellar 'Synthetic Faint A=phoenix_sol' \
  --stellar 'Synthetic Faint B=phoenix_sol'

cp -R examples/fitting/equal-binary /tmp/sdf-equal-binary
sdf-fits fit /tmp/sdf-equal-binary \
  --stellar 'Synthetic Equal A=phoenix_sol' \
  --stellar 'Synthetic Equal B=phoenix_sol'
```

The corresponding `truth.json` records the component parameter vectors. The
fast test suite verifies exact model/likelihood invariants at those values;
full MultiNest comparisons are deliberately an explicit science-regression
step rather than part of every unit-test run.

The spectrum package can be fitted directly with the package command:

```bash
cp -R examples/fitting/spectrum-star /tmp/sdf-spectrum-star
sdf-fits fit /tmp/sdf-spectrum-star \
  --stellar 'Synthetic Spectrum=phoenix_sol'
```

## Compare the pre-branch and current results

`science_regression.py` converts the old `Result` object and a current retained
candidate to the same small JSON schema. Run the legacy step with `PYTHONPATH`
pointing to a clean checkout of the pre-branch revision. Always copy the input
first because the old workflow writes MultiNest state beside it:

```bash
cp examples/fitting/single-star/single-star-rawphot.txt /tmp/legacy-rawphot.txt
PYTHONPATH=/path/to/pre-branch-sdf \
  python examples/fitting/science_regression.py legacy \
  /tmp/legacy-rawphot.txt --model phoenix_sol --output /tmp/legacy.json

sdf-fit examples/fitting/single-star/single-star-rawphot.txt \
  --output /tmp/current-single --stellar phoenix_sol
python examples/fitting/science_regression.py candidate \
  /tmp/current-single --output /tmp/current.json

python examples/fitting/science_regression.py compare \
  /tmp/legacy.json /tmp/current.json --output /tmp/comparison.json
```

The comparison reports fitted and derived stellar differences in units of the
combined posterior uncertainty, plus the evidence and chi-squared differences.
For evidence it retains the importance-sampling error but uses the more
conservative nested-sampling run uncertainty when normalizing independent
runs. MultiNest samples are intentionally not compared row by row.

For component invariants, select a binary owner and suppress system-level fit
statistics because the single and binary datasets are different:

```bash
python examples/fitting/science_regression.py candidate \
  /tmp/sdf-faint-binary --owner 'Synthetic Faint A' \
  --output /tmp/faint-a.json
python examples/fitting/science_regression.py compare \
  /tmp/current.json /tmp/faint-a.json --parameters-only
```

## Disposable comparisons for real stars

`real_single_regression.py` accepts one or more real rawphot files. It exports
the pre-branch revision into temporary storage, runs old and current exact fits
in isolated subprocesses, prints a compact normalized-difference table, and
deletes all native MultiNest and rendered results on exit. Only the optional
summary JSON is retained:

```bash
python examples/fitting/real_single_regression.py \
  /path/to/first-rawphot.txt /path/to/second-rawphot.txt \
  --output /tmp/real-single-summary.json
```

The checked-in case list spans Ross 128, HD 172051, and θ Leo. Repeat that
specific comparison against an SDB export with:

```bash
python examples/fitting/real_single_regression.py \
  --export-root /path/to/sdb/export \
  --output /tmp/real-single-summary.json
```

Package markers and spectrum pointers are removed from the temporary copies;
this check deliberately compares the photometry-only single-star likelihood.
