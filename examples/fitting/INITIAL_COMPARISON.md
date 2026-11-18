# Initial old/new fitting comparison

This records the first science-regression run made while preparing
`codex/joint-fitting` for merge. It is a smoke-test result, not a frozen set of
golden posterior samples.

- Date: 2026-11-18
- Pre-branch revision: `c8dd745` (`master`)
- Joint-fitting core parent revision: `9d45f4e`
- Sampler: MultiNest 3.10, 200 live points
- Stellar model: `phoenix_sol`
- Input: the committed deterministic `single-star` fixture

The deterministic checks first established identical model dimensions, priors,
residuals, and log likelihoods at fixed parameters. Independent old/new
MultiNest runs then gave:

| Quantity | New − old, combined posterior σ |
| --- | ---: |
| Teff | −0.039 |
| log(g) | +0.082 |
| normalization | +0.002 |
| luminosity | −0.097 |
| radius | −0.001 |

Chi-squared changed from 0.774 to 0.698. Log evidence changed by +0.323, or
0.868 of the combined nested-sampling run uncertainty. The smaller
importance-sampling evidence errors are retained in the result metadata but
are not a realistic normalization for comparing independent runs.

The current joint fitter was also run on the two binary fixtures:

| Comparison | Largest fitted/derived difference |
| --- | ---: |
| Single star vs faint-binary A | 0.268 σ |
| Single star vs equal-binary A | 0.138 σ |
| Equal-binary A vs B | 0.106 σ |

These results support the intended invariants: adding a constrained component
that is 1000 times fainter does not materially alter the bright component, and
duplicated equal components recover statistically identical stellar results.
They do not replace broader tests of calibrated input bands, spectra, automatic
model selection, or real systems.
