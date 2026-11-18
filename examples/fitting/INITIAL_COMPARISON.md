# Initial old/new fitting comparison

This records the first science-regression run made while preparing
the joint-fitting branch before its 2.0 merge. It is a smoke-test result, not a frozen set of
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

## Real single stars

The disposable real-star runner was then applied to the checked-in M/G/A case
list. No native MultiNest files, rendered pages, or complete fit results were
retained.

| Target | Teff | log(g) | norm | luminosity | radius | maximum |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Ross 128 | −0.002 σ | +0.051 σ | −0.034 σ | +0.000 σ | −0.027 σ | 0.051 σ |
| HD 172051 | −0.109 σ | −0.087 σ | +0.106 σ | −0.082 σ | +0.062 σ | 0.109 σ |
| θ Leo | +0.038 σ | −0.002 σ | −0.044 σ | −0.010 σ | −0.123 σ | 0.123 σ |

The evidence differences were at most 0.051 of the combined nested-sampling
run uncertainty, and the absolute chi-squared change was at most 0.110.

## Spectrum example

The portable spectrum package was run once through MultiNest and then deleted.
For input parameters of 6000 K, log(g) 4.5, and normalization −2.6, it recovered
Teff = 5994 (+14/−15) K, log(g) = 4.24 (+0.19/−0.31), and normalization
−2.5982 (+0.0053/−0.0047), with chi-squared 1.03. The fast regression test also
checks the combined photometry/spectrum likelihood exactly at the recorded
truth, including the spectral-normalization parameter.
