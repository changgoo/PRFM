# Core-suite correlation results

**Status:** model-level correlation audit through 2026-08-26.

## Purpose

This record quantifies the correlations summarized in
[the core-suite synthesis](core-suite-analysis-summary.md). It reports
association, not causal or independent parameter effects.

The validation quantities use all 32 models and the 400–600 Myr history
window. Density-PDF and power-spectrum quantities use the 31 clean projection
sets and temporal medians over 200–600 Myr. Each table reports Spearman's rank
coefficient `rho_s`. The `q` values are
Benjamini–Hochberg false-discovery-rate corrections across all entries in that
table.

## Design covariance

The derived variables are intentionally functions of the primary design
variables and must not be interpreted as independent axes:

| Pair | rho_s |
|---|---:|
| Sigma_star–rho_star | +0.908 |
| Omega–kappa | +0.974 |

The strong first correlation reflects the sampled `Sigma_star`
range and the definition `rho_star = Sigma_star/(2 H_star)`. The
second follows from `kappa = sqrt(2(2-q)) Omega`. Consequently,
similar correlations with the members of each pair generally identify a
stellar-gravity or rotation axis, not a unique controlling variable.

## Observational-validation outcomes

| Outcome | Predictor | rho_s | q |
|---|---|---:|---:|
| Sigma_SFR | Sigma_star | +0.387 | 0.071 |
| Sigma_SFR | H_star | -0.042 | 0.893 |
| Sigma_SFR | Omega | -0.257 | 0.310 |
| Sigma_SFR | q | -0.076 | 0.816 |
| Sigma_SFR | kappa | -0.210 | 0.409 |
| Sigma_SFR | rho_star | +0.440 | 0.047 |
| f_mol | Sigma_star | +0.085 | 0.816 |
| f_mol | H_star | +0.200 | 0.409 |
| f_mol | Omega | +0.794 | ≤0.001 |
| f_mol | q | +0.385 | 0.071 |
| f_mol | kappa | +0.684 | ≤0.001 |
| f_mol | rho_star | +0.020 | 0.914 |

The clearest `Sigma_SFR` association is with the stellar-gravity
axis, especially `rho_star`. There is no evidence in this 32-model
sample for a monotonic `Sigma_SFR` dependence on
`H_star`, rotation, or `q` considered one at a time.

The apparent `f_mol` response instead follows the rotation axis.
Because the simulated molecular fraction is low by about 1.35 dex and is not
resolution-converged, this trend is useful for targeting convergence tests but
is not yet a physical calibration.

Source data and measurement definitions are in
[validation diagnostics](validation-diagnostics.md).

## Density-PDF widths versus environment

| Width | Predictor | rho_s | q |
|---|---|---:|---:|
| sigma_delta | kappa | -0.417 | 0.034 |
| sigma_delta | rho_star | +0.670 | ≤0.001 |
| sigma_delta | Sigma_SFR | +0.552 | 0.003 |
| sigma_delta | Omega | -0.415 | 0.034 |
| sigma_delta | q | +0.163 | 0.414 |
| sigma_s | kappa | -0.328 | 0.094 |
| sigma_s | rho_star | +0.717 | ≤0.001 |
| sigma_s | Sigma_SFR | +0.557 | 0.003 |
| sigma_s | Omega | -0.325 | 0.094 |
| sigma_s | q | +0.152 | 0.414 |

Both width definitions show the same primary result: column density becomes
more intermittent in higher-`rho_star`, higher-`Sigma_SFR`
models. The linear-density width also narrows moderately along the rotation
axis. Neither width shows a supported monotonic relation with `q`
alone.

The widths correlate very strongly with absolute velocity and thermal support:
`rho_s = 0.79–0.93` for the three velocity-dispersion components,
`0.87–0.90` for the three-dimensional velocity dispersion, and
`0.89–0.91` for the thermal speed. By contrast, neither width
correlates with the three-dimensional sonic Mach number
(`rho_s = 0.09` and `0.02`). Within this suite, PDF
broadening therefore tracks the amplitudes of the support variables more
directly than sonic Mach number alone. These are co-variations, not evidence
that thermal speed independently causes the broader PDFs.

Source data and definitions are in the
[density-PDF analysis](https://github.com/changgoo/tigress_ncr_tools/blob/analysis/ncr-suite/docs/density_pdf.md).

## Power-spectrum measures versus environment

| Measure | Predictor | rho_s | q |
|---|---|---:|---:|
| L_in | kappa | -0.808 | ≤0.001 |
| L_in | rho_star | -0.439 | 0.034 |
| L_in | Sigma_SFR | +0.208 | 0.355 |
| L_in | Omega | -0.817 | ≤0.001 |
| L_in | q | +0.052 | 0.839 |
| alpha | kappa | +0.184 | 0.370 |
| alpha | rho_star | -0.300 | 0.218 |
| alpha | Sigma_SFR | -0.214 | 0.355 |
| alpha | Omega | +0.222 | 0.355 |
| alpha | q | +0.017 | 0.928 |
| A_2 | kappa | +0.575 | 0.002 |
| A_2 | rho_star | -0.199 | 0.355 |
| A_2 | Sigma_SFR | -0.634 | ≤0.001 |
| A_2 | Omega | +0.667 | ≤0.001 |
| A_2 | q | +0.281 | 0.236 |

The strongest power-spectrum result is a smaller unrestricted integral scale
at faster rotation. The integral scale also decreases with stellar midplane
density but shows no supported monotonic relation with `Sigma_SFR`
or `q`. The fitted slope has no supported correlation with any
tested environmental quantity; treating `alpha` as a primary
environmental diagnostic is therefore not justified by the present suite.

The global quadrupole is larger in faster-rotating and lower-SFR models. Its
lack of a clear `q` correlation indicates that this bivariate trend
follows `Omega`/`kappa`, rather than `q` by
itself. The absolute `L_in` and `A_2` values remain
box-size sensitive, so these statements apply to relative trends within the
common-box suite. Fixed-band scale measures and local anisotropy diagnostics
should be used to test whether the correlations survive a more robust
statistic.

Source data and definitions are in the
[power-spectrum analysis](https://github.com/changgoo/tigress_ncr_tools/blob/analysis/ncr-suite/docs/density_power_spectrum.md);
convergence limitations and alternatives are in the
[box-size study](https://github.com/changgoo/tigress_ncr_tools/blob/analysis/ncr-suite/docs/box_size_density_power_spectrum.md).

## Cross-diagnostic associations

| PDF width | Power-spectrum measure | rho_s | q |
|---|---|---:|---:|
| sigma_delta | L_in | +0.254 | 0.201 |
| sigma_delta | alpha | -0.473 | 0.011 |
| sigma_delta | A_2 | -0.623 | 0.001 |
| sigma_s | L_in | +0.192 | 0.301 |
| sigma_s | alpha | -0.513 | 0.006 |
| sigma_s | A_2 | -0.577 | 0.002 |

Models with broader column-density PDFs tend to have flatter fitted spectra
(smaller `alpha`) and weaker global quadrupole anisotropy, while PDF
width does not track the unrestricted integral scale. The slope and quadrupole
relations are useful within-suite empirical patterns, but inherit the
resolution and box-size qualifications of those measures.

## Interpretation limits

- Rank correlations test monotonic association and do not establish causality.
- FDR correction limits false discoveries among the displayed bivariate tests;
  it does not remove design covariance.
- `Sigma_star` versus `rho_star` and
  `Omega` versus `kappa` cannot be cleanly separated
  with this design.
- The 31–32 model sample is sufficient to identify strong trends but gives
  limited leverage for weak, nonlinear, or interaction effects.
- Temporal percentile ranges are not independent model points and were not
  treated as extra samples in these tests.
- Independent effects should be assessed with a multivariate emulator,
  hierarchical model, partial dependence, or controlled ablations.
