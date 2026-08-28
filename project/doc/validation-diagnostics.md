# Core-suite SFR and molecular-fraction validation

The core suite can be compared directly with its parent PHANGS gas band using
two outcomes that were deliberately excluded from the simulation design:

- star formation rate surface density, `Sigma_SFR`; and
- molecular gas fraction, `Sigma_mol / Sigma_gas`.

`project/scripts/plot_validation_diagnostics.py` plots each outcome against the
four environmental quantities varied by this fixed-gas suite
(`Sigma_star`, `H_star`, `Omega`, and `qshear`) and two derived quantities,

\[
\kappa = \sqrt{2(2-q)}\,\Omega, \qquad
\rho_* = \frac{\Sigma_*}{2H_*}.
\]

The same formulae are applied to the simulations and PHANGS parent pixels. In
particular, `rho_star` here is the density implied by the simulation input
pair (`Sigma_star`, `H_star`); it is not the PHANGS megatable's independent
`rho_star_mp` column.

## Outcome definitions

- Simulated `Sigma_SFR` uses every `sfr40` history sample in 400--600 Myr by
  default. The overlaid point is its trapezoidal time mean and the bar is its
  temporal 16th--84th percentile interval. The independent reducer value
  `sfr40_mean` is retained in the joined output as a cross-check.
- Observed `Sigma_SFR` is the PHANGS `Sigma_SFR_HaW4recal` tracer by default.
- Simulated molecular fraction uses every native volume-integrated conserved
  abundance `2 * scalar3 / mass`. `scalar3` is H2 in the five-scalar NCR build
  used by this suite. The point is the trapezoidal time mean over 400--600 Myr
  and the bar shows the temporal 16th--84th percentiles.
- Observed molecular fraction is the canonical Gaussian-aperture
  `Sigma_mol / Sigma_gas`. PHANGS molecular non-detections were encoded as zero
  when the parent table was built and are omitted from the logarithmic panels;
  the legend reports the number of positive detections.

Blue point clouds in every panel show the **full temporal distributions**, not
only their summary ranges. The molecular diagnostic uses the conserved
history integral rather than a resampled projection. This is the most direct global mass fraction and avoids
viewing-angle and pixel-remapping effects.

## Reproduce

```bash
PY=/home/changgoo/.conda/envs/pyathena/bin/python
$PY project/scripts/plot_validation_diagnostics.py \
    project/output/design_Sgas10.0_n0032.csv \
    --suite /tigress/changgoo/anvil/TIGRESS-NCR-suite \
    --summary /tigress/changgoo/anvil/TIGRESS-NCR-suite/prfm_diagnostics/prfm_model_summary.csv \
    --output-dir project/output \
    --figure-dir project/doc/paper1-suite/paper-figures
```

Products:

- `project/output/design_Sgas10.0_n0032_validation_summary.csv`;
- `project/output/design_Sgas10.0_n0032_scatter_budget.csv`;
- `project/doc/paper1-suite/paper-figures/design_Sgas10.0_n0032_sigma_sfr_vs_environment.png`;
- `project/doc/paper1-suite/paper-figures/design_Sgas10.0_n0032_fmol_vs_environment.png`.

These are validation plots, not independent one-parameter experiments. The
core design retains PHANGS covariances among its inputs, so an apparent trend
with one x-axis may partly reflect correlated variation in another input. The
within-box trends are nevertheless physically meaningful suite diagnostics;
causal attribution should use an emulator, partial dependence, or a controlled
ablation design.

## Current core-suite result

The current 32-model, 400--600 Myr comparison gives:

| Outcome | Suite median | PHANGS positive median | PHANGS scatter [dex] | Environment [dex] | Temporal [dex] | Combined suite [dex] | Suite/PHANGS variance |
|---|---:|---:|---:|---:|---:|---:|---:|
| `Sigma_SFR` | 0.00292 | 0.00280 | 0.271 | 0.378 | 0.173 | 0.415 | 2.35 |
| `Sigma_mol / Sigma_gas` | 0.0212 | 0.480 | 0.263 | 0.260 | 0.219 | 0.339 | 1.67 |

The SFR normalization is recovered very well, but the suite produces more log
scatter than the PHANGS parent band. About 83% of the modeled SFR variance is
between models and 17% is temporal. For molecular fraction, the modeled
variance divides roughly 58% between models and 42% within models. However,
the molecular normalization is low by about 1.35 dex, so matching or exceeding
the observed scatter does not constitute recovery of the molecular-fraction
distribution. The PHANGS molecular scatter here uses 746 positive detections;
45 zero-valued non-detections are omitted by the log transform.

## Scatter budget

For each positive outcome, the script works in `log10` space and gives every
model equal weight. It decomposes the suite variance as

\[
\sigma_{\rm suite}^2 =
\sigma_{\rm environment}^2 + \langle\sigma_{\rm temporal}^2\rangle,
\]

where the first term is the variance of model time means and the second is the
mean within-model temporal variance. The output table reports both components,
their quadrature sum, the PHANGS parent-band scatter, and the suite/PHANGS
variance ratio. This ratio is a scatter-coverage diagnostic, not a claim that
the model has explained individual PHANGS pixels. The finite 0.6-dex gas band
also means that the between-model term includes residual `Sigma_gas` variation
in addition to the four plotted environmental inputs. A conditional emulator
or a narrower gas slice is needed for a strictly fixed-`Sigma_gas` attribution.
