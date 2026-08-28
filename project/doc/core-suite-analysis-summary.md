# TIGRESS–PHANGS core-suite analysis: current synthesis

**Status:** core environmental suite complete; analysis summarized through 2026-08-28.

## Scope

The core suite contains 32 TIGRESS-NCR models selected around
`Sigma_gas = 10 M_sun pc^-2` (within ±0.3 dex) from a PHANGS-informed
environmental design. At fixed microphysics, the four primary environmental
variables are `Sigma_star`, `H_star`, `Omega`, and `q`; `rho_star`
and `kappa` are derived from them. The finite `Sigma_gas` width is retained
in all multivariate interpretations.

All 32 models are available for the history-based PRFM, observational
validation, phase-resolved z-profile, and projected-field analyses. The
regenerated row0000 sequence provides 601 unique maps over 0–600 Myr and
passes the total-gas and H I projection-quality checks.

## Main conclusions

| Question | Current answer | Confidence |
|---|---|---|
| Does the suite reach the observed `Sigma_SFR` regime? | Yes. The simulated and PHANGS medians agree to about 0.02 dex, and their distributions overlap substantially. | High for coverage; moderate for a distributional match. |
| What drives the simulated `Sigma_SFR` scatter? | Differences between environments dominate temporal fluctuations: 82.6% versus 17.4% of the model variance. | High within the sampled design and 400–600 Myr window. |
| Which environmental axes show the clearest associations? | `Sigma_SFR`, density-PDF width, and phase scale heights primarily track the stellar-gravity axis; power-spectrum scale/anisotropy and phase magnetic speeds primarily track rotation. Phase velocity and effective-support speeds increase most consistently with `Sigma_SFR`. No robust monotonic trend is detected for the fitted spectral slope or with `q` alone. | Moderate: these are bivariate associations in a covariant design. |
| Is vertical dynamical equilibrium recovered? | Yes. The model-mean `P_tot,2p/W` has median 1.14 and range 0.98–1.71. | High. |
| Is the PHANGS molecular fraction recovered? | No. The simulated median `f_mol` is lower by about 1.35 dex. | High for the discrepancy; low for interpreting it as a physical failure. |
| Are column-density power spectra converged with box size? | Their angle-averaged shapes agree well over common resolved wavelengths; larger boxes mainly add longer-wavelength power. | High for the common-scale spectrum. |
| Are scalar power-spectrum summaries converged? | Only conditionally. Fixed-band moments are robust, while unrestricted integral scale, a single fitted slope, and global quadrupole amplitude retain box-size or resolution sensitivity. | High based on the box-size experiment. |
| Do total gas, H I, and emission measure trace the same morphology? | Total gas and H I are tightly coupled across every measured PDF/spectrum diagnostic; EM is largely decoupled and emphasizes compact ionized structures. | High for the complete 32-model suite; H I log-width is low-column sensitive. |
| What does the phase-resolved reduction show? | Neutral gas carries most of the mass, ionized gas most of the volume, and hotter phases are progressively thicker and faster. | High for the descriptive 400–600 Myr summaries; moderate for environmental attribution. |

## PRFM and star-formation response

The suite spans `Sigma_SFR = 3.5e-4` to
`9.1e-3 M_sun kpc^-2 yr^-1`, with median `2.9e-3`. The near-unity
pressure-to-weight ratio shows that the models establish vertical dynamical
equilibrium across the sampled environments. The median total feedback yield
is 1574 km/s; very large yields occur in the lowest-SFR cases and should not
be read as a universal calibration.

The PHANGS comparison is best described as successful **coverage**, not yet a
complete distributional match. The simulated median `Sigma_SFR = 0.002919`
is close to the PHANGS median 0.002801, but the combined suite scatter is
0.415 dex compared with 0.271 dex in PHANGS. The between-environment scatter
is 0.378 dex, whereas the typical temporal scatter within a model is
0.173 dex. Thus environmental variation can readily produce the observed-order
scatter at fixed nominal `Sigma_gas`, and is more important than temporal
variability in this suite. Some low-SFR environments broaden the model
distribution, while the observed high-SFR tail is not fully reached.

This decomposition is descriptive rather than causal. The design variables
are correlated by the PHANGS-informed sampling, and the finite `Sigma_gas`
width can contribute to the apparent parameter dependence. Multivariate
modeling is required before assigning scatter to an individual environmental
variable.

## Molecular fraction

The current models give a positive-sample median `f_mol = 0.0212`, compared
with 0.4798 for PHANGS. This large offset is expected to be strongly affected
by the approximately 8 pc numerical resolution and the inability to resolve
the dense, shielded structures that host molecular gas. The native simulation
diagnostic should therefore be treated as resolution dependent, not as
evidence that the environmental model itself fails.

For completeness, the model variance in `f_mol` divides into 58.4% between
environments and 41.6% temporal variation. That decomposition is useful for
planning a resolution study, but the absolute normalization and the apparent
agreement or disagreement of scatter should not be used for physical
calibration until molecular-gas convergence is established.

## Correlation synthesis

Model-level Spearman rank tests make the principal trends more explicit:

- `Sigma_SFR` increases most clearly with `rho_star`
  (Spearman `rho_s = +0.45`, FDR-adjusted
  `q = 0.043`). Its association with `Sigma_star` is
  similar but marginal after multiple-testing correction
  (`rho_s = +0.40`, `q = 0.071`). No supported monotonic
  relation is found with `H_star`, `Omega`/`kappa`,
  or `q` in the 32-model sample.
- Both column-density PDF widths increase with `rho_star`
  (`rho_s = +0.67` and `+0.72`) and
  `Sigma_SFR` (both about `+0.57`). The linear-density
  width also decreases moderately with `Omega`/`kappa`
  (about `-0.46`).
- The unrestricted power-spectrum integral scale decreases strongly with
  `Omega`/`kappa` (about
  `rho_s = -0.79`) and more moderately with
  `rho_star` (`-0.43`). The fitted slope
  `alpha` has no statistically supported monotonic association with
  the tested environmental variables.
- Global `A_2` increases with `Omega`
  (`+0.66`) and `kappa` (`+0.57`) and
  decreases with `Sigma_SFR` (`-0.63`). This is a
  within-fixed-box association; its absolute normalization remains box-size
  sensitive.
- The resolution-limited `f_mol` diagnostic increases strongly with
  `Omega` (`+0.79`) and `kappa`
  (`+0.68`). This pattern is provisional and should not be given a
  molecular-physics interpretation before a resolution study.

These are not six independent parameter tests. In this design,
`Sigma_star` and `rho_star` have
`rho_s = 0.91`, while `Omega` and `kappa`
have `rho_s = 0.97`. The current evidence therefore supports two
broad response axes—stellar gravity and rotation—more strongly than it
supports attribution to one member of either pair. Full coefficients,
multiple-testing corrections, internal velocity correlations, and
cross-diagnostic trends are recorded in
[the correlation results](core-suite-correlations.md).

## Column-density structure

The face-on density-PDF and shear-aware 2D power-spectrum pipelines now provide
consistent measurements for all 32 projection sets. Because the main
suite shares a box size, pixel scale, Fourier lattice, radial binning, and fit
support, relative trends within this suite can represent real responses to
changes in the TIGRESS-NCR environment. They should be described as suite
associations rather than single-parameter causal effects.

The box-size experiment refines which diagnostics are safe to emphasize:

- The angle-averaged spectrum over common physical wavelengths is the primary
  robust diagnostic.
- Fixed-band variance and the fixed-band characteristic scale are much better
  converged than the unrestricted integral scale. Useful bands are 64–256 pc
  and, where supported, 64–512 pc.
- The fitted slope `alpha` uses the same 64–256 pc support in every box, but
  remains sensitive to spectral curvature and to resolution near the lower
  fit boundary.
- The unrestricted integral scale `L_in` is inherently sensitive to newly
  available long-wavelength modes.
- The global quadrupole amplitude `A_2` includes finite-mode bias and
  cancellation of differently oriented local structures. Local-patch
  quadrupoles, orientation coherence, exact-lattice null tests, and
  matched-mode subsampling are retained as more robust alternatives.

Consequently, parameter trends in the fixed-box TIGRESS-NCR suite remain
scientifically useful even where an absolute diagnostic is not box-size
converged. The strongest claims should be based on common-band spectral shape
or fixed-band moments, with `alpha`, `L_in`, and global `A_2` reported as
conditional diagnostics.

## Projected H I and emission measure

The same face-on PDF and shear-aware spectrum estimators have now been applied
to total gas, H I, and emission measure (EM) for all 32 models. Total gas
and H I are nearly rank-equivalent in linear and logarithmic PDF width
(`rho_s = 0.976` and `0.968`), unrestricted integral scale (`0.907`), and
global quadrupole amplitude (`0.992`); their fitted slopes remain substantially
coupled (`0.704`). EM is much less correlated with either neutral tracer. Its
gas-pair coefficients range from `-0.07` to `0.27`, and the largest EM pairing
is only `rho_s = 0.484` for the H I--EM slope. This supports treating EM as a
distinct compact-ionized-gas morphology rather than a substitute gas-column
tracer.

The H I logarithmic width requires a specific caveat. In a full-VTK check of a
prominent low-column bump, removing cells with `x_HI < 0.01` eliminated the
selected sightlines' tiny H I columns while removing only 0.097% of the global
H I column. The linear width was unchanged, but the positive-column log width
fell from 3.67 to 2.45. Production measurements remain uncensored for a common
definition; future observational comparisons should report censored area and
use an explicit H I column threshold or a robust conditional width.

## Phase-resolved structure and dynamics

The completed z-profile reduction tracks CNM+CMM, UNM, WNM, WIM, WHIM, and HIM
without requiring full VTK snapshots. UIM is preserved as an explicit closure
residual, and neutral, ionized, and true whole-gas aggregates are reconstructed
before temporal statistics. The ensemble medians below are medians across the
32 models of each model's 400–600 Myr temporal median:

| Phase | Mass fraction | Volume fraction | Mass height (pc) | `sigma_3D` (km/s) | `sigma_eff,z` (km/s) |
|---|---:|---:|---:|---:|---:|
| CNM+CMM | 0.175 | 0.0016 | 77 | 8.1 | 5.8 |
| UNM | 0.248 | 0.028 | 138 | 11.0 | 9.8 |
| WNM | 0.408 | 0.147 | 249 | 14.8 | 15.0 |
| WIM | 0.112 | 0.131 | 396 | 23.5 | 21.1 |
| WHIM | 0.0045 | 0.164 | 811 | 61.9 | 60.7 |
| HIM | 0.0018 | 0.482 | 1127 | 135.3 | 182.9 |

The reduced neutral component has median mass fraction 0.867 but volume
fraction 0.199; the ionized component has mass fraction 0.118 but volume
fraction 0.784. Phase ordering therefore cleanly separates the mass-dominant
neutral disk from the volume-filling ionized atmosphere.

The strongest phase associations reinforce the two broad environmental axes.
Neutral mass and volume scale heights decrease with `Sigma_star`
(`rho_s = -0.927` and `-0.933`) and `rho_star` (`-0.845` and `-0.869`).
Neutral effective vertical support increases with mean `Sigma_SFR`
(`+0.831`), while whole-gas support correlates with both mean `Sigma_SFR`
(`+0.744`) and stellar gravity (`+0.715` with `Sigma_star`, `+0.740` with
`rho_star`). Mean-field Alfvén speed tracks rotation across phases; for the
whole gas its coefficient with `Omega` is `+0.850`, while the ionized
perturbed-field speed reaches `+0.887`. These unadjusted rank coefficients are
exploratory, not independent causal tests.

## Overall assessment

The core suite is already sufficient to support five conclusions:

1. TIGRESS-NCR recovers the central PHANGS `Sigma_SFR` regime while
   maintaining vertical dynamical equilibrium.
2. Across the sampled environments, variations between models explain more of
   the simulated `Sigma_SFR` scatter than temporal fluctuations within a
   model.
3. Column-density structure can be compared systematically within the
   fixed-box suite, provided that common-band diagnostics are prioritized.
4. Total-gas and H I fractional morphology are tightly coupled, whereas EM
   provides a substantially different view of compact ionized structure.
5. The phase-resolved products recover an ordered neutral-disk/ionized-halo
   structure and expose distinct stellar-gravity, rotation, and SFR response
   axes in phase thicknesses and characteristic speeds.

It does **not** yet support a molecular-fraction calibration, a causal
one-variable attribution of the observed scatter, or universal absolute
values for box-sensitive power-spectrum summaries. It also does not make the
native H I low-column tail directly observation-ready or turn the phase rank
correlations into controlled one-parameter experiments.

## Priorities

1. Fit a multivariate emulator or hierarchical variance model using
   `Sigma_gas`, the four primary variables, and the two derived variables.
   Use conditional or partial-dependence diagnostics to distinguish individual
   environmental responses from design covariance.
2. Carry the full temporal distributions into posterior-predictive comparisons
   rather than reducing each model to a median and percentile range.
3. Include the tracer and phase summaries in the multivariate analysis, using
   explicit censoring/selection models for H I and observation-matched
   synthetic treatment for EM.
4. Promote common-band power-spectrum shape, fixed-band variance, and
   fixed-band characteristic scale to the primary density diagnostics.
5. Test molecular-gas convergence with higher spatial resolution and, if
   needed, a synthetic observation/selection treatment before drawing
   conclusions from `f_mol`.
6. Complete the planned prior augmentation in `H_star` and `Omega`, then add
   controlled physics-parameter variations to separate environmental from
   microphysical effects.

## Detailed records

The purpose of this document is synthesis. Definitions, implementation
choices, reproduction commands, and figure-level notes remain in the focused
records below.

- Suite design: [sampling design](../TIGRESS-PHANGS-sampling-design.md) and
  [Sobol sampling workflow](../TIGRESS-PHANGS-sobol-sampling.md)
- Observational comparison and variance decomposition:
  [validation diagnostics](validation-diagnostics.md)
- Model-level environmental and cross-diagnostic associations:
  [correlation results](core-suite-correlations.md)
- Proposed design expansion: [prior augmentation](prior-augmentation.md)
- Proposed physics experiments:
  [physics-parameter extension](physics-parameter-extension.md)
- Project workflow and provenance: [workflow](../WORKFLOW.md)
- PRFM definitions and conventions:
  [tigress_ncr_tools PRFM analysis](https://github.com/changgoo/tigress_ncr_tools/blob/analysis/ncr-suite/docs/prfm_analysis.md)
- Density-PDF analysis:
  [tigress_ncr_tools density PDF](https://github.com/changgoo/tigress_ncr_tools/blob/analysis/ncr-suite/docs/density_pdf.md)
- 2D density power spectra:
  [tigress_ncr_tools power-spectrum analysis](https://github.com/changgoo/tigress_ncr_tools/blob/analysis/ncr-suite/docs/density_power_spectrum.md)
- Box-size convergence and diagnostic alternatives:
  [tigress_ncr_tools box-size study](https://github.com/changgoo/tigress_ncr_tools/blob/analysis/ncr-suite/docs/box_size_density_power_spectrum.md)
- H I/EM estimators, low-column test, and cross-tracer results:
  [tigress_ncr_tools projected-species analysis](https://github.com/changgoo/tigress_ncr_tools/blob/analysis/ncr-suite/docs/projected_species_analysis.md)
- Phase definitions, estimators, correlations, and products:
  [tigress_ncr_tools phase-resolved analysis](https://github.com/changgoo/tigress_ncr_tools/blob/analysis/ncr-suite/docs/phase_resolved_analysis.md)
