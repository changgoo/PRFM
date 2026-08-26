# TIGRESS–PHANGS core-suite analysis: current synthesis

**Status:** core environmental suite complete; analysis summarized through 2026-08-26.

## Scope

The core suite contains 32 TIGRESS-NCR models selected around
`Sigma_gas = 10 M_sun pc^-2` (within ±0.3 dex) from a PHANGS-informed
environmental design. At fixed microphysics, the four primary environmental
variables are `Sigma_star`, `H_star`, `Omega`, and `q`; `rho_star`
and `kappa` are derived from them. The finite `Sigma_gas` width is retained
in all multivariate interpretations.

All 32 models are available for the history-based PRFM and observational
validation analyses. Column-density statistics use 31 models because the
row-0000 projection sequence is corrupted.

## Main conclusions

| Question | Current answer | Confidence |
|---|---|---|
| Does the suite reach the observed `Sigma_SFR` regime? | Yes. The simulated and PHANGS medians agree to about 0.02 dex, and their distributions overlap substantially. | High for coverage; moderate for a distributional match. |
| What drives the simulated `Sigma_SFR` scatter? | Differences between environments dominate temporal fluctuations: 82.6% versus 17.4% of the model variance. | High within the sampled design and 400–600 Myr window. |
| Which environmental axes show the clearest associations? | `Sigma_SFR` and density-PDF width primarily track the stellar-gravity axis; power-spectrum integral scale and anisotropy primarily track the rotation axis. No robust monotonic trend is detected for the fitted spectral slope or with `q` alone. | Moderate: these are bivariate associations in a covariant design. |
| Is vertical dynamical equilibrium recovered? | Yes. The model-mean `P_tot,2p/W` has median 1.14 and range 0.98–1.71. | High. |
| Is the PHANGS molecular fraction recovered? | No. The simulated median `f_mol` is lower by about 1.35 dex. | High for the discrepancy; low for interpreting it as a physical failure. |
| Are column-density power spectra converged with box size? | Their angle-averaged shapes agree well over common resolved wavelengths; larger boxes mainly add longer-wavelength power. | High for the common-scale spectrum. |
| Are scalar power-spectrum summaries converged? | Only conditionally. Fixed-band moments are robust, while unrestricted integral scale, a single fitted slope, and global quadrupole amplitude retain box-size or resolution sensitivity. | High based on the box-size experiment. |

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
0.414 dex compared with 0.271 dex in PHANGS. The between-environment scatter
is 0.376 dex, whereas the typical temporal scatter within a model is
0.172 dex. Thus environmental variation can readily produce the observed-order
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
  (Spearman `rho_s = +0.44`, FDR-adjusted
  `q = 0.047`). Its association with `Sigma_star` is
  similar but marginal after multiple-testing correction
  (`rho_s = +0.39`, `q = 0.071`). No supported monotonic
  relation is found with `H_star`, `Omega`/`kappa`,
  or `q` in the 32-model sample.
- Both column-density PDF widths increase with `rho_star`
  (`rho_s = +0.67` and `+0.72`) and
  `Sigma_SFR` (both about `+0.55`). The linear-density
  width also decreases moderately with `Omega`/`kappa`
  (about `-0.42`).
- The unrestricted power-spectrum integral scale decreases strongly with
  `Omega`/`kappa` (about
  `rho_s = -0.82`) and more moderately with
  `rho_star` (`-0.44`). The fitted slope
  `alpha` has no statistically supported monotonic association with
  the tested environmental variables.
- Global `A_2` increases with `Omega`
  (`+0.67`) and `kappa` (`+0.58`) and
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
consistent measurements for the 31 clean projection sets. Because the main
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

## Overall assessment

The core suite is already sufficient to support three conclusions:

1. TIGRESS-NCR recovers the central PHANGS `Sigma_SFR` regime while
   maintaining vertical dynamical equilibrium.
2. Across the sampled environments, variations between models explain more of
   the simulated `Sigma_SFR` scatter than temporal fluctuations within a
   model.
3. Column-density structure can be compared systematically within the
   fixed-box suite, provided that common-band diagnostics are prioritized.

It does **not** yet support a molecular-fraction calibration, a causal
one-variable attribution of the observed scatter, or universal absolute
values for box-sensitive power-spectrum summaries.

## Priorities

1. Fit a multivariate emulator or hierarchical variance model using
   `Sigma_gas`, the four primary variables, and the two derived variables.
   Use conditional or partial-dependence diagnostics to distinguish individual
   environmental responses from design covariance.
2. Carry the full temporal distributions into posterior-predictive comparisons
   rather than reducing each model to a median and percentile range.
3. Promote common-band power-spectrum shape, fixed-band variance, and
   fixed-band characteristic scale to the primary density diagnostics.
4. Test molecular-gas convergence with higher spatial resolution and, if
   needed, a synthetic observation/selection treatment before drawing
   conclusions from `f_mol`.
5. Complete the planned prior augmentation in `H_star` and `Omega`, then add
   controlled physics-parameter variations to separate environmental from
   microphysical effects.
6. Repair or regenerate the row-0000 projection products if a complete
   32-model column-density sample becomes important.

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
