# PHANGS-Informed Sampling Design For TIGRESS-PHANGS-ILI

## Purpose

The sampling step defines which TIGRESS-NCR simulations should be run first. Its role is to turn the PHANGS aperture-level observations into a compact, physically valid, and inference-useful set of local ISM environments.

The design is motivated by the same logic as CAMELS: construct a structured simulation suite that spans physically relevant parameters, supports emulator and ILI training, and can be compared directly with observations. For this project, the target is not cosmological parameter inference, but local star-forming ISM inference: given observed PHANGS environments and observables, infer which TIGRESS/PRFM physical parameters reproduce the observed gas, star formation, pressure, and feedback statistics.

## Statistical Target

Let each PHANGS aperture be represented by a vector

```text
x = (Sigma_gas, Sigma_atom, Sigma_mol, Sigma_star, H_star, Omega_d, y)
```

where `y` denotes one or more SFR tracers, for example

```text
y in {Sigma_SFR_HaW4recal, Sigma_SFR_FUVW4recal, Sigma_SFR_Hacorr}.
```

The first sampling objective is conditional on total gas surface density:

```text
x ~ p(x | log10 Sigma_gas in [log10 Sigma_gas,target - Delta,
                              log10 Sigma_gas,target + Delta]).
```

This isolates local environments at approximately fixed gas surface density, then samples the remaining environmental degrees of freedom that enter TIGRESS-NCR initial conditions.

The main simulation-design environmental parameters are

```text
theta_design = (Sigma_gas, Sigma_star, H_star, Omega_d).
```

Gas phase information is retained for validation through

```text
Sigma_gas = Sigma_atom + Sigma_mol,
f_mol = Sigma_mol / Sigma_gas.
```

For the current TIGRESS-NCR setup, `f_mol`, `Sigma_atom`, and `Sigma_mol` should not be treated as independent simulation input parameters. They are validation targets for the simulated phase balance.

The TIGRESS-NCR input mapping is approximately

```text
(Sigma_gas, Sigma_star, H_star, Omega_d) ->
(Sigma, Sigma_star, z_star, Omega, rho_dm),
```

where the precise mapping must be documented in the simulation setup layer.

## Design Fields And Validation Fields

The sampling workflow should distinguish between fields that define simulation inputs and fields used to validate simulation outputs.

Simulation-design fields:

```text
D = (Sigma_gas, Sigma_star, H_star, Omega_d).
```

Validation fields:

```text
V = (Sigma_atom, Sigma_mol, f_mol, Sigma_SFR, P_DE, sigma_eff, feedback yields, ...).
```

The design fields determine which TIGRESS-NCR simulations are run. The validation fields determine whether those simulations reproduce PHANGS gas partitioning, star formation, and PRFM-related observables.

## Physical Constraints For Validation Quantities

Gas components must obey

```text
Sigma_atom >= 0,
Sigma_mol  >= 0,
Sigma_atom <= Sigma_gas,
Sigma_mol  <= Sigma_gas,
Sigma_atom + Sigma_mol = Sigma_gas.
```

A naive KDE over `(Sigma_gas, Sigma_atom, Sigma_mol)` can violate these constraints because the three quantities are sampled as correlated but independent coordinates. The current implementation avoids this by transforming gas validation variables to

```text
Sigma_gas,
f_mol = Sigma_mol / Sigma_gas,
```

with

```text
0 <= f_mol <= 1.
```

For KDE fitting, the molecular fraction is represented by a logit-like transform:

```text
u = log10(f_mol / (1 - f_mol)).
```

The inverse transform is

```text
f_mol = 10^u / (1 + 10^u),
Sigma_mol  = f_mol Sigma_gas,
Sigma_atom = Sigma_gas - Sigma_mol.
```

Small floors are used only for numerical stability when `f_mol` or `Sigma_mol` is zero. Floor values should not be interpreted as physical detections and are excluded from visualization when appropriate.

## Selection Of Reference Pixels

For a target gas surface density `Sigma_gas,target` and a half-width `Delta` in dex, define the mask

```text
M_gas = |log10 Sigma_gas - log10 Sigma_gas,target| <= Delta.
```

Additional validity requirements are imposed on the fields needed for sampling:

```text
Sigma_gas  > 0,
Sigma_atom > 0,
Sigma_star > 0,
H_star     > 0,
Omega_d    > 0.
```

The selected reference set is

```text
D_ref = {x_i : M_gas(i) and validity cuts are satisfied}.
```

This reference set is the empirical PHANGS distribution against which samples are judged.

## Observed-Pixel Latin Hypercube Sampling

The most conservative sampling method selects actual PHANGS pixels. It preserves observed correlations and attached quantities such as SFR, radius, galaxy identity, and any unmodeled observational systematics.

For selected fields

```text
theta = (Sigma_star, H_star, Omega_d),
```

first transform each field to log space:

```text
z_j = log10 theta_j.
```

Then convert each dimension to empirical rank coordinates:

```text
r_ij = rank(z_ij) / N,
```

where `N` is the number of valid reference pixels. This maps the observed cloud into the unit cube `[0,1]^d` while preserving marginal ranks.

A Latin hypercube design draws `n` target points

```text
u_k in [0,1]^d,  k = 1,...,n,
```

with stratification in each dimension. Each LHS point is matched to the nearest available observed point in rank space. The result is an observed-pixel sample

```text
S_obs subset D_ref.
```

Advantages:

- physically observed combinations only;
- real SFR values are retained;
- no generative-model assumptions;
- useful as a baseline.

Limitations:

- cannot create new environments between observed points;
- may require large `n` for fair coverage;
- rare regions and incomplete SFR tracers can dominate the fairness metric.

## Synthetic KDE-LHS Sampling

The synthetic method fits an approximate continuous distribution to the reference pixels, then samples from it. The current backend is a log-space Gaussian KDE with LHS selection from an oversampled KDE candidate pool.

Let the fitted vector be

```text
w = (log10 Sigma_gas,
     logit10 f_mol,
     log10 Sigma_star,
     log10 H_star,
     log10 Omega_d,
     log10 y)
```

where `y` is an optional SFR field used for joint validation-oriented sampling. For simulation-design sampling, the LHS coordinates should normally be restricted to the design fields rather than `f_mol` or SFR.

The KDE model is

```text
p_hat(w) = (1/N) sum_i K_H(w - w_i),
```

where `K_H` is a multivariate Gaussian kernel with bandwidth matrix `H` selected by the KDE implementation.

The algorithm is:

1. Fit `p_hat(w)` to complete-case reference pixels.
2. Draw a large candidate pool from the KDE.
3. Convert selected LHS fields to empirical rank coordinates.
4. Select candidates nearest to an LHS design.
5. Transform back to physical variables.
6. Reconstruct `Sigma_atom` and `Sigma_mol` from `Sigma_gas` and `f_mol` for validation plots and diagnostics.

Advantages:

- can generate compact samples much smaller than the observed-pixel sample;
- supports interpolation between observed environments;
- naturally produces emulator/ILI design points;
- can include SFR or gas phase quantities in the joint fitted distribution for validation diagnostics.

Limitations:

- KDE can smooth or bias sparse distributions;
- high-dimensional KDE is data-hungry;
- complete-case filtering can change the effective target distribution;
- if only one SFR tracer is fitted, that sampled SFR should not be assumed to represent all SFR tracers.


## Proposed Synthetic Backends

The current KDE-LHS implementation is a useful baseline, but the sampling class should support multiple synthetic backends. The goal is not to choose one method permanently, but to compare how backend assumptions affect the resulting TIGRESS simulation design.

### Backend 1: Constrained Gaussian KDE

This is the current implementation. It fits a Gaussian KDE in transformed variables:

```text
(log10 Sigma_gas,
 logit10 f_mol,
 log10 Sigma_star,
 log10 H_star,
 log10 Omega_d,
 log10 y).
```

It is simple, non-parametric, and easy to inspect. It is useful for first-pass interpolation inside the PHANGS distribution and for debugging the sampling workflow.

Expected strengths:

- quick to implement and run;
- no need to choose the number of mixture components;
- gives a smooth approximation to the observed distribution;
- works naturally with the existing LHS candidate-selection step.

Expected weaknesses:

- high-dimensional KDE is noisy when the complete-case sample is small;
- Gaussian kernels can oversmooth sharp edges and tails;
- bandwidth selection can strongly affect SFR distributions;
- complete-case filtering can make the fitted distribution differ from the marginal observed distribution.

This backend should remain the baseline because it is transparent and provides a reference for more structured methods.

### Backend 2: Gaussian Mixture Model With Physical Transforms

A Gaussian mixture model (GMM) would fit the transformed variable vector as

```text
p(w) = sum_k pi_k N(w | mu_k, Sigma_k).
```

The same physical transforms should be used as in the KDE backend: log variables for positive quantities and `logit10 f_mol` for the gas partition. The number of mixture components can be selected with BIC, cross-validation, or posterior predictive diagnostics.

Expected strengths:

- can represent multimodal structure better than a single global KDE bandwidth;
- can generate large candidate pools cheaply after fitting;
- gives interpretable mixture components that may correspond to radial, galaxy, or phase regimes;
- can regularize covariance matrices when the reference sample is small.

Expected weaknesses:

- requires choosing or validating the number of components;
- Gaussian components can still generate unphysical tails unless transformed variables and clipping are handled carefully;
- small SFR-complete samples may lead to unstable mixture fits;
- component interpretation can be misleading if completeness effects dominate.

This backend is attractive if PHANGS environments separate into recognizable clusters, for example atom-dominated and molecule-dominated regimes or inner- and outer-disk regimes.

### Backend 3: Copula Or Conditional Model

A copula-based backend separates marginal distributions from dependence structure. One possible construction is:

```text
u_j = F_j(w_j),
```

where each `F_j` is an empirical or KDE-smoothed marginal CDF, followed by a model for the joint distribution of `u` in the unit cube. This could use a Gaussian copula, vine copula, or another dependence model.

A related option is a conditional SFR model:

```text
p(y | Sigma_gas, f_mol, Sigma_star, H_star, Omega_d),
```

combined with a separate sampler for the environmental variables. This is useful if the simulation design should sample environments first and treat SFR as an observed validation target rather than as a design coordinate.

Expected strengths:

- preserves observed one-dimensional marginals more directly than KDE or GMM;
- can compare different assumptions about dependence separately from marginal completeness;
- conditional modeling makes the role of SFR explicit;
- can support multiple SFR tracers as alternative conditionals.

Expected weaknesses:

- more implementation complexity;
- empirical CDFs need careful handling of missing values and nondetections;
- copula assumptions may fail in sparse tails;
- conditional SFR models require separate validation and may hide selection effects.

This backend is likely the best long-term direction if the KDE and GMM backends show persistent marginal bias, especially in SFR tracers.

### Backend 4: Normalizing Flow

A normalizing flow would learn an invertible transformation between a simple latent distribution and the transformed PHANGS variables:

```text
z ~ N(0, I),
w = f_phi(z),
```

where `w` is the physically transformed vector, for example

```text
w = (log10 Sigma_gas,
     logit10 f_mol,
     log10 Sigma_star,
     log10 H_star,
     log10 Omega_d,
     log10 y).
```

The learned density is

```text
p(w) = p(z) |det d f_phi^{-1} / d w|.
```

Like the KDE and GMM backends, the flow should operate in transformed coordinates and then reconstruct physical gas components from `Sigma_gas` and `f_mol`. This keeps the generated samples physically valid while allowing a more flexible distribution than KDE or GMM.

Expected strengths:

- can represent nonlinear, non-Gaussian, and multimodal distributions;
- can scale better than KDE in moderately high dimensions;
- can generate large candidate pools efficiently after training;
- can support conditional variants, for example `p(environment | Sigma_gas)` or `p(SFR | environment)`;
- provides an explicit density that can be useful for ILI diagnostics and out-of-distribution checks.

Expected weaknesses:

- needs substantially more training and validation machinery;
- can overfit small complete-case samples, especially for SFR fields;
- generated marginals can look plausible while tails or conditional structure are wrong;
- requires choices about architecture, regularization, training epochs, and validation metrics;
- less transparent than KDE or GMM for early-stage debugging.

A normalizing flow is most attractive after the simpler backends have clarified the failure modes. It may be especially useful for a larger simulation-design stage with many `Sigma_gas` slices, additional environmental variables, or field-level summary statistics. For the current PHANGS sampling problem, a conditional flow may be more useful than an unconditional flow:

```text
p(Sigma_star, H_star, Omega_d, f_mol, y | Sigma_gas)
```

or

```text
p(y | Sigma_gas, f_mol, Sigma_star, H_star, Omega_d).
```

The first form supports validation-aware simulation-design sampling; the second form directly targets the SFR-bias problem. In both cases, the actual TIGRESS design coordinates should remain the simulation input fields unless the simulation setup changes.

### Backend Comparison Criteria

Each backend should be evaluated with the same diagnostics:

- design-field quantile fairness for `Sigma_gas`, `Sigma_star`, `H_star`, and `Omega_d`;
- physical gas constraints are exactly satisfied for validation quantities;
- validation-field fairness for `f_mol`, `Sigma_atom`, `Sigma_mol`, and SFR tracers;
- SFR quantile fairness for one or multiple SFR tracers;
- comparison of observed, fitted-model, and finite-sample distributions;
- held-out likelihood or posterior-predictive checks where the backend supports density evaluation;
- sensitivity to `Sigma_gas,target`, `Delta`, and completeness cuts;
- stability under different random seeds;
- ability to generate a compact but representative TIGRESS simulation design.

## SFR Fields And Complete-Case Bias

SFR tracers have different validity masks. For example, the reference set may contain many valid `Sigma_SFR_FUVW4recal` values but far fewer valid `Sigma_SFR_Hacorr` values. A joint KDE fitted with all SFR fields targets

```text
p(x | all selected SFR fields are valid),
```

not

```text
p(x | each field is valid independently).
```

These distributions can differ substantially.

To avoid an overly restrictive complete-case cut, the current workflow allows one SFR field to be used for joint KDE sampling:

```text
sfr_fields = "Sigma_SFR_Hacorr"
```

while comparing the sampled SFR distribution against multiple observed SFR tracers:

```text
comparison_sfr_fields = [
    "Sigma_SFR_HaW4recal",
    "Sigma_SFR_FUVW4recal",
    "Sigma_SFR_Hacorr",
]
```

This makes the interpretation explicit: the sample contains one synthetic SFR variable, and the plots ask whether that synthetic distribution resembles each observed tracer.

## Fairness Metric

Sampling fairness is currently evaluated using quantile mismatch in log space. For field `a`, define reference and sample log-values:

```text
A_ref = log10 a_ref,
A_sam = log10 a_sam.
```

For quantiles

```text
q in {0.1, 0.2, ..., 0.9},
```

compute

```text
epsilon_a = max_q |Q_q(A_sam) - Q_q(A_ref)|.
```

For a set of fields `F`, the total mismatch is

```text
epsilon_F = max_{a in F} epsilon_a.
```

A sample is considered fair if

```text
epsilon_F <= epsilon_tol.
```

Typical tolerances explored so far are

```text
0.10 dex, 0.15 dex, 0.20 dex.
```

A tolerance of `0.1 dex` is strict: it corresponds to roughly 26 percent in linear space. For sparse SFR-complete samples, this can force sample sizes close to the full observed set. A tolerance of `0.15-0.20 dex` may be more practical for simulation design.

## Direct SFR Assignment For Primary Samples

If a sample is generated only in primary environmental space, it does not have an SFR value unless one is assigned afterward. The current helper assigns SFR using distance-weighted nearest neighbors in gas + primary space.

For a synthetic point `s`, find nearby observed points using

```text
(log10 Sigma_gas,
 logit10 f_mol,
 log10 Sigma_star,
 log10 H_star,
 log10 Omega_d).
```

Then draw an SFR value from the nearest observed neighbors with weights approximately

```text
w_i proportional to 1 / d_i^2.
```

This is useful for diagnostics but should not be treated as a physical SFR model. If SFR is a target observable, the preferred route is joint sampling or a dedicated conditional SFR model.

## Diagnostic Plots

The sampling notebook uses three classes of diagnostic plots:

1. Correlation matrix with selected `Sigma_gas` bands overplotted.
2. Observed-vs-sample marginal distribution overlays.
3. Observed-vs-sample pair distribution overlays.

For KDE samples, the marginal overlay can also show the KDE model distribution itself. The visual comparison then separates three objects:

- observed reference distribution;
- actual finite sample selected from the design;
- fitted KDE model distribution.

This is important when the sample appears biased. The bias may come from the KDE fit, from the LHS selection step, from complete-case filtering, or from comparing against a different observed validity mask.

## Recommended Workflow

1. Choose `Sigma_gas,target` and `Delta`.
2. Inspect the selected PHANGS pixels in the correlation matrix.
3. Decide which SFR tracer, if any, defines the joint sampling target.
4. Generate an observed-pixel LHS baseline.
5. Generate a KDE-LHS synthetic sample.
6. Verify design-field fairness in `Sigma_gas`, `Sigma_star`, `H_star`, and `Omega_d`.
7. Verify gas phase and SFR validation distributions, paying attention to completeness masks.
8. Decide whether the sample is intended for interpolation, extrapolation, or inference.
9. Convert only the design fields into TIGRESS-NCR simulation parameters.
10. Track validation fields separately from simulation input fields in the design table.

## Open Methodological Questions

- Should the primary sampling space include radius or galaxy-level stratification?
- Should SFR be part of the design space or only a validation observable?
- Should KDE be replaced by a Gaussian mixture, normalizing flow, or copula model?
- What tolerance is scientifically sufficient for simulation design?
- How should nondetections be represented in a likelihood-free inference pipeline?
- Should the first simulation suite prioritize broad PHANGS coverage or dense local coverage around a few representative environments?
