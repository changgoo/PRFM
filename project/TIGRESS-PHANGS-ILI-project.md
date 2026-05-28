# TIGRESS-PHANGS-ILI Project

## Overview

This project connects TIGRESS simulation outputs, PHANGS observational constraints, and simulation-based inference/implicit likelihood inference (ILI) to calibrate and test PRFM-related models of star formation and feedback in nearby galaxies.

The central idea is to build an observationally anchored simulation and inference program for the star-forming interstellar medium. PHANGS defines an informed prior over real galactic environments, TIGRESS-NCR provides controlled numerical experiments in those environments, and ILI provides a practical way to infer physical parameters and quantify model uncertainty when explicit likelihoods are unavailable. The PHANGS sampling step is therefore a simulation-design tool: it should preserve the main observed correlations well enough to avoid implausible environmental parameter combinations, while still covering the large input space needed for emulator and inference work.

## Goals

1. Define a PHANGS-informed prior over environmental input parameters for PRFM/TIGRESS comparisons.
2. Build sampling strategies that preserve observed environmental correlations without requiring exact reproduction of the PHANGS sample.
3. Generate a TIGRESS simulation suite over the informed-prior parameter space.
4. Train and validate ILI models that infer physical parameters from observables.
5. Quantify where the PRFM/TIGRESS model family agrees with or fails against PHANGS constraints.

Finally, this will build a simulation suite like CAMELS (https://www.camel-simulations.org/) that enables ILI for the star-forming ISM.

## Motivation

Bring state-of-the-art simulations and high-fidelity observational surveys together to provide a systematic inference framework for ISM physics. This has not yet been done in a way that simultaneously uses resolved galaxy surveys, controlled numerical experiments, and likelihood-free inference for the multiphase star-forming ISM.

The broader motivation follows the spirit of CAMELS: build a large, structured simulation dataset that maps model parameters to observable statistics, then use machine learning and simulation-based inference to calibrate uncertain physics against observations. CAMELS targets cosmology and galaxy formation by spanning cosmological and astrophysical parameter variations, producing theory predictions, marginalizing over nuisance baryonic effects, mapping between simulation classes, and calibrating subgrid parameters to observations. Here, the analogous goal is to span local ISM environmental and physics parameters, generate TIGRESS-based predictions for star formation and feedback observables, and infer which physical parameter combinations are consistent with PHANGS.

The project should also help extract more information from existing and future surveys. Instead of comparing a small number of hand-picked simulations to a few summary plots, the simulation suite should be designed from the start as an inference-ready dataset: broad enough to cover the relevant environmental input volume, structured enough to train emulators, and constrained enough by PHANGS correlations to avoid wasting simulations on unrealistic parameter combinations.

## Data Products

- PHANGS aperture-level tables with derived PRFM inputs.
- PHANGS-informed environmental parameter priors for TIGRESS-NCR simulation design.
- A TIGRESS-NCR simulation suite spanning PHANGS-informed environmental parameters.
- Extended TIGRESS suites spanning uncertain physics parameters.
- TIGRESS simulation summaries mapped onto comparable observables.
- Emulator or ILI training datasets.
- Diagnostic plots comparing observed, sampled, and model-predicted distributions.
- Synthetic multiwavelength observations and derived summary statistics.
- Documentation of priors, parameter transformations, completeness cuts, and validation diagnostics.

## Sub Projects

### 1. Create TIGRESS-NCR simulation suites for environmental parameter variations

This is the first CAMELS-like step: build a controlled simulation suite whose parameter design is motivated by the observed PHANGS distribution rather than by an ad hoc grid. The intent is not to reproduce PHANGS aperture counts exactly, but to use PHANGS as an informed prior for plausible environmental inputs.

- The initial suite targets one or a small number of fixed `Sigma_gas` values.
- Draw plausible simulation-design parameters from PHANGS data, especially `Sigma_star`, `Omega`, and `H_star`.
- Convert PHANGS-derived design fields into TIGRESS-NCR simulation inputs: `Sigma_star`, `z_star`, `Omega`, and `rho_dm`.
- Treat `Sigma_atom`, `Sigma_mol`, and `f_mol` as validation quantities rather than independent simulation inputs for the current TIGRESS-NCR setup.
- Run a suite of simulations with sampled environmental parameters, including `qshear`, while holding physics choices such as `Zgas`, `Zdust`, and feedback/population-synthesis assumptions fixed for the initial suite.
- Assess the resulting `Sigma_SFR`, `Sigma_atom`, `Sigma_mol`, and `f_mol` against PHANGS validation distributions.
- Compare sub-kpc statistical properties: PDFs, two-point statistics, wavelet scattering coefficients, morphology, and eventually field-level statistics.
- Develop emulators for ILI and PRFM parameter inference.
- Predict `Sigma_SFR`, `P_DE`, `sigma_eff`, gas phase balance, feedback yields, and related PRFM quantities.
- Use this suite as concept verification for a larger TIGRESS-PHANGS simulation program.

The design should explicitly track whether the suite is optimized for interpolation within the PHANGS-informed prior, extrapolation to rare environments, or parameter inference. These goals imply different sampling choices and different tolerance for mismatch with the observed marginal distributions.

### 2. Extend the simulation suite for physics parameter variations

After the environmental suite is established, extend the design to uncertain physics parameters. This is the closer analog of CAMELS varying astrophysical subgrid parameters.

- Metallicity variations.
- Dust-to-gas ratio variations.
- Dust model variations.
- Population synthesis model variations.
- IMF variations.
- Feedback strength and feedback channel variations.

The goal is to quantify how much PHANGS observables constrain environmental parameters versus uncertain microphysics and feedback prescriptions.

### 3. Develop multiwavelength synthetic observation pipeline

- Build synthetic emission, dust, and gas tracers from TIGRESS outputs.
- Generate mock JWST dust maps and compare them with observed dust structures.
- Generate mock H-alpha, FUV, IR, HI, and CII observables where possible.
- Apply observational resolution, noise, masking, and aperture definitions consistently.
- Develop field-level inference tests using images or maps.
- Compare the information content of field-level statistics against aperture-level summary statistics.

This pipeline is important because PHANGS and future surveys constrain more than one-point distributions. Spatial structure may be essential for distinguishing feedback models that produce similar global SFRs.

## Core Tasks

### 1. PHANGS Data Preparation

- Load PHANGS megatables for selected aperture definitions.
- Compute derived quantities such as `Sigma_gas`, `Omega`, and `H_star`.
- Define quality cuts for gas, stellar, dynamical, and SFR fields.
- Document completeness and selection effects for each SFR tracer.

### 2. Observed Parameter-Space Sampling

- Select PHANGS pixels in target `Sigma_gas` bands.
- Build informed priors for the simulation-design fields: `Sigma_gas`, `Sigma_star`, `H_star`, `Omega`, and `qshear`.
- Compare observed-pixel LHS sampling with KDE-based and other synthetic sampling backends.
- Preserve observed correlations well enough to avoid implausible simulation inputs.
- Keep `Sigma_atom`, `Sigma_mol`, `f_mol`, and SFR tracers as validation fields unless the simulation setup changes.
- Validate design-field coverage and validation-field behavior separately.

### 3. TIGRESS Model Mapping

- Identify TIGRESS simulation parameters corresponding to PHANGS-derived inputs.
- Build summary statistics for TIGRESS outputs.
- Match simulation outputs to PHANGS observables and units.
- Track model assumptions, priors, and unsupported regions of parameter space.

### 4. ILI / Emulator Development

- Choose inference targets and observable summaries.
- Build training and validation splits.
- Train baseline density-estimation or emulator models.
- Test calibration, posterior coverage, and sensitivity to sampling choices.
- Compare direct emulation of observables with posterior estimation for physical parameters.
- Test whether summary statistics are sufficient for recovering known simulation inputs.
- Quantify posterior degeneracies between environmental parameters and physics parameters.
- Build diagnostics for simulation coverage, out-of-distribution PHANGS pixels, and emulator uncertainty.

### 5. Validation And Diagnostics

- Compare PHANGS-informed priors, sampled simulation inputs, validation fields, and model outputs.
- Check bias by galaxy, radius, gas phase, and SFR tracer.
- Quantify residuals in SFR, pressure, scale height, and feedback-related quantities.
- Identify parameter regions requiring more simulations or better observational cuts.

## Open Questions

- Which SFR tracer should define the primary validation target?
- Should the informed prior be conditioned on galaxy/radius in addition to `Sigma_gas`?
- How should nondetections and missing molecular gas be handled in inference?
- Which TIGRESS summary statistics are most informative for PHANGS observables?
- What level of sampling fairness is scientifically sufficient for an informed prior: 0.1, 0.15, or 0.2 dex?

## Near-Term Milestones

1. Finalize PHANGS sampling class and notebook examples.
2. Produce diagnostic plots for several `Sigma_gas` targets and SFR tracers.
3. Define the first TIGRESS-to-PHANGS observable mapping table.
4. Create a small pilot ILI dataset.
5. Run a first inference sanity check on mock data.

