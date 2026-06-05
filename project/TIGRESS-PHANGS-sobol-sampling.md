# Sobol-Based Sampling Design For TIGRESS-PHANGS-ILI

## Motivation: Why Sobol Over Latin Hypercube

The Latin Hypercube Sampling (LHS) method used in the initial implementation has a fundamental limitation for an iteratively grown simulation suite: it provides no extensibility guarantee. To add simulations to an existing LHS design while maintaining its space-filling properties, the entire design must be resampled from scratch. This is incompatible with the CAMELS-style philosophy of the TIGRESS-PHANGS program, where the simulation suite is expected to grow progressively — starting from a small pilot set and expanding as computational resources allow or as new scientific questions arise.

Sobol sequences are quasi-random low-discrepancy sequences that solve this problem directly. A Sobol sequence of length $2^{k+1}$ contains the $2^k$-point sequence as its first half. Adding simulations never displaces or invalidates earlier design points. The resulting suite is therefore coherent across all stages of the project, and emulators or inference models trained at each stage remain compatible with the accumulated dataset.

## Environmental Parameter Space For The First Suite

The first simulation suite samples the environmental input parameters for TIGRESS-NCR. The design vector is restricted to the five fields that enter directly as simulation inputs:

```text
theta_env = (Sigma_gas, Sigma_star, H_star, Omega, qshear)
```

The molecular fraction `f_mol` and SFR tracers (`Sigma_SFR`) are **excluded** from the sampling input for two reasons:

1. They are not independent TIGRESS-NCR inputs. `f_mol`, `Sigma_atom`, and `Sigma_mol` are emergent outputs of the simulation, determined by the gas thermochemistry given `Sigma_gas`. `Sigma_SFR` is similarly an emergent quantity.
2. Including them in the KDE fit would require rows to have valid SFR and molecular gas detections, creating a complete-case bias that biases the reference set toward molecular-rich, actively star-forming environments and reduces the number of usable reference pixels.

Both `f_mol` and `Sigma_SFR` are used only for **validation** after TIGRESS-NCR simulations are run: the emergent values are compared to the observed PHANGS distributions as the primary test of the simulation suite.

## Physical Constraints On Design Fields

All five design fields must be positive. In addition, the shear parameter satisfies a physical bound:

```text
Sigma_gas  > 0
Sigma_star > 0
H_star     > 0
Omega      > 0
0 < qshear <= 1.5
```

The shear parameter `q = -d ln Omega / d ln R` has a natural range set by galactic rotation curves:

- `q = 0`: solid-body rotation (innermost regions, bars)
- `q = 1`: flat rotation curve (typical disk)
- `q = 1.5`: Keplerian rotation (outer disk or point-mass dominated)

Values above 1.5 are super-Keplerian and unphysical for a rotationally supported galactic disk. Reference pixels with `q > 1.5` are excluded from the reference set `D_ref` before KDE fitting. Generated samples that map to `q > 1.5` through the KDE quantile function are rejected and replaced by the next sequential Sobol point.

## KDE Prior On Design Fields

Let the reference set `D_ref` be the PHANGS apertures satisfying the Sigma_gas band mask and all validity cuts. Define the log-transformed design vector:

```text
w = (log10 Sigma_gas,
     log10 Sigma_star,
     log10 H_star,
     log10 Omega,
     log10 qshear)
```

This is a 5-dimensional vector. No logit transform is needed since `f_mol` is excluded. Fit a multivariate Gaussian KDE:

```text
p_hat(w) = (1/N) sum_{i=1}^{N} K_H(w - w_i)
```

where `N = |D_ref|` and `K_H` is a Gaussian kernel with bandwidth matrix `H` selected by the KDE implementation (e.g., Scott's rule or cross-validation). The KDE is fit entirely in the 5D design space, requiring only that the five design fields are valid and positive.

From the fitted KDE, extract marginal CDFs `F_j(w_j)` and their inverses (quantile functions) `Q_j(u) = F_j^{-1}(u)` for each dimension `j = 1, ..., 5`. In practice these are evaluated numerically from a large sample of the fitted KDE.

## KDE-Sobol Sampling Algorithm

Given a target sample size `n`:

1. Generate a scrambled Sobol sequence of length `n'` >= `n` in `[0,1]^5`:
   ```python
   sampler = scipy.stats.qmc.Sobol(d=5, scramble=True, seed=42)
   u = sampler.random_base2(m)  # n = 2^m; or sampler.random(n) for arbitrary n
   ```
   Using `random_base2(m)` with `n = 2^m` guarantees optimal uniformity properties.

2. Map each Sobol point `u_k = (u_{k,1}, ..., u_{k,5})` through the KDE marginal quantile functions:
   ```text
   w_{k,j} = Q_j(u_{k,j})  for j = 1, ..., 5
   ```
   This is the probability integral transform (inverse CDF mapping), which converts uniform Sobol points into samples from the KDE-estimated prior.

3. Convert from log space back to physical variables:
   ```text
   theta_{k,j} = 10^{w_{k,j}}  for each dimension j
   ```

4. Apply physical bounds: reject any point with `qshear > 1.5` or any negative physical value, and replace with the next sequential Sobol point. Track the number of replacements; if the rejection fraction exceeds ~10%, investigate whether the KDE is generating excessive unphysical tails.

5. The resulting sample `S_Sobol = {theta_k : k = 1, ..., n}` is the simulation design.

Note: unlike the KDE-LHS method, no second nearest-neighbor selection step is needed. The Sobol-through-quantile approach maps directly to physical space without requiring a separate candidate pool.

## Progressive Extension Property

Let `S(n)` denote the Sobol-based sample of size `n`. The key property is:

```text
S(2^k) subset S(2^{k+1})  for all k >= 0
```

That is, the first `2^k` points of the sequence are exactly the `2^k`-point design. Adding `2^k` new simulations to reach `2^{k+1}` requires running only the new points; the existing simulations are not displaced. This enables:

- Starting with a pilot suite of `n = 64` simulations for initial PRFM validation
- Extending to `n = 128` by running 64 additional simulations that fill coverage gaps
- Further extension to `n = 256, 512, ...` following the same pattern
- Merging datasets across stages for emulator training without re-weighting

Canonical suite sizes are powers of two: `n in {64, 128, 256, 512, 1024}`. Arbitrary `n` is permitted but produces slightly less uniform coverage at non-power-of-two sizes. The Sobol sequence optimality guarantee applies at `n = 2^m` for integer `m`.

The rejection of unphysical `q > 1.5` points uses sequential replacement: if point `k` is rejected, point `k + n_extra` is used instead, where `n_extra` is the count of extra points generated beyond `n`. As long as the rejection fraction is small (expected to be <1% given the observed PHANGS distribution of `q`), the resulting sample still has near-optimal uniformity at each canonical size.

## Block-Structured Extension To Physics Parameters

The first suite fixes physics parameters (metallicity `Z'`, dust-to-gas ratio, IMF, feedback parameterization) at fiducial values and varies only environmental inputs `theta_env`. When the program extends to physics parameter variations, a separate Sobol block is used:

```text
Block 1 (environmental): Sobol(d=5) over theta_env
Block 2 (physics):       Sobol(d=d_phys) over theta_phys
```

where `theta_phys` might include `Z'`, dust-to-gas ratio normalization, IMF upper-mass cutoff, feedback yield scaling, etc.

The combined simulation design pairs environmental and physics points:

```text
(theta_env_i, theta_phys_j)  for selected index pairs (i, j)
```

Each block can be extended independently. For example:
- Extend the environmental block from 64 to 128 while holding physics fixed
- For a fixed environmental point, run a Sobol sweep over physics at `d_phys` dimensions

This mirrors the CAMELS approach: CAMELS varies cosmological and astrophysical parameters in separate Latin Hypercube suites. Here, environmental and physics blocks play the analogous roles.

An alternative is Option A: pre-allocate the full `d = d_env + d_phys` dimensional Sobol sequence from the start, using the first `d_env` dimensions for the environmental suite and activating physics dimensions later. This provides fully nested extensibility in both `n` and `d` but requires knowing the complete parameter space upfront.

## Fairness Diagnostic

Sampling fairness is evaluated on the five design fields using quantile mismatch between the KDE-Sobol sample and the PHANGS reference set. For design field `a`:

```text
epsilon_a = max_{q in {0.1, 0.2, ..., 0.9}} |Q_q(log10 a_sample) - Q_q(log10 a_ref)|
```

The overall design fairness is:

```text
epsilon_env = max_{a in theta_env} epsilon_a
```

A design is considered fair if `epsilon_env <= epsilon_tol` with `epsilon_tol` typically in the range 0.10-0.20 dex.

Validation fields (`f_mol`, `Sigma_SFR`) are NOT included in this fairness criterion. They are checked separately after TIGRESS simulations are run by comparing the distribution of emergent simulation outputs against the PHANGS observed distributions.

## Comparison With KDE-LHS

The KDE-LHS method (previous implementation) and KDE-Sobol differ only in the final point-selection step:

- **KDE-LHS**: draw a large candidate pool from the KDE, convert to empirical rank space, select candidates nearest to LHS target points
- **KDE-Sobol**: map Sobol points directly through KDE marginal quantile functions (no candidate pool, no nearest-neighbor selection)

For a fixed `n`, both methods produce samples with similar one-dimensional coverage properties. The critical difference is extensibility: KDE-LHS has none, KDE-Sobol has the nested property.

For diagnostic purposes, observed-pixel LHS remains a useful baseline: it selects actual PHANGS pixels (no generative model), preserving real SFR and gas phase information, and can verify that the KDE-Sobol prior generates physically realistic parameter combinations.

## Implementation Notes

```python
import numpy as np
from scipy.stats import qmc
from scipy.interpolate import interp1d

# Generate Sobol sequence
sampler = qmc.Sobol(d=5, scramble=True, seed=42)
u = sampler.random_base2(m=6)  # n = 64 = 2^6

# Map through KDE marginal quantile functions (precomputed from fitted KDE)
# kde_marginal_quantile[j] is a callable Q_j: [0,1] -> log-space value
w = np.column_stack([kde_marginal_quantile[j](u[:, j]) for j in range(5)])

# Convert to physical space
theta = 10**w

# Apply physical bounds: qshear is column index 4
mask_valid = (theta[:, 4] <= 1.5) & np.all(theta > 0, axis=1)
# Replace invalid points from extended Sobol sequence as needed

# Coverage diagnostic
discrepancy = qmc.discrepancy(u)
```

The scrambled Sobol implementation in `scipy.stats.qmc.Sobol` uses the Owen (1998) scrambling by default, which breaks the correlation structure at higher dimensions while preserving the progressive nesting property. Always use `scramble=True` and fix `seed` for reproducibility.

For `n = 2^m` use `.random_base2(m)` which generates exactly `2^m` points with guaranteed optimal uniformity. For arbitrary `n` use `.random(n)`.
