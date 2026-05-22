# API Reference

## `prfm.prfm` — Core PRFM model

The main computation module. All physical quantities use CGS units internally.

### `PRFM` class

Object-oriented wrapper around the functional layer with automatic unit conversion (CGS ↔ astronomical).

```python
from prfm import PRFM

model = PRFM()
result = model.get_self_consistent_solution(Sigma_gas, Sigma_star, Omega_d, H_star)
```

### Scale height solvers

| Function | Description |
|----------|-------------|
| `get_scale_height(Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff, ...)` | Dispatcher — selects analytic or numerical solver based on which weight terms are active |
| `get_scale_height_analytic(...)` | Full analytic solution of the cubic vertical equilibrium equation (all three weight terms) |
| `get_scale_height_numerical(...)` | Numerical solution via Brent's method (vectorized via `np.vectorize`) |
| `get_scale_height_thick(Sigma_gas, rho_star, sigma_eff)` | Approximate solution for H★ ≫ H_gas using stellar volume density |
| `get_scale_height_thin(Sigma_gas, Sigma_star, sigma_eff)` | Approximate solution for H★ ≪ H_gas |
| `get_scale_height_gas_only(...)` | Analytic limit: gas self-gravity only |
| `get_scale_height_star_only(...)` | Analytic limit: stellar gravity only |
| `get_scale_height_dm_only(...)` | Analytic limit: dark matter gravity only |
| `get_scale_height_star_gas(...)` | Analytic solution neglecting dark matter |
| `get_scale_height_dm_gas(...)` | Analytic solution neglecting stars |

### Weight functions

| Function | Description |
|----------|-------------|
| `get_weight_gas(Sigma_gas)` | Gas self-gravity weight: π G Σ²/2 |
| `get_weight_star(Sigma_gas, H_gas, Sigma_star, H_star)` | Stellar gravity weight |
| `get_weight_star_gaussian(Sigma_gas, H_gas, Sigma_star, H_star)` | Stellar weight assuming Gaussian profiles |
| `get_weight_dm(Sigma_gas, H_gas, Omega_d, zeta_d)` | Dark matter gravity weight |
| `get_pressure(Sigma_gas, H_gas, sigma_eff)` | Midplane pressure P = σ² Σ / (2 H) |
| `get_weights(Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff, ...)` | Compute H_gas and all three weight terms |
| `get_weight_contribution(Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff, ...)` | Return fractional weight contributions (f_gas, f_star, f_dm) |
| `get_weights_thick(Sigma_gas, rho_star, sigma_eff)` | Weights in the thick-stellar-disk limit |
| `get_weights_thin(Sigma_gas, Sigma_star, sigma_eff)` | Weights in the thin-stellar-disk limit |

### Feedback models

| Function | Description |
|----------|-------------|
| `get_sigma_eff(P, Z, model)` | Pressure-dependent effective velocity dispersion |
| `get_feedback_yield(P, Z, model)` | Total feedback yield Υ(P) [km s⁻¹] |
| `get_feedback_yield_comp(P, Z, comp, model)` | Per-component feedback yield (thermal/turbulent/magnetic) |
| `get_sfr(P, Z, Ytot)` | SFR surface density from pressure and yield |
| `get_self_consistent_solution(Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff, ...)` | Iteratively solve for converged (W_tot, H_gas, σ_eff) |

---

## `prfm.phangs` — PHANGS megatable utilities

Tools for downloading, loading, and analysing the PHANGS megatable (Sun et al. 2022/2023 v4.0).

### Data access

| Function | Description |
|----------|-------------|
| `list_galaxies()` | Return list of galaxy names in the v4.0 release |
| `list_files(data_dir, aperture)` | List ECSV files for a given aperture type |
| `download(data_dir, aperture, overwrite)` | Download PHANGS megatable files from CANFAR |
| `load(fpath)` | Load a single ECSV file as an Astropy Table |
| `load_all(data_dir, aperture)` | Load all ECSV files and vstack into one Table |

### PRFM input computation

| Function | Description |
|----------|-------------|
| `compute_prfm_inputs(table)` | Derive Σ_gas, Ω_d, H★ from raw megatable columns |
| `valid_rows(table, cols, rel_error)` | Boolean mask selecting rows with finite, positive PRFM inputs |
| `run_prfm(table, prfm_cols, sigma_eff_model, yield_model, zprime_col, variation)` | Apply PRFM model to every row; adds P_weight, H_gas, sigma_eff_sol, Sigma_SFR_pred columns |
| `get_weights(table, variation)` | Return (f_gas, f_star, f_dm) weight fractions for rows with PRFM solutions |

---

## `prfm.phangs_plot` — PHANGS plotting utilities

Centralized plotting helpers for PHANGS analysis notebooks.

| Function | Description |
|----------|-------------|
| `scatter_plot(t, xcol, ycol, ...)` | Log–log scatter plot with optional error bars, slices, and galaxy colouring |
| `hist_plot(t, cols, ...)` | Diagonal-and-off-diagonal scatter/histogram matrix |
| `plot_weights(tbl, ax, variation)` | Scatter + binned-mean weight fractions vs P_DE |
| `col_label(col)` | Human-readable axis label for a megatable column name |
| `resolve_columns(t, xcol, ycol, ecol)` | Resolve column names and extract data arrays with units |
| `get_symmetric_log_errorbars(linear_mean, linear_std)` | Convert linear-space std to visually symmetric log-axis error bars |

---

## `prfm.simulations` — Simulation data utilities

Utilities for loading and plotting TIGRESS simulation results.

### `PRFM_data` class

Container for simulation pressure and SFR data from a single paper/run.

```python
data = PRFM_data("KO15")
data.SFR = ...
data.convert_linear_log()
data.get_Ptotal()
data.get_yield()
```

**Methods:**

| Method | Description |
|--------|-------------|
| `convert_log_linear()` | Convert log10-stored fields to linear; propagate uncertainties |
| `convert_linear_log()` | Convert linear fields to log10; propagate uncertainties |
| `get_Ptotal()` | Compute Ptot, Pnonth, Pimag from component pressures |
| `get_yield()` | Compute feedback yield Y = P/SFR [km s⁻¹] for each pressure field |

### Data loaders

| Function | Description |
|----------|-------------|
| `load_sim_data()` | Load all simulation datasets (pre-TIGRESS, classic, GC, NCR) and return as a dict |
| `load_pretigress()` | Load KO15, KOK13, KKO11 pre-TIGRESS datasets |
| `load_classic_data(from_table)` | Load TIGRESS-classic simulation data |
| `load_gc_data(from_table)` | Load TIGRESS galactic-centre simulation data |
| `load_ncr_data(fname)` | Load TIGRESS-NCR data from NetCDF file |

### Plotting utilities

| Function | Description |
|----------|-------------|
| `add_one_sim(data, xf, yf, ms, ...)` | Plot a single simulation dataset with error bars |
| `add_ncr_sim(data, xf, yf, ms, legend)` | Plot NCR simulation with metallicity colour coding |
| `add_PSFR_model_line(Wmin, Wmax, Z, model, ...)` | Overlay a P–SFR model line |
| `add_PSFR_model_lines(Wmin, Wmax, model)` | Overlay multiple P–SFR model lines |
| `add_PSFR_ncr_model_lines(Wmin, Wmax, ...)` | Overlay NCR P–SFR model lines for multiple metallicities |
| `add_yield_model_line(Wmin, Wmax, Z, comp, model, ...)` | Overlay a feedback yield model line |
| `add_yield_model_lines(Wmin, Wmax, comp)` | Overlay multiple yield model lines |
| `add_yield_ncr_model_lines(Wmin, Wmax, comp)` | Overlay NCR yield model lines for multiple metallicities |
| `setup_axes(nrow, figsize, width_ratios)` | Create figure with main axes and side-panel grid |
| `get_log_errorbars(log_mean, log_std)` | Convert log-space statistics to matplotlib asymmetric error bars |
| `get_ncr_color(Z, cmap, Zmin, Zmax)` | Map metallicity value to colour via LogNorm |
