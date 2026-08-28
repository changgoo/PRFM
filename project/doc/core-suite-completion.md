# Core (Paper I) suite — completion & analysis record

**Status:** completed — 32/32 models run and analyzed.
**Date:** 2026-08-21
**Branch:** `project/tigress-phangs-ili`

This is the master hand-off record for the first TIGRESS-PHANGS environmental
suite (the "core suite" designed in `project/doc/paper1-suite/`). It ties
together the design provenance, the run location, the design↔model mapping, the
external analysis tooling, and the key results, so the loop from
*PHANGS-informed design* → *TIGRESS-NCR run* → *PRFM analysis* is reproducible
from one page.

Related notes: [`prior-augmentation.md`](prior-augmentation.md),
[`physics-parameter-extension.md`](physics-parameter-extension.md); program
context in [`../TIGRESS-PHANGS-ILI-project.md`](../TIGRESS-PHANGS-ILI-project.md)
and [`../WORKFLOW.md`](../WORKFLOW.md).

---

## 1. What the suite is

A single-`Sigma_gas` environmental suite drawn from the PHANGS-informed prior
by KDE-Sobol conditional sampling. Physics parameters (`Z_gas`, `Z_dust`,
`xi_CR_amp`, feedback / population-synthesis assumptions) are held fixed at
fiducial (solar) values; only the **five environmental design fields** vary.

| Property | Value |
|---|---|
| Base model | `R8_8pc` (TIGRESS-NCR, 8 pc) |
| Target `Sigma_gas` | 10 M⊙ pc⁻² |
| Band half-width `delta` | 0.3 dex |
| Sample size `n` | 32 |
| Design fields | `Sigma_gas, Sigma_star, H_star, Omega, qshear` |
| Sampler | `synthesize_kde_sobol` (fixed Sobol seed, `kde_bandwidth_factor=1.0`, no augmentation) |
| Design CSV | `project/output/design_Sgas10.0_n0032.csv` |

Design-context figures (design over the PHANGS prior) live in
`project/doc/paper1-suite/paper-figures/design_Sgas10.0_n0032_{corner,marginals,selection}.png`.

> Note: this submitted design is tightly confined to the observed PHANGS
> marginals. The `H_star`/`Omega` coverage shortfall and the augmentation-block
> remedy are documented in [`prior-augmentation.md`](prior-augmentation.md); the
> augmentation block was **not** part of this completed core run.

## 2. Where the run lives

```
/tigress/changgoo/anvil/TIGRESS-NCR-suite/
├── R8_8pc_NCR_row0000/ … R8_8pc_NCR_row0031/   # 32 model directories
├── prfm_diagnostics/                            # PRFM analysis products (see §4)
├── hst_sfr_grid.png, hst_summary.png, …         # suite-level history figures
├── surface_density_evolution_*/                 # per-suite evolution figures
└── …
```

Model directory `R8_8pc_NCR_row<NNNN>` corresponds to **row `<NNNN>` (0-based)**
of the design CSV — see §3.

## 3. Design ↔ model mapping

The integer suffix in the model name equals the design CSV row index:

```
R8_8pc_NCR_row0002   ⟷   design_Sgas10.0_n0032.csv  row index 2
```

The join is exact (float round-off only), verified three independent ways
against `prfm_model_summary.csv`:

- `qshear`        matches the summary `qshear`            (Δ ≲ 5e-5);
- `Omega/1000`    matches the summary `omega` [km/s/pc]   (Δ ≲ 5e-8);
- `Sigma_star/(2 H_star)` matches `stellar_midplane_density` (Δ ≲ 3e-7).

Reusable join code: **`project/scripts/suite_analysis.py`**
(`load_model_summary`, `join_design_summary`, `load_joined`). It returns the
design augmented with `model`, `row`, `Sigma_SFR` (aliased to the chosen SFR
column), `rho_star`, and selected PRFM summary columns.

## 4. Analysis tooling & products

Analysis is done in the separate repo **`~/tigress_ncr_tools`** (not part of
`PRFM`). The relevant entry point and its full method write-up:

- CLI: `plot-suite-prfm /tigress/changgoo/anvil/TIGRESS-NCR-suite`
- Implementation: `src/tigress_ncr_tools/plot_suite_prfm.py`
- **Method reference: `~/tigress_ncr_tools/docs/prfm_analysis.md`** (authoritative
  for definitions, operators, and units — summarized below).

Products in `/tigress/changgoo/anvil/TIGRESS-NCR-suite/prfm_diagnostics/`:

| File | Contents |
|---|---|
| `prfm_model_summary.csv` | one row per model; snapshot mean + p16/p50/p84 for SFR, pressures, weights, yields; plus `omega`, `stellar_midplane_density`, `qshear` |
| `prfm_time_series.csv` | one row per (model, z-profile dump) |
| `prfm_vertical_profiles.csv` | one row per (model, height) |
| `model_sfr_colors.csv` | model rank / plotted color |
| `prfm_pressure_weight_relations*.png` | P_tot vs W, P vs Σ_SFR, W vs Σ_SFR |
| `prfm_pressure_components_yields*.png` | 4 pressure + 4 yield panels |
| `prfm_delta_pressure_weight_relations*.png` | pressure-drop variant |
| `prfm_vertical_profiles*.png` | density + stress profiles (2p and total gas) |
| `prfm_pressure_weight_time_evolution.png` | 200–600 Myr tracks |

`*` = also `_color_by_omega`, `_color_by_stellar_midplane_density`,
`_color_by_qshear` variants.

**Conventions to remember** (from `prfm_analysis.md`):

- Two-phase (warm+cold, phases 7/11/12/13) midplane pressure vs **whole-gas**
  vertical weight — deliberately different phase selections.
- `pressure_total` = turbulent + thermal + δB + mean-field magnetic; **excludes**
  radiation, cosmic-ray, and hot-gas pressure.
- Pressure/weight columns are `P/k_B` in K cm⁻³; yields Υ in km/s.
- Yields use `sfr40`; model color uses time-weighted `sfr10` (200–600 Myr).
- Cache window 200–600 Myr; summary statistics window **400–600 Myr**.
- Snapshot summary means are **not** cadence-weighted.

## 5. New figure — resulting Σ_SFR over the PHANGS parameter space

**Script:** `project/scripts/plot_sfr_corner.py`
**Figure:** `project/doc/paper1-suite/paper-figures/design_Sgas10.0_n0032_sfr_corner.png`

A **6×6** log-space corner that closes the design loop and validates the suite
against observations. The six fields are the five native design fields **plus
`Sigma_SFR`** as a dimension:

- **gray cloud / histograms** — the PHANGS reference band (`Sigma_gas ≈ 10`,
  791 pixels) the design was drawn from. On the `Sigma_SFR` axis this is the
  *observed* PHANGS `Sigma_SFR` (default Hα+W4 recalibrated,
  `Sigma_SFR_HaW4recal`; 358 finite pixels in-band);
- **colored points** — the 32 design models, positioned at their **simulated**
  `Sigma_SFR` (`sfr40_mean`, 400–600 Myr mean) and colored by the same value;
- **star** — the R8 fiducial (design fields only; no fiducial SFR is drawn).

The bottom row (`Sigma_SFR` vs each parameter) and the bottom-right diagonal
(simulated points on the observed distribution) turn this into a direct
simulated-vs-observed SFR comparison, not just a design-placement figure.

Reproduce:

```bash
PY=~/.conda/envs/pyathena/bin/python   # any env with pandas+matplotlib+scipy+astropy
$PY project/scripts/plot_sfr_corner.py \
    project/output/design_Sgas10.0_n0032.csv \
    --summary /tigress/changgoo/anvil/TIGRESS-NCR-suite/prfm_diagnostics/prfm_model_summary.csv \
    --output-dir project/doc/paper1-suite/paper-figures
# --sfr-col sfr10_mean      color/position by the 10 Myr SFR instead of 40 Myr
# --obs-sfr-col Sigma_SFR_FUVW4recal   use a different observed SFR tracer
```

The PHANGS background requires the megatable (`config/phangs_prfm.yml` →
`data/phangs_megatable/`, hexagon + gauss apertures):
`python scripts/download_phangs.py --aperture hexagon` and `--aperture gauss`.

**What it shows:** the simulated `Sigma_SFR` overlaps the observed PHANGS
distribution at `Sigma_gas ≈ 10` (simulated median log ≈ −2.53 vs observed
Hα+W4 ≈ −2.55) and follows the observed `Sigma_gas`–`Sigma_SFR` and
`Sigma_star`–`Sigma_SFR` trends; the suite does not reach the observed
high-SFR tail (log `Sigma_SFR ≳ −2.0`), and a few low-`rho_star` models fall at
or just below the observed low-SFR edge.

## 6. Key results (summary window 400–600 Myr, n = 32)

| Quantity | min | median | max |
|---|---|---|---|
| `Sigma_SFR` (sfr40) [M⊙ kpc⁻² yr⁻¹] | 3.5e-4 | 2.9e-3 | 9.1e-3 |
| `P_tot,2p / W` (pressure–weight ratio) | 0.98 | 1.14 | 1.71 |
| `Upsilon_tot` [km/s] | 866 | 1574 | 8986 |
| `Omega` [km/s/pc] | 0.012 | 0.029 | 0.078 |
| `rho_star` [M⊙ pc⁻³] | 0.016 | 0.093 | 0.645 |
| `qshear` | 0.44 | 0.91 | 1.28 |

- Resulting `Sigma_SFR` spans **~1.4 dex** across the fixed-`Sigma_gas` band —
  driven by the variation in `Sigma_star`/`H_star` (→ `rho_star`), `Omega`, and
  `qshear`, not by `Sigma_gas`.
- Vertical dynamical equilibrium holds: `P_tot,2p ≈ W` (median ratio 1.14).
- Component yield medians [km/s]: turbulent 546, thermal 312, δB magnetic 258,
  mean-field magnetic 369.
- The large `Upsilon_tot` tail (up to ~9000) comes from the lowest-SFR models,
  where `P/Σ_SFR` is inflated; inspect medians, not means, for those.

## 7. Open items / next steps

- Prior augmentation (`H_star`, `Omega` coverage) — block generated but not yet
  run; see [`prior-augmentation.md`](prior-augmentation.md).
- Physics-parameter extension (`Z_gas`, `Z_dust`, `xi_CR_amp`) — Phase 1
  ablation designed; see [`physics-parameter-extension.md`](physics-parameter-extension.md).
- Validation of `Sigma_atom`, `Sigma_mol`, `f_mol` against PHANGS (held as
  validation quantities, not inputs).
- Feed the joined design + outcome table into emulator / ILI training
  (Sub-project goals in the program doc).
