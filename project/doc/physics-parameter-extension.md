# Extending the suite with physics parameters (Z_gas, Z_dust, xi_CR_amp)

**Status:** Phase 1 (ablation) implemented; Phase 2 (joint SBI design) planned.
**Date:** 2026-07-29

The first TIGRESS-PHANGS suite varies five *environmental* parameters
(`Sigma_gas, Sigma_star, H_star, Omega, qshear`) and holds all *physics*
parameters fixed at fiducial (solar) values. This note records how and why we
extend the suite to vary three physics parameters:

| Parameter    | Range     | athinput key        | Meaning                                   |
|--------------|-----------|---------------------|-------------------------------------------|
| `Z_gas`      | 0.1 – 1   | `problem/Z_gas`     | Gas-phase metallicity (solar-normalized)  |
| `Z_dust`     | derived   | `problem/Z_dust`    | Dust abundance = `f_dtm · Z_gas`          |
| `xi_CR_amp`  | 0.1 – 1   | `problem/xi_CR_amp` | CR ionization-rate amplitude (fiducial 1) |

with the **dust-to-metal ratio** `f_dtm = Z_dust / Z_gas ∈ [0.1, 1]` as the
actual free parameter. All three are dimensionless and solar-normalized; the
fiducial suite sits at `(Z_gas, f_dtm, xi_CR_amp) = (1, 1, 1)`.

## Why these are sampled differently from the five environmental parameters

The five environmental fields are **observationally correlated** — the pipeline
fits a KDE to the PHANGS pixels in a `Sigma_gas` band and maps a Sobol sequence
through the marginal quantile functions (`phangs_sampling.synthesize_kde_sobol`).
The three physics parameters are a different kind of quantity:

- **`Z_gas`** — metallicity *is* observed (`Zprime`), but at `Sigma_gas ≈ 10`
  the PHANGS pixels are nearly all near-solar (massive spirals). A KDE
  conditioned on PHANGS would never reach `Z ~ 0.1`. To explore a full decade
  sub-solar we must use an **independent (expanded) prior**, not the observed
  correlation.
- **`Z_dust`** — not free; parameterized through the ratio `f_dtm`, which makes
  `Z_dust ≤ Z_gas` automatic and matches the physically meaningful knob.
- **`xi_CR_amp`** — a pure model knob with no PHANGS constraint.

All three span a decade, so we sample them **log-uniform** over `[0.1, 1]`.

## Design: product of two low-discrepancy blocks

We keep the environmental block exactly as-is and **append** an independent
block for the physics parameters:

```
environmental design  =  synthesize_kde_sobol(d=5, seed=42)          [UNCHANGED]
physics design        =  independent Sobol(d=3, log-uniform priors)  [NEW]
suite                 =  column-concat  ->  Z_dust = f_dtm * Z_gas
```

A *product* (rather than one joint 8-D KDE-Sobol) is chosen because:

1. **It preserves the existing environmental design byte-for-byte.** The `d=5`
   KDE-Sobol call is untouched (same seed, same qshear rejection path), so its
   rows — and the already-submitted simulations — remain valid.
2. **It matches the two natures of the parameters:** observationally-conditioned
   environment vs. independent model priors. For inference we *want* the physics
   marginals flat and decoupled; the KDE would fight that.
3. Each block stays independently **nestable in `n`** (both use `random_base2`
   with a fixed seed), so the suite can still grow n = 64 → 128 → 256.

The existing 32 runs (all at `(1,1,1)`) become the **fiducial-metallicity anchor
slice**: for any environmental point they provide the solar-metallicity, fiducial
reference against which varied-Z/dust/CR runs are compared.

## Two phases

### Phase 1 — ablation (implemented)

A one-at-a-time (OAT) grid over `(Z_gas, f_dtm, xi_CR_amp)` at a few
**environmental anchor points drawn from the existing design**:

- **Anchors:** rows spanning the physical regime (default: 3 rows spanning
  `Sigma_star`). Their `(1,1,1)` base runs already exist, so each ablation curve
  gets its reference point for free, and the anchors coincide with real runs.
- **Sweep:** hold two knobs at 1.0, step the third through the sub-fiducial
  levels `{0.1, 0.22, 0.46}` (log-spaced; 1.0 is the existing base).
- **Cost:** 3 params × 3 levels × 3 anchors = **27 runs** (bases reused).
- **Interactions** (e.g. low-Z *and* low-dust together) are intentionally *not*
  in the OAT grid — those come for free in Phase 2.

Job IDs encode provenance, e.g. `R8_8pc_NCR_a0019_Z0p10` = anchor row 19, `Z_gas`
sweep, level 0.10.

Tooling: `project/scripts/make_ablation_grid.py` emits the ablation CSV;
`project/scripts/csv_to_slurm_yaml.py` turns it into a suite YAML (SLURM as
usual). Generated Phase-1 outputs for the `Sgas10, n0032` design live at
`project/output/ablation_design_Sgas10.0_n0032.csv` and
`project/suites/ablation_design_Sgas10.0_n0032.yml`.

### Phase 2 — joint SBI design (planned)

Extend `phangs_sampling` with independent log-uniform priors and a
`synthesize_independent_sobol` that column-concatenates onto the untouched 5-D
KDE-Sobol block (the product design above), at larger `n` (target 128–256, since
going 5 → 8 dimensions triples the volume).

**Phase 1 reuse:** because both phases draw the physics parameters from the *same
log-uniform support*, the ablation runs are valid samples from the joint prior
and fold into the Phase-2 SBI training set. If the ablation anchors are included
as environmental points in the product design, the two phases join seamlessly.

## Tooling changes (both phases)

`csv_to_slurm_yaml.py` gained an `OPTIONAL_PARAM_MAP`: if `Z_gas`, `Z_dust`, or
`xi_CR_amp` columns are present in a design CSV they vary per row (and their
fixed defaults are dropped); if absent, the fixed defaults apply. It also honors
an optional `suffix` column for descriptive job IDs. A physical-only 5-column
CSV reproduces the previous output byte-for-byte.
