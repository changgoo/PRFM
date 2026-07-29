# Augmenting the environmental prior (H_star, Omega)

**Status:** implemented (augmentation-block generator + Sgas10 n=32 block).
**Date:** 2026-07-29

## The issue

The submitted suite was drawn with `synthesize_kde_sobol` at
`kde_bandwidth_factor = 1.0` and **no augmentation** — `run_sampling.py` never
enabled the "expanded" path. As a result the design marginals are tightly
confined to the observed PHANGS distribution. In particular `H_star` and
`Omega` needed roughly a factor of two more coverage than PHANGS provides at
`Sigma_gas ≈ 10`.

Note that the pre-existing "expanded" machinery does **not** do what we want:

- `synthesize_expanded_kde_sobol` broadens via **KDE bandwidth**, which widens
  *all five* fields (and would break the `Sigma_gas` band and the `q ≤ 1.5`
  bound).
- The per-field `expanded_prior_lower_dex` config was only wired into the
  **LHS** path (`synthesize_expanded_kde_lhs`), not the Sobol design.

## Mechanism: selective per-field marginal stretch

Because the Sobol design maps each field through its **independent marginal
quantile function**, we can broaden *only* the chosen fields. `synthesize_kde_sobol`
now takes `augment_dex` (default `config.augment_prior_dex`): each named field's
quantile function is stretched about its log-median so the tails extend by
`~delta` dex on each side, leaving the other fields — and thus the `Sigma_gas`
band and `q` bound — untouched.

For this suite: `augment_dex = {"H_star": 0.3, "Omega": 0.3}` (≈ ×2 per tail).
Empty `augment_dex` reproduces the previous design byte-for-byte.

## Can new samples be added to the existing suite?

Yes — but **not** by "increasing n". Two distinct operations:

- **Densify** (same narrow prior, larger n, same seed via `random_base2`):
  genuinely nested, `S(32) ⊂ S(64)`, old points preserved — but only packs more
  points into the same PHANGS-limited range.
- **Augment** (broader H_star/Omega): changes the marginal quantile functions,
  hence the Sobol→physical map. A fresh augmented design's first 32 points are
  **not** the existing 32, so the current runs are not a subset of it.

So augmentation is added as a **separate block**: keep all 32 runs, and draw N
new samples from the augmented prior with an **independent Sobol seed** (so they
don't overlap the existing points), appended as new rows (`row0032 …`). The
existing simulations stay valid; the block extends coverage into the broadened
tails.

**Inference caveat.** The union is then a *mixture* of two proposals (narrow
core + augmented block). For amortized SBI this is fine — what matters is that
the training set *covers* the augmented prior support; the extra density in the
core is harmless. This mirrors CAMELS combining 1P + LH suites. If instead a
clean i.i.d. draw from a single prior is required, the only alternative is
re-drawing the whole design from the augmented prior and re-running (discarding
the 32).

## Tooling

- `project/scripts/augment_design.py` — draws the augmentation block (independent
  seed, continued row suffixes) and writes both the block CSV (for new SLURM
  rows) and a combined CSV (for plots/record).
- `project/scripts/csv_to_slurm_yaml.py` — turns the block CSV into a suite YAML
  containing only the new rows.
- `project/scripts/plot_augmentation_coverage.py` — narrow vs block vs PHANGS
  marginals, with axes spanning the full design range so the broadened tails are
  visible (`plot_design_from_csv.py` locks axes to the PHANGS range and clips
  them).

Generated Phase artifacts (Sgas10, n=32 block):
`project/output/augment_design_Sgas10.0_n0032_n0032.csv`,
`project/output/design_Sgas10.0_n0032_plus_aug0032.csv`,
`project/suites/augment_design_Sgas10.0_n0032_n0032.yml`.

## Preventing recurrence

`run_sampling.py` still produces *unaugmented* core designs. To make future
production designs augmented, set `SamplingConfig.augment_prior_dex` (e.g.
`{"H_star": 0.3, "Omega": 0.3}`) — a `--augment` flag on `run_sampling.py` would
expose this on the command line (not yet added).
