#!/usr/bin/env python
"""Run KDE-Sobol sampling for one or more Sigma_gas targets.

Usage:
    python project/scripts/run_sampling.py
    python project/scripts/run_sampling.py --sigma-gas 10 15 20
    python project/scripts/run_sampling.py --n-samples 64 128 256 --sigma-gas 10

Outputs:
    project/output/design_Sgasxx.x_nYYY.csv   one CSV per (target, n) pair
    project/output/fairness_summary.csv        per-field quantile mismatch table
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Ensure project root is on sys.path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from prfm.phangs import load_configured_phangs
from prfm.phangs_sampling import PHANGSSamplingDesigner, SamplingConfig

OUTPUT_DIR = ROOT / "project" / "output"
CONFIG_PATH = ROOT / "config" / "phangs_prfm.yml"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--sigma-gas", nargs="+", type=float,
        default=[10.0],
        help="Target Sigma_gas values in M_sun/pc^2 (default: 10)",
    )
    p.add_argument(
        "--delta", type=float, default=0.3,
        help="Half-width of log10 Sigma_gas band in dex (default: 0.3)",
    )
    p.add_argument(
        "--n-samples", nargs="+", type=int,
        default=[64, 128, 256],
        help="Sample sizes to generate (default: 64 128 256)",
    )
    p.add_argument(
        "--config", type=Path, default=CONFIG_PATH,
        help="Path to phangs_prfm.yml config file",
    )
    p.add_argument(
        "--seed", type=int, default=42,
        help="Sobol seed (default: 42)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading PHANGS data from config: {args.config}")
    result = load_configured_phangs(str(args.config), base_dir=str(ROOT))
    table = result["table"]
    print(f"  Loaded {len(table)} apertures")

    cfg = SamplingConfig(
        delta_sigma_gas=args.delta,
        sobol_seed=args.seed,
        kde_aux_sample_size=100_000,
    )
    designer = PHANGSSamplingDesigner(table, config=cfg)

    fairness_rows = []

    for sigma_gas_target in args.sigma_gas:
        print(f"\n=== Sigma_gas target = {sigma_gas_target:.1f} M_sun/pc^2 ===")
        designer = designer.with_config(target_sigma_gas=sigma_gas_target)
        reference = designer.select_reference_pixels()
        print(f"  Reference pixels: {len(reference)}")

        if len(reference) < 20:
            print(f"  WARNING: fewer than 20 reference pixels — skipping.")
            continue

        for n in args.n_samples:
            print(f"  Sampling n={n} ...", end=" ", flush=True)
            sample = designer.synthesize_kde_sobol(reference, n_samples=n)
            n_extra = sample.attrs.get("n_extra", 0)
            print(f"done (n_extra={n_extra})")

            # Fairness metric
            max_err, per_field = designer.quantile_error_report(
                reference, sample, fields=cfg.design_fields
            )
            print(f"    fairness max_err = {max_err:.3f} dex")

            for _, row in per_field.iterrows():
                fairness_rows.append({
                    "sigma_gas_target": sigma_gas_target,
                    "n_samples": n,
                    "field": row["field"],
                    "max_quantile_err_dex": row["max_quantile_error_dex"],
                })

            # Save design table
            out_path = OUTPUT_DIR / (
                f"design_Sgas{sigma_gas_target:.1f}_n{n:04d}.csv"
            )
            sample.to_csv(out_path, index=False)
            print(f"    saved → {out_path.relative_to(ROOT)}")

    # Save fairness summary
    if fairness_rows:
        fair_df = pd.DataFrame(fairness_rows)
        fair_path = OUTPUT_DIR / "fairness_summary.csv"
        fair_df.to_csv(fair_path, index=False)
        print(f"\nFairness summary → {fair_path.relative_to(ROOT)}")
        print(fair_df.to_string(index=False))


if __name__ == "__main__":
    main()
