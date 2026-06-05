#!/usr/bin/env python
"""Generate KDE-Sobol sampling diagnostic figures.

Usage:
    python project/scripts/plot_sampling_diagnostics.py
    python project/scripts/plot_sampling_diagnostics.py --sigma-gas 10 --n-samples 64 128

Outputs (in figures/phangs/sobol/Sgas<X>/) :
    correlation_matrix.png    PHANGS distribution with Sigma_gas band overlay
    distributions_n<N>.png    1-D marginal overlays: reference vs Sobol sample
    pairs_n<N>.png            Pairwise scatterplot overlays
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from prfm.phangs import load_configured_phangs
from prfm.phangs_sampling import PHANGSSamplingDesigner, SamplingConfig

CONFIG_PATH = ROOT / "config" / "phangs_prfm.yml"
FIGURE_DIR = ROOT / "figures" / "phangs" / "sobol"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--sigma-gas", nargs="+", type=float, default=[10.0],
        help="Target Sigma_gas values (default: 10)",
    )
    p.add_argument(
        "--delta", type=float, default=0.3,
        help="Half-width of log10 Sigma_gas band in dex (default: 0.3)",
    )
    p.add_argument(
        "--n-samples", nargs="+", type=int, default=[64, 128],
        help="Sample sizes to plot (default: 64 128)",
    )
    p.add_argument(
        "--config", type=Path, default=CONFIG_PATH,
    )
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--dpi", type=int, default=150)
    return p.parse_args()


def save(fig: plt.Figure, path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved → {path.relative_to(ROOT)}")


def main() -> None:
    args = parse_args()

    print("Loading PHANGS data ...")
    result = load_configured_phangs(str(args.config), base_dir=str(ROOT))
    table = result["table"]
    print(f"  {len(table)} apertures loaded")

    cfg = SamplingConfig(
        delta_sigma_gas=args.delta,
        sobol_seed=args.seed,
        kde_aux_sample_size=100_000,
    )
    designer = PHANGSSamplingDesigner(table, config=cfg)

    for sigma_gas_target in args.sigma_gas:
        print(f"\n=== Sigma_gas = {sigma_gas_target:.1f} ===")
        tag = f"Sgas{sigma_gas_target:.1f}"
        out_dir = FIGURE_DIR / tag
        designer = designer.with_config(target_sigma_gas=sigma_gas_target)
        reference = designer.select_reference_pixels()
        print(f"  {len(reference)} reference pixels")

        # 1. Correlation matrix with Sigma_gas band overlay
        print("  Plotting correlation matrix ...")
        fig, axes, _ = designer.plot_correlation_matrix(
            targets=[sigma_gas_target],
            delta_sigma_gas=args.delta,
        )
        save(fig, out_dir / "correlation_matrix.png", args.dpi)

        for n in args.n_samples:
            print(f"  Sampling n={n} ...")
            sample = designer.synthesize_kde_sobol(reference, n_samples=n)

            # 2. Distribution overlays (design fields only)
            print(f"  Plotting distributions n={n} ...")
            fig, axes = designer.plot_distribution_overlay(
                reference, sample,
                fields=cfg.design_fields,
            )
            fig.suptitle(
                rf"KDE-Sobol $n={n}$, "
                rf"$\Sigma_{{\rm gas}}={sigma_gas_target:.1f}\,M_\odot\,{{\rm pc}}^{{-2}}$",
                fontsize="medium",
            )
            save(fig, out_dir / f"distributions_n{n:04d}.png", args.dpi)

            # 3. Pairwise overlays
            print(f"  Plotting pairs n={n} ...")
            fig, axes = designer.plot_pair_overlay(
                reference, sample,
                fields=cfg.design_fields,
            )
            fig.suptitle(
                rf"KDE-Sobol $n={n}$, "
                rf"$\Sigma_{{\rm gas}}={sigma_gas_target:.1f}\,M_\odot\,{{\rm pc}}^{{-2}}$",
                fontsize="medium",
            )
            save(fig, out_dir / f"pairs_n{n:04d}.png", args.dpi)

    print("\nDone.")


if __name__ == "__main__":
    main()
