#!/usr/bin/env python
"""Plot KDE-Sobol sample coverage of the 5 environmental design fields.

Produces three figures per Sigma_gas target:
  design_fields_selection.png        Reduced correlation matrix: full PHANGS
                                     (gray) + selected band (blue)
  design_fields_corner_n<N>.png      5x5 corner: reference + standard samples
  design_fields_corner_expanded.png  5x5 corner: reference + expanded samples
  design_fields_distributions.png    1-D marginals: standard vs expanded

Usage:
    python project/scripts/plot_design_fields.py
    python project/scripts/plot_design_fields.py --sigma-gas 10 --n-samples 32 64 128 256
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
plt.style.use(ROOT / "project" / "scripts" / "prfm.mplstyle")
sys.path.insert(0, str(ROOT))

from prfm.phangs import load_configured_phangs
from prfm.phangs_sampling import PHANGSSamplingDesigner, SamplingConfig

CONFIG_PATH = ROOT / "config" / "phangs_prfm.yml"
FIGURE_DIR  = ROOT / "figures" / "phangs" / "sobol"

FIELD_LABELS = {
    "Sigma_gas":  r"$\log\,\Sigma_{\rm gas}$"  + "\n" + r"[$M_\odot\,{\rm pc}^{-2}$]",
    "Sigma_star": r"$\log\,\Sigma_\star$"       + "\n" + r"[$M_\odot\,{\rm pc}^{-2}$]",
    "H_star":     r"$\log\,H_\star$"            + "\n" + r"[pc]",
    "Omega":      r"$\log\,\Omega$"             + "\n" + r"[km s$^{-1}$ kpc$^{-1}$]",
    "qshear":     r"$\log\,q$",
}

# Colors: gray=all PHANGS, blue=selected band, orange/green=samples
COLOR_ALL  = "0.78"
COLOR_REF  = "#3a86ff"   # selected band
COLOR_STD  = ["#ff6b35", "#06d6a0", "#118ab2", "#ffd166"]  # n=32,64,128,256
COLOR_EXP  = "#9b2226"   # expanded prior (single color, one n)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sigma-gas",  nargs="+", type=float, default=[10.0])
    p.add_argument("--delta",      type=float, default=0.3)
    p.add_argument("--n-samples",  nargs="+", type=int,   default=[32, 64, 128, 256])
    p.add_argument("--n-expanded", type=int,  default=64,
                   help="Sample size for expanded-prior comparison (default: 64)")
    p.add_argument("--config", type=Path, default=CONFIG_PATH)
    p.add_argument("--seed",   type=int,  default=42)
    p.add_argument("--dpi",    type=int,  default=150)
    return p.parse_args()


def save(fig: plt.Figure, path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved → {path.relative_to(ROOT)}")


# ── helpers ──────────────────────────────────────────────────────────────────

def _log(table_or_df, field: str) -> np.ndarray:
    """Return log10 values of a field, masking non-positive/non-finite."""
    vals = np.asarray(table_or_df[field], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        lv = np.log10(vals)
    lv[~np.isfinite(lv)] = np.nan
    return lv


def _axis_range(log_vals: np.ndarray,
                pad: float = 0.05) -> tuple[float, float]:
    lo, hi = np.nanpercentile(log_vals, [0.5, 99.5])
    span = hi - lo
    return lo - pad * span, hi + pad * span


def _corner_axes(fields: list[str],
                 size: float = 2.2) -> tuple[plt.Figure, np.ndarray]:
    n = len(fields)
    fig, axes = plt.subplots(n, n, figsize=(size * n, size * n),
                             constrained_layout=False)
    fig.subplots_adjust(hspace=0.08, wspace=0.08,
                        left=0.12, bottom=0.10, right=0.98, top=0.95)
    return fig, axes


def _decorate_corner(axes: np.ndarray, fields: list[str]) -> None:
    """Apply axis labels and hide upper triangle."""
    n = len(fields)
    for row in range(n):
        for col in range(n):
            ax = axes[row, col]
            if col > row:
                ax.set_visible(False)
                continue
            ax.tick_params(labelsize="small")
            if row == n - 1:
                ax.set_xlabel(FIELD_LABELS[fields[col]], fontsize="small")
                ax.xaxis.set_major_locator(ticker.MaxNLocator(3, prune="both"))
            else:
                ax.tick_params(labelbottom=False)
            if col == 0 and row > 0:
                ax.set_ylabel(FIELD_LABELS[fields[row]], fontsize="small")
                ax.yaxis.set_major_locator(ticker.MaxNLocator(3, prune="both"))
            else:
                ax.tick_params(labelleft=False)


# ── figure 1: reduced correlation matrix ──────────────────────────────────────

def plot_selection(
    full_table,
    reference,
    fields: list[str],
    sigma_gas_target: float,
    delta: float,
    dpi: int,
    out_dir: Path,
) -> None:
    """Reduced correlation matrix: all PHANGS (gray) + selected band (blue)."""
    n_fields = len(fields)
    fig, axes = _corner_axes(fields)

    # Pre-compute log arrays for full table and reference
    log_full = {f: _log(full_table, f) for f in fields}
    log_ref  = {f: _log(reference,   f) for f in fields}

    # Ranges from full table
    ranges = {f: _axis_range(log_full[f]) for f in fields}
    n_bins = 30

    for row, yf in enumerate(fields):
        for col, xf in enumerate(fields):
            ax = axes[row, col]
            if col > row:
                ax.set_visible(False)
                continue

            xlo, xhi = ranges[xf]
            bins_x = np.linspace(xlo, xhi, n_bins + 1)

            if col == row:
                ax.hist(log_full[xf], bins=bins_x, density=True,
                        color=COLOR_ALL, edgecolor="none")
                ax.hist(log_ref[xf],  bins=bins_x, density=True,
                        histtype="step", color=COLOR_REF, linewidth=1.5)
                ax.set_yticks([])
            else:
                ylo, yhi = ranges[yf]
                ax.scatter(log_full[xf], log_full[yf],
                           s=0.8, c=COLOR_ALL, rasterized=True, zorder=1)
                ax.scatter(log_ref[xf],  log_ref[yf],
                           s=1.5, c=COLOR_REF, rasterized=True, zorder=2,
                           alpha=0.7)
                ax.set_ylim(ylo, yhi)

            ax.set_xlim(xlo, xhi)

    _decorate_corner(axes, fields)

    # Legend via proxy
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor=COLOR_ALL, label=f"PHANGS ({len(full_table):,})"),
        Line2D([0], [0], color=COLOR_REF, linewidth=2,
               label=(rf"$|\log\Sigma_{{\rm gas}}-\log{sigma_gas_target:.0f}|"
                      rf"\leq{delta}\,\rm dex$ ({len(reference):,})")),
    ]
    axes[0, 0].legend(handles=handles, fontsize="small",
                      loc="upper right", framealpha=0.85)

    fig.suptitle(
        rf"Design-field distribution, "
        rf"$\Sigma_{{\rm gas}}={sigma_gas_target:.0f}\,M_\odot\,{{\rm pc}}^{{-2}}$",
        fontsize="medium", y=0.99,
    )
    save(fig, out_dir / "design_fields_selection.png", dpi)


# ── figure 2: corner plot with samples ────────────────────────────────────────

def plot_corner(
    reference,
    samples: dict[int, pd.DataFrame],
    fields: list[str],
    sigma_gas_target: float,
    label_suffix: str,
    colors: list[str],
    dpi: int,
    out_dir: Path,
    fname: str,
) -> None:
    """5×5 corner: PHANGS reference (gray) + Sobol samples (colored)."""
    n_fields = len(fields)
    fig, axes = _corner_axes(fields)

    log_ref = {f: _log(reference, f) for f in fields}
    ranges  = {f: _axis_range(log_ref[f]) for f in fields}
    n_bins  = 25

    for row, yf in enumerate(fields):
        for col, xf in enumerate(fields):
            ax = axes[row, col]
            if col > row:
                ax.set_visible(False)
                continue

            xlo, xhi = ranges[xf]
            bins_x = np.linspace(xlo, xhi, n_bins + 1)

            if col == row:
                ax.hist(log_ref[xf], bins=bins_x, density=True,
                        color=COLOR_ALL, edgecolor="none", label="PHANGS ref")
                for i, (n, sdf) in enumerate(sorted(samples.items())):
                    ax.hist(np.log10(sdf[xf].values), bins=bins_x, density=True,
                            histtype="step",
                            color=colors[i % len(colors)], linewidth=1.4,
                            label=f"$n={n}${label_suffix}")
                ax.set_yticks([])
            else:
                ylo, yhi = ranges[yf]
                ax.scatter(log_ref[xf], log_ref[yf],
                           s=1.2, c=COLOR_ALL, rasterized=True, zorder=1)
                for i, (n, sdf) in enumerate(sorted(samples.items())):
                    ms = max(5, 140 // n)
                    ax.scatter(np.log10(sdf[xf].values),
                               np.log10(sdf[yf].values),
                               s=ms, c=colors[i % len(colors)],
                               zorder=2 + i, alpha=0.85)
                ax.set_ylim(ylo, yhi)

            ax.set_xlim(xlo, xhi)

    _decorate_corner(axes, fields)
    axes[0, 0].legend(fontsize="small", loc="upper right", framealpha=0.85)
    fig.suptitle(
        rf"KDE-Sobol design, "
        rf"$\Sigma_{{\rm gas}}={sigma_gas_target:.0f}\,M_\odot\,{{\rm pc}}^{{-2}}$",
        fontsize="medium", y=0.99,
    )
    save(fig, out_dir / fname, dpi)


# ── figure 3: marginal distributions standard vs expanded ─────────────────────

def plot_distributions(
    reference,
    std_samples: dict[int, pd.DataFrame],
    exp_sample: pd.DataFrame | None,
    fields: list[str],
    sigma_gas_target: float,
    n_expanded: int,
    dpi: int,
    out_dir: Path,
) -> None:
    """1-D marginals: standard (solid) vs expanded (dashed) at matched n."""
    ncols = 3
    nrows = (len(fields) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.0 * ncols, 3.2 * nrows))
    n_bins = 30

    for fi, f in enumerate(fields):
        ax = list(axes.flat)[fi]
        log_ref = _log(reference, f)
        xlo, xhi = _axis_range(log_ref)
        bins = np.linspace(xlo, xhi, n_bins + 1)

        ax.hist(log_ref, bins=bins, density=True,
                color=COLOR_ALL, edgecolor="none", label="PHANGS ref")

        for i, (n, sdf) in enumerate(sorted(std_samples.items())):
            ax.hist(np.log10(sdf[f].values), bins=bins, density=True,
                    histtype="step",
                    color=COLOR_STD[i % len(COLOR_STD)],
                    linewidth=1.4, label=f"$n={n}$")

        if exp_sample is not None:
            ax.hist(np.log10(exp_sample[f].values), bins=bins, density=True,
                    histtype="step", color=COLOR_EXP, linewidth=1.6,
                    linestyle="--",
                    label=rf"$n={n_expanded}$ (expanded)")

        ax.set_xlabel(FIELD_LABELS[f], fontsize="medium")
        ax.set_ylabel("PDF", fontsize="medium")
        ax.tick_params(labelsize="small")
        ax.set_xlim(xlo, xhi)

    for ax in list(axes.flat)[len(fields):]:
        ax.set_visible(False)

    handles, labels = list(axes.flat)[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower right", fontsize="medium",
               framealpha=0.85, bbox_to_anchor=(1.0, 0.02))
    fig.suptitle(
        rf"Marginal distributions, "
        rf"$\Sigma_{{\rm gas}}={sigma_gas_target:.0f}\,M_\odot\,{{\rm pc}}^{{-2}}$",
        fontsize="large",
    )
    fig.tight_layout()
    save(fig, out_dir / "design_fields_distributions.png", dpi)


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    print("Loading PHANGS data ...")
    result = load_configured_phangs(str(args.config), base_dir=str(ROOT))
    table = result["table"]
    print(f"  {len(table)} apertures")

    cfg = SamplingConfig(
        delta_sigma_gas=args.delta,
        sobol_seed=args.seed,
        kde_aux_sample_size=100_000,
    )
    designer = PHANGSSamplingDesigner(table, config=cfg)
    fields = cfg.design_fields

    for sigma_gas_target in args.sigma_gas:
        print(f"\n=== Sigma_gas = {sigma_gas_target:.0f} ===")
        tag     = f"Sgas{sigma_gas_target:.0f}"
        out_dir = FIGURE_DIR / tag
        designer  = designer.with_config(target_sigma_gas=sigma_gas_target)
        reference = designer.select_reference_pixels()
        print(f"  {len(reference)} reference pixels")

        # ── figure 1: reduced correlation matrix ──────────────────────────
        print("  Plotting selection (reduced correlation matrix) ...")
        plot_selection(table, reference, fields,
                       sigma_gas_target, args.delta, args.dpi, out_dir)

        # ── standard samples ──────────────────────────────────────────────
        std_samples: dict[int, pd.DataFrame] = {}
        for n in sorted(args.n_samples):
            print(f"  Standard KDE-Sobol n={n} ...", end=" ", flush=True)
            std_samples[n] = designer.synthesize_kde_sobol(reference, n_samples=n)
            print("done")

        n_str = "_".join(str(n) for n in sorted(args.n_samples))
        print("  Plotting standard corner ...")
        plot_corner(reference, std_samples, fields,
                    sigma_gas_target, label_suffix="",
                    colors=COLOR_STD, dpi=args.dpi, out_dir=out_dir,
                    fname=f"design_fields_corner_n{n_str}.png")

        # ── expanded sample ───────────────────────────────────────────────
        print(f"  Expanded KDE-Sobol n={args.n_expanded} ...", end=" ", flush=True)
        exp_sample = designer.synthesize_expanded_kde_sobol(
            reference, n_samples=args.n_expanded
        )
        print(f"done (bw_factor={exp_sample.attrs['bandwidth_factor']})")

        print("  Plotting expanded corner ...")
        # Show standard n=args.n_expanded vs expanded at same n
        n_match = args.n_expanded
        cmp_samples = {
            n_match: std_samples.get(
                n_match,
                designer.synthesize_kde_sobol(reference, n_samples=n_match)
            ),
        }
        plot_corner(reference, cmp_samples, fields,
                    sigma_gas_target,
                    label_suffix=" (std)",
                    colors=COLOR_STD, dpi=args.dpi, out_dir=out_dir,
                    fname="design_fields_corner_expanded_cmp.png")

        # Corner with only the expanded sample
        plot_corner(reference, {n_match: exp_sample}, fields,
                    sigma_gas_target,
                    label_suffix=" (exp)",
                    colors=[COLOR_EXP], dpi=args.dpi, out_dir=out_dir,
                    fname="design_fields_corner_expanded.png")

        # ── figure 3: distributions standard vs expanded ──────────────────
        print("  Plotting distributions (standard + expanded) ...")
        plot_distributions(reference, std_samples, exp_sample,
                           fields, sigma_gas_target, args.n_expanded,
                           args.dpi, out_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
