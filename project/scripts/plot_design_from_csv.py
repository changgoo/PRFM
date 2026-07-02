#!/usr/bin/env python
"""Plot a design CSV overlaid on the PHANGS reference distribution.

Unlike plot_design_fields.py (which re-samples internally), this script reads
the exact CSV that will be handed to TIGRESS-NCR and shows it in context.

Outputs (next to the CSV, unless --output-dir is given):
    <stem>_selection.png     PHANGS full distribution + Sigma_gas band highlight
                              + design points overlaid
    <stem>_corner.png         5x5 corner: PHANGS reference (gray) + design (red)
    <stem>_marginals.png      1-D marginals: PHANGS reference (gray filled) +
                              design (red step)

Usage:
    python project/scripts/plot_design_from_csv.py \\
        project/output/design_Sgas10.0_n0032.csv

    python project/scripts/plot_design_from_csv.py \\
        project/output/design_Sgas10.0_n0032.csv \\
        --sigma-gas 10 --delta 0.3 \\
        --output-dir figures/phangs/sobol/Sgas10/pilot_n32
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
plt.style.use(ROOT / "project" / "scripts" / "prfm.mplstyle")

from prfm.phangs import load_configured_phangs
from prfm.phangs_sampling import PHANGSSamplingDesigner, SamplingConfig

CONFIG_PATH = ROOT / "config" / "phangs_prfm.yml"

DESIGN_FIELDS = ["Sigma_gas", "Sigma_star", "H_star", "Omega", "qshear"]

FIELD_LABELS = {
    "Sigma_gas":  r"$\log\,\Sigma_{\rm gas}$"  + "\n" + r"[$M_\odot\,{\rm pc}^{-2}$]",
    "Sigma_star": r"$\log\,\Sigma_\star$"       + "\n" + r"[$M_\odot\,{\rm pc}^{-2}$]",
    "H_star":     r"$\log\,H_\star$"            + "\n" + r"[pc]",
    "Omega":      r"$\log\,\Omega$"             + "\n" + r"[km s$^{-1}$ kpc$^{-1}$]",
    "qshear":     r"$\log\,q$",
}

COLOR_ALL    = "0.78"
COLOR_REF    = "#3a86ff"
COLOR_DESIGN = "#e63946"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("csv", type=Path, help="Design CSV to plot")
    p.add_argument("--sigma-gas", type=float, default=None,
                   help="Target Sigma_gas for band highlight (default: parsed from CSV)")
    p.add_argument("--delta", type=float, default=0.3,
                   help="Half-width of log10 Sigma_gas band in dex (default: 0.3)")
    p.add_argument("--config", type=Path, default=CONFIG_PATH,
                   help="PHANGS config YAML")
    p.add_argument("--output-dir", type=Path, default=None,
                   help="Where to save figures (default: alongside the CSV)")
    p.add_argument("--dpi", type=int, default=150)
    return p.parse_args()


def _log(vals) -> np.ndarray:
    arr = np.asarray(vals, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        lv = np.log10(arr)
    lv[~np.isfinite(lv)] = np.nan
    return lv


def _range(log_vals: np.ndarray, pad: float = 0.05) -> tuple[float, float]:
    lo, hi = np.nanpercentile(log_vals, [0.5, 99.5])
    span = hi - lo
    return lo - pad * span, hi + pad * span


def _corner_axes(size: float = 2.2):
    n = len(DESIGN_FIELDS)
    fig, axes = plt.subplots(n, n, figsize=(size * n, size * n))
    fig.subplots_adjust(hspace=0.08, wspace=0.08,
                        left=0.12, bottom=0.10, right=0.98, top=0.95)
    return fig, axes


def _decorate(axes: np.ndarray) -> None:
    n = len(DESIGN_FIELDS)
    for row in range(n):
        for col in range(n):
            ax = axes[row, col]
            if col > row:
                ax.set_visible(False)
                continue
            ax.tick_params(labelsize="small")
            if row == n - 1:
                ax.set_xlabel(FIELD_LABELS[DESIGN_FIELDS[col]], fontsize="small")
                ax.xaxis.set_major_locator(ticker.MaxNLocator(3, prune="both"))
            else:
                ax.tick_params(labelbottom=False)
            if col == 0 and row > 0:
                ax.set_ylabel(FIELD_LABELS[DESIGN_FIELDS[row]], fontsize="small")
                ax.yaxis.set_major_locator(ticker.MaxNLocator(3, prune="both"))
            else:
                ax.tick_params(labelleft=False)


def _save(fig, path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    print(f"  saved -> {path}")


def _title(sigma_gas: float, delta: float, n: int) -> str:
    return (rf"$\Sigma_{{\rm gas}}={sigma_gas:.0f}\,M_\odot\,{{\rm pc}}^{{-2}}$"
            rf" ($\Delta={delta}\,$dex, $n={n}$)")


def plot_selection(full_table, reference, design: pd.DataFrame,
                   sigma_gas: float, delta: float, dpi: int, out: Path) -> None:
    """Corner plot: full PHANGS + selected band + design points."""
    fig, axes = _corner_axes()
    log_full = {f: _log(full_table[f]) for f in DESIGN_FIELDS}
    log_ref  = {f: _log(reference[f])  for f in DESIGN_FIELDS}
    log_des  = {f: _log(design[f])     for f in DESIGN_FIELDS}
    ranges   = {f: _range(log_full[f]) for f in DESIGN_FIELDS}
    n_bins = 30

    for row, yf in enumerate(DESIGN_FIELDS):
        for col, xf in enumerate(DESIGN_FIELDS):
            ax = axes[row, col]
            if col > row:
                ax.set_visible(False)
                continue
            xlo, xhi = ranges[xf]
            bins_x = np.linspace(xlo, xhi, n_bins + 1)
            if col == row:
                ax.hist(log_full[xf], bins=bins_x, density=True,
                        color=COLOR_ALL, edgecolor="none")
                ax.hist(log_ref[xf], bins=bins_x, density=True,
                        histtype="step", color=COLOR_REF, linewidth=1.5)
                ax.hist(log_des[xf], bins=bins_x, density=True,
                        histtype="step", color=COLOR_DESIGN, linewidth=1.5)
                ax.set_yticks([])
            else:
                ylo, yhi = ranges[yf]
                ax.scatter(log_full[xf], log_full[yf],
                           s=0.8, c=COLOR_ALL, rasterized=True, zorder=1)
                ax.scatter(log_ref[xf], log_ref[yf],
                           s=1.5, c=COLOR_REF, rasterized=True, zorder=2, alpha=0.7)
                ax.scatter(log_des[xf], log_des[yf],
                           s=18, c=COLOR_DESIGN, edgecolor="white", linewidth=0.4,
                           zorder=3)
                ax.set_ylim(ylo, yhi)
            ax.set_xlim(xlo, xhi)

    _decorate(axes)
    handles = [
        Patch(facecolor=COLOR_ALL, label=f"PHANGS ({len(full_table):,})"),
        Line2D([0], [0], color=COLOR_REF, linewidth=2,
               label=(rf"$|\log\Sigma_{{\rm gas}}-\log{sigma_gas:.0f}|"
                      rf"\leq{delta}\,$dex ({len(reference):,})")),
        Line2D([0], [0], color=COLOR_DESIGN, marker="o", linestyle="none",
               markersize=6, markeredgecolor="white",
               label=f"design ({len(design):,})"),
    ]
    axes[0, 0].legend(handles=handles, fontsize="small",
                      loc="upper right", framealpha=0.85)
    fig.suptitle(_title(sigma_gas, delta, len(design)),
                 fontsize="medium", y=0.99)
    _save(fig, out, dpi)


def plot_corner(reference, design: pd.DataFrame,
                sigma_gas: float, delta: float, dpi: int, out: Path) -> None:
    """Reference (gray) + design (red). Zoomed to reference range."""
    fig, axes = _corner_axes()
    log_ref = {f: _log(reference[f]) for f in DESIGN_FIELDS}
    log_des = {f: _log(design[f])    for f in DESIGN_FIELDS}
    ranges  = {f: _range(log_ref[f]) for f in DESIGN_FIELDS}
    n_bins = 25

    for row, yf in enumerate(DESIGN_FIELDS):
        for col, xf in enumerate(DESIGN_FIELDS):
            ax = axes[row, col]
            if col > row:
                ax.set_visible(False)
                continue
            xlo, xhi = ranges[xf]
            bins_x = np.linspace(xlo, xhi, n_bins + 1)
            if col == row:
                ax.hist(log_ref[xf], bins=bins_x, density=True,
                        color=COLOR_ALL, edgecolor="none")
                ax.hist(log_des[xf], bins=bins_x, density=True,
                        histtype="step", color=COLOR_DESIGN, linewidth=1.5)
                ax.set_yticks([])
            else:
                ylo, yhi = ranges[yf]
                ax.scatter(log_ref[xf], log_ref[yf],
                           s=1.2, c=COLOR_ALL, rasterized=True, zorder=1)
                ax.scatter(log_des[xf], log_des[yf],
                           s=22, c=COLOR_DESIGN, edgecolor="white",
                           linewidth=0.5, zorder=2)
                ax.set_ylim(ylo, yhi)
            ax.set_xlim(xlo, xhi)

    _decorate(axes)
    handles = [
        Patch(facecolor=COLOR_ALL, label=f"PHANGS ref ({len(reference):,})"),
        Line2D([0], [0], color=COLOR_DESIGN, marker="o", linestyle="none",
               markersize=6, markeredgecolor="white",
               label=f"design ({len(design):,})"),
    ]
    axes[0, 0].legend(handles=handles, fontsize="small",
                      loc="upper right", framealpha=0.85)
    fig.suptitle(_title(sigma_gas, delta, len(design)),
                 fontsize="medium", y=0.99)
    _save(fig, out, dpi)


def plot_marginals(reference, design: pd.DataFrame,
                   sigma_gas: float, delta: float, dpi: int, out: Path) -> None:
    """1-D marginals for the 5 design fields."""
    ncols = 3
    nrows = (len(DESIGN_FIELDS) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.2 * nrows))
    n_bins = 30

    for fi, f in enumerate(DESIGN_FIELDS):
        ax = list(axes.flat)[fi]
        log_ref = _log(reference[f])
        log_des = _log(design[f])
        xlo, xhi = _range(log_ref)
        bins = np.linspace(xlo, xhi, n_bins + 1)
        ax.hist(log_ref, bins=bins, density=True,
                color=COLOR_ALL, edgecolor="none", label="PHANGS ref")
        ax.hist(log_des, bins=bins, density=True,
                histtype="step", color=COLOR_DESIGN,
                linewidth=1.5, label="design")
        ax.set_xlabel(FIELD_LABELS[f], fontsize="medium")
        ax.set_ylabel("PDF", fontsize="medium")
        ax.set_xlim(xlo, xhi)

    for ax in list(axes.flat)[len(DESIGN_FIELDS):]:
        ax.set_visible(False)

    handles, labels = list(axes.flat)[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower right",
               fontsize="medium", framealpha=0.85, bbox_to_anchor=(1.0, 0.02))
    fig.suptitle(_title(sigma_gas, delta, len(design)), fontsize="large")
    fig.tight_layout()
    _save(fig, out, dpi)


def _parse_sigma_gas_from_stem(stem: str) -> float | None:
    # Design CSVs are named like: design_Sgas10.0_n0032
    import re
    m = re.search(r"Sgas([0-9.]+)", stem)
    return float(m.group(1)) if m else None


def main() -> None:
    args = parse_args()
    if not args.csv.exists():
        print(f"ERROR: CSV not found: {args.csv}", file=sys.stderr)
        sys.exit(1)

    design = pd.read_csv(args.csv)
    missing = [c for c in DESIGN_FIELDS if c not in design.columns]
    if missing:
        print(f"ERROR: CSV missing columns: {missing}", file=sys.stderr)
        sys.exit(1)

    sigma_gas = args.sigma_gas
    if sigma_gas is None:
        sigma_gas = _parse_sigma_gas_from_stem(args.csv.stem)
        if sigma_gas is None:
            print("ERROR: --sigma-gas not given and could not parse from CSV name",
                  file=sys.stderr)
            sys.exit(1)

    print(f"Loading PHANGS data ...")
    result = load_configured_phangs(str(args.config), base_dir=str(ROOT))
    full_table = result["table"]
    print(f"  {len(full_table)} apertures loaded")

    designer = PHANGSSamplingDesigner(
        full_table,
        config=SamplingConfig(target_sigma_gas=sigma_gas,
                              delta_sigma_gas=args.delta),
    )
    reference = designer.select_reference_pixels()
    print(f"  {len(reference)} reference pixels in Sigma_gas={sigma_gas} band")
    print(f"  {len(design)} design points to plot")

    out_dir = args.output_dir or args.csv.parent
    stem = args.csv.stem

    plot_selection(full_table, reference, design, sigma_gas, args.delta,
                   args.dpi, out_dir / f"{stem}_selection.png")
    plot_corner(reference, design, sigma_gas, args.delta,
                args.dpi, out_dir / f"{stem}_corner.png")
    plot_marginals(reference, design, sigma_gas, args.delta,
                   args.dpi, out_dir / f"{stem}_marginals.png")

    print("Done.")


if __name__ == "__main__":
    main()
