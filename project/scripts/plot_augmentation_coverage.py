#!/usr/bin/env python
"""Plot augmented-prior coverage: narrow core vs augmentation block vs PHANGS.

Shows the five design-field marginals on axes spanning the *full* design range,
so the broadened H_star and Omega tails are visible extending beyond the PHANGS
reference (unlike plot_design_from_csv.py, whose axes are locked to the PHANGS
range and therefore clip the augmented points). Gray = PHANGS reference band,
blue = existing narrow design, red = augmentation block.

Usage:
    python project/scripts/plot_augmentation_coverage.py \
        project/output/design_Sgas10.0_n0032.csv \
        project/output/augment_design_Sgas10.0_n0032_n0032.csv
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
plt.style.use(ROOT / "project" / "scripts" / "prfm.mplstyle")
sys.path.insert(0, str(ROOT))

DESIGN_FIELDS = ["Sigma_gas", "Sigma_star", "H_star", "Omega", "qshear"]
AUGMENTED_DEFAULT = ["H_star", "Omega"]
FIELD_LABELS = {
    "Sigma_gas":  r"$\log\,\Sigma_{\rm gas}$ [$M_\odot\,{\rm pc}^{-2}$]",
    "Sigma_star": r"$\log\,\Sigma_\star$ [$M_\odot\,{\rm pc}^{-2}$]",
    "H_star":     r"$\log\,H_\star$ [pc]",
    "Omega":      r"$\log\,\Omega$ [km s$^{-1}$ kpc$^{-1}$]",
    "qshear":     r"$\log\,q$",
}
R8_FIDUCIAL = {
    "Sigma_gas": 12.0, "Sigma_star": 42.0, "H_star": 245.0,
    "Omega": 28.0, "qshear": 1.0,
}
COLOR_REF, COLOR_NARROW, COLOR_BLOCK, COLOR_R8 = "0.8", "#3a86ff", "#e63946", "#ffb000"


def build_coverage_figure(
    reference, narrow, block,
    fields: list[str] = DESIGN_FIELDS,
    augmented: list[str] = AUGMENTED_DEFAULT,
):
    """Return a 2x3 marginals figure with union-spanning axes.

    Each of ``reference``, ``narrow``, ``block`` is a mapping (DataFrame/Table)
    with the design-field columns. Axes span the union of all three so the
    augmented tails remain visible.
    """
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.2))
    axes = axes.ravel()
    for ax, f in zip(axes, fields):
        rlog = np.log10(np.asarray(reference[f], dtype=float))
        nlog = np.log10(np.asarray(narrow[f], dtype=float))
        blog = np.log10(np.asarray(block[f], dtype=float))
        lo = min(rlog.min(), nlog.min(), blog.min())
        hi = max(rlog.max(), nlog.max(), blog.max())
        pad = 0.03 * (hi - lo)
        bins = np.linspace(lo - pad, hi + pad, 26)
        ax.hist(rlog, bins=bins, density=True, color=COLOR_REF, label="PHANGS ref")
        ax.hist(nlog, bins=bins, density=True, histtype="step",
                color=COLOR_NARROW, lw=1.4, label="narrow (core)")
        ax.hist(blog, bins=bins, density=True, histtype="step",
                color=COLOR_BLOCK, lw=1.4, label="augment block")
        ax.axvline(np.log10(R8_FIDUCIAL[f]), color=COLOR_R8, ls="--", lw=1.0)
        label = FIELD_LABELS[f] + ("  (augmented)" if f in augmented else "")
        ax.set_xlabel(label, fontsize="small")
        ax.set_yticks([])
        ax.tick_params(labelsize="x-small")

    for ax in axes[len(fields):]:
        ax.axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[-1].legend(handles, labels, loc="center", fontsize="small", frameon=False)
    fig.tight_layout()
    return fig


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("narrow_csv", type=Path, help="Existing (narrow) design CSV")
    p.add_argument("block_csv", type=Path, help="Augmentation-block CSV")
    p.add_argument("--sigma-gas", type=float, default=None)
    p.add_argument("--delta", type=float, default=0.3)
    p.add_argument("--output", "-o", type=Path, default=None)
    p.add_argument("--dpi", type=int, default=150)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    # Reuse augment_design's PHANGS loader / band selection.
    spec = importlib.util.spec_from_file_location(
        "augment_design", Path(__file__).with_name("augment_design.py")
    )
    ad = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ad)

    narrow = pd.read_csv(args.narrow_csv)
    block = pd.read_csv(args.block_csv)
    sigma_gas = args.sigma_gas or ad.parse_sigma_gas_from_stem(args.narrow_csv.stem)

    _, reference = ad.select_reference(sigma_gas, args.delta, ad.CONFIG_PATH, 100_000)
    reference = {f: np.asarray(reference[f], dtype=float) for f in DESIGN_FIELDS}

    fig = build_coverage_figure(reference, narrow, block)
    out = args.output or (ROOT / "figures" / "phangs" / "sobol" / f"Sgas{sigma_gas}"
                          / "augmentation_coverage.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
