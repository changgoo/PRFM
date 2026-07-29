#!/usr/bin/env python
"""Build a one-at-a-time (OAT) physics-ablation design CSV.

Starting from an existing physical design CSV (the columns
    Sigma_gas, Sigma_star, H_star, Omega, qshear
produced by run_sampling.py), this generates an ablation grid that sweeps the
physics parameters
    Z_gas, f_dtm (dust-to-metal ratio), xi_CR_amp
one at a time about their fiducial value (1.0), at a handful of physical anchor
points drawn from the design. Dust metallicity is derived as
    Z_dust = f_dtm * Z_gas
so the constraint Z_dust <= Z_gas is automatic.

Each swept row holds two of the three physics parameters at the base level and
steps the third through the sub-fiducial levels; the base (1, 1, 1) corner is
NOT emitted because those runs already exist in the original physical suite.
The resulting CSV carries a descriptive `suffix` column (consumed by
csv_to_slurm_yaml.py) plus provenance columns.

See project/doc/physics-parameter-extension.md for the design rationale.

Usage:
    python project/scripts/make_ablation_grid.py \
        project/output/design_Sgas10.0_n0032.csv
    python project/scripts/make_ablation_grid.py <design.csv> \
        --anchors 7,15,23 --levels 0.1,0.22,0.46,1.0 -o <out.csv>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

PHYSICAL_FIELDS = ["Sigma_gas", "Sigma_star", "H_star", "Omega", "qshear"]

# Sweep parameter -> short tag used in the job-id suffix.
SWEEP_TAGS = {"Z_gas": "Z", "f_dtm": "fd", "xi_CR_amp": "cr"}

DEFAULT_LEVELS = [0.1, 0.22, 0.46, 1.0]
DEFAULT_BASE_LEVEL = 1.0


def _fmt_level(level: float) -> str:
    """Format a level for a job-id suffix, e.g. 0.1 -> '0p10'."""
    return f"{level:.2f}".replace(".", "p")


def select_anchor_indices(
    df: pd.DataFrame, n_anchors: int, by: str = "Sigma_star"
) -> list[int]:
    """Pick ``n_anchors`` original row indices spanning the design in ``by``.

    Rows are ordered by ``by`` and sampled at evenly spaced quantiles, so
    n_anchors=3 returns the minimum, median, and maximum rows. The returned
    indices are the *original* DataFrame indices (i.e. the row#### ids of the
    physical suite), so the anchors coincide with existing simulations.
    """
    if by not in df.columns:
        raise ValueError(f"Cannot span anchors by missing column: {by!r}")
    if n_anchors < 1:
        raise ValueError("n_anchors must be >= 1")
    order = df[by].to_numpy().argsort(kind="mergesort")
    positions = np.linspace(0, len(order) - 1, n_anchors).round().astype(int)
    # Deduplicate while preserving order (small designs / large n_anchors).
    seen: set[int] = set()
    picked: list[int] = []
    for p in positions:
        idx = int(df.index[order[p]])
        if idx not in seen:
            seen.add(idx)
            picked.append(idx)
    return picked


def build_ablation_grid(
    df: pd.DataFrame,
    anchor_indices: list[int],
    sweep_params: list[str] = list(SWEEP_TAGS),
    levels: list[float] = DEFAULT_LEVELS,
    base_level: float = DEFAULT_BASE_LEVEL,
    anchor_width: int = 4,
) -> pd.DataFrame:
    """Return the OAT ablation design as a DataFrame.

    For each anchor and each sweep parameter, one row is produced per level
    that differs from ``base_level``; the other two physics parameters are held
    at ``base_level``. Physics columns emitted are ``Z_gas``, ``Z_dust``
    (= f_dtm * Z_gas), and ``xi_CR_amp``.
    """
    missing = [c for c in PHYSICAL_FIELDS if c not in df.columns]
    if missing:
        raise ValueError(f"Design CSV missing physical columns: {missing}")

    sub_levels = [lvl for lvl in levels if lvl != base_level]
    rows: list[dict] = []
    for anchor in anchor_indices:
        phys = {f: float(df.at[anchor, f]) for f in PHYSICAL_FIELDS}
        anchor_suffix = f"row{anchor:0{anchor_width}d}"
        for param in sweep_params:
            for level in sub_levels:
                # Base values for the three physics knobs, then override one.
                z_gas = base_level
                f_dtm = base_level
                xi_cr = base_level
                if param == "Z_gas":
                    z_gas = level
                elif param == "f_dtm":
                    f_dtm = level
                elif param == "xi_CR_amp":
                    xi_cr = level
                else:
                    raise ValueError(f"Unknown sweep parameter: {param!r}")

                rows.append(
                    {
                        **phys,
                        "Z_gas": z_gas,
                        "Z_dust": f_dtm * z_gas,
                        "xi_CR_amp": xi_cr,
                        "suffix": f"a{anchor:0{anchor_width}d}_{SWEEP_TAGS[param]}{_fmt_level(level)}",
                        # provenance (ignored by csv_to_slurm_yaml)
                        "anchor_row": anchor,
                        "anchor_suffix": anchor_suffix,
                        "sweep_param": param,
                        "level": level,
                        "f_dtm": f_dtm,
                    }
                )

    return pd.DataFrame(rows)


def default_output_path(design_csv: Path) -> Path:
    return ROOT / "project" / "output" / f"ablation_{design_csv.stem}.csv"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("design_csv", type=Path, help="Physical design CSV")
    p.add_argument(
        "--anchors",
        default=None,
        help="Comma-separated original row indices to use as anchors. If"
        " omitted, --n-anchors rows spanning --by are chosen.",
    )
    p.add_argument("--n-anchors", type=int, default=3, help="Number of anchors (default 3)")
    p.add_argument(
        "--by", default="Sigma_star", help="Column to span anchors by (default Sigma_star)"
    )
    p.add_argument(
        "--levels",
        default=",".join(str(x) for x in DEFAULT_LEVELS),
        help="Comma-separated sweep levels (default 0.1,0.22,0.46,1.0)",
    )
    p.add_argument(
        "--base-level", type=float, default=DEFAULT_BASE_LEVEL,
        help="Fiducial level held for un-swept parameters (default 1.0)",
    )
    p.add_argument(
        "--sweep",
        default=",".join(SWEEP_TAGS),
        help="Comma-separated sweep parameters (default Z_gas,f_dtm,xi_CR_amp)",
    )
    p.add_argument("--output", "-o", type=Path, default=None, help="Output CSV path")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if not args.design_csv.exists():
        print(f"ERROR: design CSV not found: {args.design_csv}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(args.design_csv)
    levels = [float(x) for x in args.levels.split(",")]
    sweep_params = [s.strip() for s in args.sweep.split(",")]

    if args.anchors is not None:
        anchor_indices = [int(x) for x in args.anchors.split(",")]
    else:
        anchor_indices = select_anchor_indices(df, args.n_anchors, by=args.by)

    grid = build_ablation_grid(
        df, anchor_indices, sweep_params=sweep_params,
        levels=levels, base_level=args.base_level,
    )

    out_path = args.output or default_output_path(args.design_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    grid.to_csv(out_path, index=False)

    print(f"Anchors (row indices): {anchor_indices}")
    print(f"Sweep params: {sweep_params}")
    print(f"Sub-fiducial levels: {[lvl for lvl in levels if lvl != args.base_level]}")
    print(f"Wrote {len(grid)} ablation rows → {out_path}")


if __name__ == "__main__":
    main()
