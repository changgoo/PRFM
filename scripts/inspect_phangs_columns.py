"""
Quick inspection of PHANGS megatable columns, dtypes, units, and NaN coverage.
Output is printed to stdout.

Usage
-----
python scripts/inspect_phangs_columns.py [--aperture annulus|gauss|hexagon]
"""

import argparse
from pathlib import Path

import numpy as np
from prfm import phangs

REPO_ROOT = Path(__file__).parent.parent

parser = argparse.ArgumentParser()
parser.add_argument("--aperture", default="annulus",
                    choices=["annulus", "gauss", "hexagon"])
args = parser.parse_args()

t = phangs.load_all(REPO_ROOT / "data/phangs_megatable", aperture=args.aperture)
print(f"Aperture   : {args.aperture}")
t = phangs.compute_prfm_inputs(t)

print(f"Total rows : {len(t)}")
print(f"Galaxies   : {len(set(t['GALAXY']))}")
print()
print(f"{'Column':<40} {'Unit':<25} {'dtype':<10} {'NaN%':>6}  {'min':>12}  {'max':>12}")
print("-" * 115)
for col in t.colnames:
    c = t[col]
    unit = str(getattr(c, "unit", "")) or ""
    dtype = str(c.dtype)
    try:
        vals = np.asarray(c, dtype=float)
        nan_pct = 100.0 * np.sum(~np.isfinite(vals)) / len(vals)
        finite = vals[np.isfinite(vals)]
        vmin = f"{finite.min():.3g}" if len(finite) else "—"
        vmax = f"{finite.max():.3g}" if len(finite) else "—"
    except (TypeError, ValueError):
        nan_pct = float("nan")
        vmin = vmax = "—"
    print(f"{col:<40} {unit:<25} {dtype:<10} {nan_pct:>6.1f}%  {vmin:>12}  {vmax:>12}")
