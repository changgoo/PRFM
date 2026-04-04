"""
Correlation exploration of PHANGS megatable key fields.

Produces two figures saved under figures/phangs/<aperture>/:
  correlation_matrix.png   -- Spearman rank correlation heatmap
  scatter_matrix.png       -- scatter matrix colored by r_gal

Usage
-----
python scripts/phangs_explore_correlations.py [--aperture annulus|gauss|hexagon]
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from scipy import stats

from prfm import phangs
from prfm.phangs_plot import build_matrix, resolve_columns, to_log_or_raw

parser = argparse.ArgumentParser()
parser.add_argument("--aperture", default="annulus",
                    choices=["annulus", "gauss", "hexagon"])
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Base columns (suffix applied at load time): (base_name, x_label, log_scale)
_BASE_COLS: list[tuple[str, str, bool]] = [
    ("Sigma_gas",           r"$\Sigma_\mathrm{gas}$",  True),
    ("Sigma_mol",           r"$\Sigma_\mathrm{mol}$",  True),
    ("Sigma_atom",          r"$\Sigma_\mathrm{atom}$", True),
    ("Sigma_star",          r"$\Sigma_\star$",          True),
    ("H_star",              r"$H_\star$",               True),
    ("Omega_d",             r"$\Omega_d$",              True),
    ("Zprime",              r"$Z'$",                    False),
    ("Sigma_SFR_HaW4recal", r"$\Sigma_\mathrm{SFR}$",  True),
]

REPO_ROOT = Path(__file__).parent.parent
DATA_DIR = REPO_ROOT / "data/phangs_megatable"
OUT_DIR = REPO_ROOT / "figures/phangs" / args.aperture

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

print(f"Loading PHANGS {args.aperture} data…")
t = phangs.load_all(DATA_DIR, aperture=args.aperture)
t = phangs.vstack_tables(t)
if args.aperture == "annulus":
    t = phangs.compute_prfm_inputs(t)

# Resolve actual column names, keep only those present in this aperture's table
COLS = resolve_columns(_BASE_COLS, t, aperture=args.aperture)

r_gal_all = np.asarray(t["r_gal"], dtype=float)

OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Figure 1: Spearman correlation heatmap
# ---------------------------------------------------------------------------

# Compute pairwise Spearman correlations on all valid rows per pair
n = len(COLS)
rho_mat = np.zeros((n, n))
pval_mat = np.zeros((n, n))

for i, (ci, _, logi) in enumerate(COLS):
    vi = to_log_or_raw(t, ci, logi)
    for j, (cj, _, logj) in enumerate(COLS):
        vj = to_log_or_raw(t, cj, logj)
        valid = np.isfinite(vi) & np.isfinite(vj)
        if valid.sum() >= 10:
            r, p = stats.spearmanr(vi[valid], vj[valid])
        else:
            r, p = np.nan, np.nan
        rho_mat[i, j] = r
        pval_mat[i, j] = p

labels = [lbl for _, lbl, _ in COLS]

fig, ax = plt.subplots(figsize=(8, 7))
im = ax.imshow(rho_mat, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label(r"Spearman $\rho$", fontsize=10)

ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
ax.set_yticklabels(labels, fontsize=9)

# Annotate with rho values; mark significant ones
for i in range(n):
    for j in range(n):
        r = rho_mat[i, j]
        if not np.isfinite(r):
            continue
        sig = "**" if pval_mat[i, j] < 0.01 else ("*" if pval_mat[i, j] < 0.05 else "")
        color = "white" if abs(r) > 0.6 else "black"
        ax.text(j, i, f"{r:.2f}{sig}", ha="center", va="center",
                fontsize=7, color=color)

ax.set_title("Spearman rank correlation — PHANGS annulus\n"
             "(log-transformed where marked; * p<0.05, ** p<0.01)", fontsize=10)
fig.tight_layout()
out = OUT_DIR / "correlation_matrix.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {out}")

# ---------------------------------------------------------------------------
# Figure 2: scatter matrix colored by r_gal
# ---------------------------------------------------------------------------

data, labels = build_matrix(t, COLS)
n_valid = data.shape[0]

# Get r_gal for the same valid rows
arrays_tmp = [to_log_or_raw(t, col, log) for col, _, log in COLS]
stacked_tmp = np.column_stack(arrays_tmp)
valid_mask = np.all(np.isfinite(stacked_tmp), axis=1)
r_valid = r_gal_all[valid_mask]

norm = mcolors.LogNorm(vmin=max(r_valid.min(), 0.1), vmax=r_valid.max())
cmap = plt.cm.plasma_r

n = len(COLS)
fig, axes = plt.subplots(n, n, figsize=(n * 2.0, n * 2.0))

for i in range(n):
    for j in range(n):
        ax = axes[i, j]
        if i == j:
            # diagonal: histogram of column i
            ax.hist(data[:, i], bins=30, color="gray", alpha=0.7, density=True)
            ax.set_yticks([])
        else:
            sc = ax.scatter(data[:, j], data[:, i],
                            c=r_valid, cmap=cmap, norm=norm,
                            s=1, alpha=0.4, rasterized=True, linewidths=0)
        # axis labels on edges only
        if i == n - 1:
            ax.set_xlabel(labels[j], fontsize=7)
        else:
            ax.set_xticklabels([])
        if j == 0:
            ax.set_ylabel(labels[i], fontsize=7)
        else:
            ax.set_yticklabels([])
        ax.tick_params(labelsize=6)

# Colorbar on the right
fig.subplots_adjust(right=0.88, hspace=0.05, wspace=0.05)
cax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = fig.colorbar(sm, cax=cax)
cb.set_label(r"$r_\mathrm{gal}$ [kpc]", fontsize=9)

fig.suptitle(f"PHANGS scatter matrix — {n_valid} valid rows, colored by $r_{{\\rm gal}}$",
             fontsize=10, y=0.995)
out = OUT_DIR / "scatter_matrix.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {out}")

print("Done.")
