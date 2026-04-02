"""
Histogram exploration of PHANGS megatable key fields.

Produces three figures saved under figures/phangs/:
  field_distributions.png       -- overall distribution of each field (all data)
  field_distributions_bygal.png -- per-galaxy distributions overlaid
  field_distributions_byrgal.png -- distributions split by galactocentric radius bin

Usage
-----
python scripts/phangs_explore_fields.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from prfm import phangs

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Key columns to examine: (column_name, x_label, log_scale)
FIELDS: list[tuple[str, str, bool]] = [
    ("r_gal",               r"$r_\mathrm{gal}$ [kpc]",                    False),
    ("Sigma_mol",           r"$\Sigma_\mathrm{mol}$ [$M_\odot$ pc$^{-2}$]", True),
    ("Sigma_atom",          r"$\Sigma_\mathrm{atom}$ [$M_\odot$ pc$^{-2}$]", True),
    ("Sigma_gas",           r"$\Sigma_\mathrm{gas}$ [$M_\odot$ pc$^{-2}$]", True),
    ("Sigma_star",          r"$\Sigma_\star$ [$M_\odot$ pc$^{-2}$]",       True),
    ("rho_star_mp",         r"$\rho_\star$ [$M_\odot$ pc$^{-3}$]",         True),
    ("H_star",              r"$H_\star$ [pc]",                              True),
    ("V_circ_CO21_URC",     r"$V_\mathrm{circ}$ [km s$^{-1}$]",            False),
    ("Omega_d",             r"$\Omega_d$ [km s$^{-1}$ kpc$^{-1}$]",        True),
    ("Zprime",              r"$Z'$ [$Z_\odot$]",                            False),
    ("Sigma_SFR_HaW4recal", r"$\Sigma_\mathrm{SFR}^\mathrm{Ha+W4}$ [$M_\odot$ yr$^{-1}$ kpc$^{-2}$]", True),
    ("Sigma_SFR_FUVW4recal",r"$\Sigma_\mathrm{SFR}^\mathrm{FUV+W4}$ [$M_\odot$ yr$^{-1}$ kpc$^{-2}$]", True),
    ("alpha_CO21_S20",      r"$\alpha_\mathrm{CO}$ [S20]",                  True),
]

# Radial bins for the by-r_gal figure
RGAL_BINS = [(0, 2), (2, 5), (5, 10), (10, 999)]
RGAL_LABELS = [r"$r<2$ kpc", r"$2-5$ kpc", r"$5-10$ kpc", r"$r>10$ kpc"]
RGAL_COLORS = ["#c0392b", "#e67e22", "#27ae60", "#2980b9"]

REPO_ROOT = Path(__file__).parent.parent
DATA_DIR = REPO_ROOT / "data/phangs_megatable"
OUT_DIR = REPO_ROOT / "figures/phangs"
N_BINS = 40

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def clean(vals: np.ndarray, log: bool) -> np.ndarray:
    """Return finite, positive (if log) values as a plain float array."""
    v = np.asarray(vals, dtype=float)
    mask = np.isfinite(v)
    if log:
        mask &= v > 0
    return v[mask]


def bin_edges(vals: np.ndarray, log: bool, n: int = N_BINS) -> np.ndarray:
    if log:
        return np.logspace(np.log10(vals.min()), np.log10(vals.max()), n + 1)
    return np.linspace(vals.min(), vals.max(), n + 1)


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

print("Loading PHANGS annulus data…")
t = phangs.load_all(DATA_DIR, aperture="annulus")
t = phangs.compute_prfm_inputs(t)

galaxies = sorted(set(t["GALAXY"]))
gal_cmap = plt.cm.tab20
gal_colors = {g: gal_cmap(i / len(galaxies)) for i, g in enumerate(galaxies)}
r_gal_all = np.asarray(t["r_gal"], dtype=float)

OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Figure 1: overall distributions
# ---------------------------------------------------------------------------

ncols = 4
nrows = int(np.ceil(len(FIELDS) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 2.8))
axes = axes.flatten()

for ax, (col, xlabel, log) in zip(axes, FIELDS):
    v = clean(np.asarray(t[col], dtype=float), log)
    if len(v) == 0:
        ax.set_visible(False)
        continue
    edges = bin_edges(v, log)
    ax.hist(v, bins=edges, color="steelblue", edgecolor="none", alpha=0.85, density=True)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel("PDF", fontsize=8)
    if log:
        ax.set_xscale("log")
    nan_pct = 100.0 * (1 - len(v) / len(t))
    ax.set_title(f"{col}  (NaN {nan_pct:.0f}%)", fontsize=7, pad=3)
    ax.tick_params(labelsize=7)

for ax in axes[len(FIELDS):]:
    ax.set_visible(False)

fig.suptitle("PHANGS megatable — key field distributions (all data)", fontsize=11, y=1.01)
fig.tight_layout()
out = OUT_DIR / "field_distributions.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {out}")

# ---------------------------------------------------------------------------
# Figure 2: per-galaxy distributions
# ---------------------------------------------------------------------------

fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 2.8))
axes = axes.flatten()

for ax, (col, xlabel, log) in zip(axes, FIELDS):
    v_all = clean(np.asarray(t[col], dtype=float), log)
    if len(v_all) == 0:
        ax.set_visible(False)
        continue
    edges = bin_edges(v_all, log)
    # per-galaxy thin lines
    for gal in galaxies:
        mask = np.asarray(t["GALAXY"]) == gal
        v_g = clean(np.asarray(t[col][mask], dtype=float), log)
        if len(v_g) < 3:
            continue
        counts, _ = np.histogram(v_g, bins=edges, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        ax.plot(centers, counts, lw=0.6, alpha=0.5, color=gal_colors[gal])
    # overall thick line
    counts_all, _ = np.histogram(v_all, bins=edges, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    ax.plot(centers, counts_all, lw=2, color="black", label="all")
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel("PDF", fontsize=8)
    if log:
        ax.set_xscale("log")
    ax.set_title(col, fontsize=7, pad=3)
    ax.tick_params(labelsize=7)

for ax in axes[len(FIELDS):]:
    ax.set_visible(False)

fig.suptitle("PHANGS — per-galaxy distributions (thin) vs. all (black)", fontsize=11, y=1.01)
fig.tight_layout()
out = OUT_DIR / "field_distributions_bygal.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {out}")

# ---------------------------------------------------------------------------
# Figure 3: distributions split by r_gal bin
# ---------------------------------------------------------------------------

fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 2.8))
axes = axes.flatten()

for ax, (col, xlabel, log) in zip(axes, FIELDS):
    v_all = clean(np.asarray(t[col], dtype=float), log)
    if len(v_all) == 0:
        ax.set_visible(False)
        continue
    edges = bin_edges(v_all, log)
    centers = 0.5 * (edges[:-1] + edges[1:])
    for (rmin, rmax), label, color in zip(RGAL_BINS, RGAL_LABELS, RGAL_COLORS):
        rmask = (r_gal_all >= rmin) & (r_gal_all < rmax)
        v_r = clean(np.asarray(t[col][rmask], dtype=float), log)
        if len(v_r) < 3:
            continue
        counts, _ = np.histogram(v_r, bins=edges, density=True)
        ax.plot(centers, counts, lw=1.5, color=color, label=label)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel("PDF", fontsize=8)
    if log:
        ax.set_xscale("log")
    ax.set_title(col, fontsize=7, pad=3)
    ax.tick_params(labelsize=7)

# shared legend on first visible axis
axes[0].legend(fontsize=6, loc="upper right")

for ax in axes[len(FIELDS):]:
    ax.set_visible(False)

fig.suptitle(r"PHANGS — distributions split by $r_\mathrm{gal}$", fontsize=11, y=1.01)
fig.tight_layout()
out = OUT_DIR / "field_distributions_byrgal.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {out}")

print("Done.")
