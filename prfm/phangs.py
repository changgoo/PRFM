"""
PHANGS megatable data handling for PRFM theory applications.

Data source: Sun et al. 2022/2023 (v4.0 public release)
  https://www.canfar.net/storage/vault/list/phangs/RELEASES/Sun_etal_2022

Key references:
  Sun et al. 2022 — https://ui.adsabs.harvard.edu/abs/2022AJ....164...43S
  Sun et al. 2023 — https://ui.adsabs.harvard.edu/abs/2023ApJ...945L..19S
"""

import warnings
from pathlib import Path
from typing import Optional

import astropy.constants as ac
import astropy.units as au
import numpy as np
from astropy.table import Table, vstack

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_CANFAR_BASE = (
    "https://ws-cadc.canfar.net/vault/files/phangs/RELEASES"
    "/Sun_etal_2022/v4p0_public_release"
)

_VERSION = "v4p0"

# Aperture type → suffix in filename
_APERTURE_SUFFIX = {
    "annulus": "annulus_0p5kpc",
    "gauss": "gauss_1p5kpc",
    "hexagon": "hexagon_1p5kpc",
}

# All 90 galaxies in the v4.0 release (alphabetical order)
_GALAXIES = [
    "ESO097-013", "IC1954", "IC5273", "IC5332",
    "NGC0253", "NGC0300", "NGC0628", "NGC0685",
    "NGC1087", "NGC1097", "NGC1300", "NGC1317", "NGC1365", "NGC1385",
    "NGC1433", "NGC1511", "NGC1512", "NGC1546", "NGC1559", "NGC1566",
    "NGC1637", "NGC1792", "NGC1809",
    "NGC2090", "NGC2283", "NGC2566", "NGC2775", "NGC2835", "NGC2903",
    "NGC2997",
    "NGC3059", "NGC3137", "NGC3239", "NGC3351", "NGC3489", "NGC3507",
    "NGC3511", "NGC3521", "NGC3596", "NGC3599", "NGC3621", "NGC3626",
    "NGC3627",
    "NGC4254", "NGC4293", "NGC4298", "NGC4303", "NGC4321", "NGC4457",
    "NGC4459", "NGC4476", "NGC4477", "NGC4496A", "NGC4535", "NGC4536",
    "NGC4540", "NGC4548", "NGC4569", "NGC4571", "NGC4596", "NGC4689",
    "NGC4731", "NGC4781", "NGC4826", "NGC4941", "NGC4951",
    "NGC5042", "NGC5068", "NGC5128", "NGC5134", "NGC5236", "NGC5248",
    "NGC5530", "NGC5643",
    "NGC6300", "NGC6744",
    "NGC7456", "NGC7496", "NGC7743", "NGC7793",
]


# ---------------------------------------------------------------------------
# Catalogue helpers
# ---------------------------------------------------------------------------


def list_galaxies():
    """Return the list of galaxy names in the v4.0 release."""
    return list(_GALAXIES)


def filename(galaxy, aperture="annulus"):
    """Return the ECSV filename for a given galaxy and aperture type.

    Parameters
    ----------
    galaxy : str
        Galaxy name, e.g. ``"NGC0628"``.
    aperture : str
        One of ``"annulus"``, ``"gauss"``, ``"hexagon"``.

    Returns
    -------
    str
    """
    if aperture not in _APERTURE_SUFFIX:
        raise ValueError(
            f"Unknown aperture '{aperture}'. "
            f"Choose from {list(_APERTURE_SUFFIX)}"
        )
    return f"{galaxy}_{_APERTURE_SUFFIX[aperture]}.ecsv"


def list_files(aperture=None):
    """Return all ECSV filenames in the release.

    Parameters
    ----------
    aperture : str or None
        If given, filter to that aperture type
        (``"annulus"``, ``"gauss"``, or ``"hexagon"``).

    Returns
    -------
    list of str
    """
    apertures = [aperture] if aperture else list(_APERTURE_SUFFIX)
    files = []
    for ap in apertures:
        for galaxy in _GALAXIES:
            files.append(filename(galaxy, aperture=ap))
    return files


def build_url(fname):
    """Return the download URL for a given filename."""
    return f"{_CANFAR_BASE}/{fname}"


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------


def download(dest_dir, aperture="annulus", galaxies=None, force=False):
    """Download PHANGS megatable ECSV files to *dest_dir*.

    Parameters
    ----------
    dest_dir : str or Path
        Directory to save files into (created if absent).
    aperture : str or ``"all"``
        Aperture type to download, or ``"all"`` for all three types.
    galaxies : list of str or None
        Subset of galaxies.  ``None`` downloads all 90.
    force : bool
        Re-download even if the file already exists locally.

    Returns
    -------
    list of Path
        Paths of downloaded files.
    """
    import urllib.request

    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    apertures = list(_APERTURE_SUFFIX) if aperture == "all" else [aperture]
    gals = galaxies if galaxies is not None else _GALAXIES

    downloaded = []
    for ap in apertures:
        for gal in gals:
            fname = filename(gal, aperture=ap)
            dest = dest_dir / fname
            if dest.exists() and not force:
                downloaded.append(dest)
                continue
            url = build_url(fname)
            try:
                urllib.request.urlretrieve(url, dest)
                downloaded.append(dest)
            except Exception as exc:
                print(f"Warning: could not download {fname}: {exc}")
    return downloaded


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------


def load(path):
    """Load a single PHANGS megatable ECSV file.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    `~astropy.table.Table`
        Table with units attached and ``GALAXY`` stored in ``.meta``.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    return Table.read(path, format="ascii.ecsv")


def load_all(data_dir, aperture="annulus"):
    """Load and vertically stack all ECSV files in *data_dir*.

    A ``GALAXY`` column is added (taken from table metadata) before stacking.

    Parameters
    ----------
    data_dir : str or Path
    aperture : str
        Filter to this aperture type.

    Returns
    -------
    `~astropy.table.Table`
    """
    data_dir = Path(data_dir)
    pattern = f"*_{_APERTURE_SUFFIX[aperture]}.ecsv"
    files = sorted(data_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(
            f"No files matching '{pattern}' found in {data_dir}"
        )
    tables = []
    for f in files:
        t = load(f)
        galaxy = t.meta.get("GALAXY", f.stem.split("_")[0])
        t["GALAXY"] = galaxy
        tables.append(t)
    return vstack(tables, metadata_conflicts="silent")


# ---------------------------------------------------------------------------
# Derived PRFM quantities
# ---------------------------------------------------------------------------


def compute_prfm_inputs(table):
    """Add PRFM-relevant derived columns to a megatable.

    Derived columns added
    ---------------------
    ``Sigma_gas``
        Total gas surface density = ``Sigma_mol + Sigma_atom``
        [M_sun / pc^2].
    ``Omega_d``
        Angular velocity = ``V_circ_CO21_URC / r_gal`` [km / s / kpc].
    ``H_star``
        Stellar scale height = ``Sigma_star / (2 * rho_star_mp)`` [pc].

    Parameters
    ----------
    table : `~astropy.table.Table`

    Returns
    -------
    `~astropy.table.Table`
        A copy of the input table with the new columns appended.
    """
    t = table.copy()

    Sigma_mol = t["Sigma_mol"].to(au.M_sun / au.pc**2)
    Sigma_atom = t["Sigma_atom"].to(au.M_sun / au.pc**2)
    t["Sigma_gas"] = (Sigma_mol + Sigma_atom).to(au.M_sun / au.pc**2)
    t["Sigma_gas"].description = "Total gas surface density (mol + atom)"

    V_circ = t["V_circ_CO21_URC"].to(au.km / au.s)
    r_gal = t["r_gal"].to(au.kpc)
    t["Omega_d"] = (V_circ / r_gal).to(au.km / au.s / au.kpc)
    t["Omega_d"].description = "Angular velocity V_circ / r_gal"

    Sigma_star = t["Sigma_star"].to(au.M_sun / au.pc**2)
    rho_star = t["rho_star_mp"].to(au.M_sun / au.pc**3)
    t["H_star"] = (Sigma_star / (2.0 * rho_star)).to(au.pc)
    t["H_star"].description = "Stellar scale height Sigma_star / (2 * rho_star_mp)"

    return t


# ---------------------------------------------------------------------------
# Quality masking
# ---------------------------------------------------------------------------


def valid_rows(table, cols=None):
    """Boolean mask selecting rows that are finite and positive.

    Parameters
    ----------
    table : `~astropy.table.Table`
    cols : list of str or None
        Columns to check.  ``None`` checks all numeric columns.

    Returns
    -------
    numpy.ndarray of bool
    """
    if cols is None:
        cols = [
            c for c in table.colnames
            if table[c].dtype.kind in ("f", "i")
        ]
    mask = np.ones(len(table), dtype=bool)
    for col in cols:
        vals = np.asarray(table[col], dtype=float)
        mask &= np.isfinite(vals) & (vals > 0)
    return mask


# ---------------------------------------------------------------------------
# PRFM model application
# ---------------------------------------------------------------------------

# Unit conversion factors (astro → CGS) matching prfm.prfm conventions
_surf_cgs = (ac.M_sun / ac.pc**2).cgs.value      # M_sun/pc² → g/cm²
_pc_cgs = ac.pc.cgs.value                          # pc → cm
_kms_kpc_cgs = 1.0e5 / (3.085677581e21)            # km/s/kpc → 1/s
_kbol_cgs = ac.k_B.cgs.value                       # erg/K
_sfr_cgs = (ac.M_sun / ac.kpc**2 / au.yr).cgs.value  # M_sun/kpc²/yr → g/cm²/s


def run_prfm(
    table: Table,
    sigma_eff_model: str = "tigress-ncr-avg",
    yield_model: str = "tigress-ncr-decomp-all",
    zprime_col: Optional[str] = "Zprime",
) -> Table:
    """Apply the PRFM model to every row of a PHANGS megatable.

    Runs :func:`compute_prfm_inputs` if ``Sigma_gas`` is not already present,
    then calls the self-consistent PRFM solver (vectorized over all rows).
    Invalid rows (NaN or non-positive PRFM inputs) receive ``NaN`` outputs.

    Columns added
    -------------
    ``P_weight``
        Dynamical equilibrium pressure (total gas weight)
        in units of :math:`k_\\mathrm{B}` K cm\\ :sup:`-3`.
    ``H_gas``
        Self-consistent gas scale height [pc].
    ``sigma_eff_sol``
        Self-consistent effective velocity dispersion [km s\\ :sup:`-1`].
    ``Sigma_SFR_pred``
        Predicted SFR surface density [M_sun kpc\\ :sup:`-2` yr\\ :sup:`-1`].

    Parameters
    ----------
    table : `~astropy.table.Table`
        PHANGS megatable, optionally pre-processed by :func:`compute_prfm_inputs`.
    sigma_eff_model : str
        Velocity dispersion model name (key in ``prfm._sigma_eff_models``).
    yield_model : str
        Feedback yield model name (key in ``prfm._yield_models``).
    zprime_col : str or None
        Column name for metallicity (relative to solar).  Pass ``None`` to
        ignore metallicity dependence.

    Returns
    -------
    `~astropy.table.Table`
        Copy of the input table with four new columns appended.
    """
    from prfm.prfm import get_self_consistent_solution, get_sfr

    t = table.copy()

    # Ensure derived PRFM inputs exist
    if "Sigma_gas" not in t.colnames:
        t = compute_prfm_inputs(t)

    # Identify valid rows
    prfm_cols = ["Sigma_gas", "Sigma_star", "Omega_d", "H_star"]
    mask = valid_rows(t, cols=prfm_cols)
    n_invalid = (~mask).sum()
    if n_invalid > 0:
        warnings.warn(
            f"run_prfm: {n_invalid} row(s) with invalid PRFM inputs will be NaN.",
            stacklevel=2,
        )

    # Initialise output arrays with NaN
    n = len(t)
    P_weight = np.full(n, np.nan)
    H_gas = np.full(n, np.nan)
    sigma_sol = np.full(n, np.nan)
    Sigma_SFR_pred = np.full(n, np.nan)

    if mask.any():
        # Extract valid rows and convert to CGS
        sg = np.asarray(t["Sigma_gas"][mask], dtype=float) * _surf_cgs
        ss = np.asarray(t["Sigma_star"][mask], dtype=float) * _surf_cgs
        od = np.asarray(t["Omega_d"][mask], dtype=float) * _kms_kpc_cgs
        hs = np.asarray(t["H_star"][mask], dtype=float) * _pc_cgs

        # Metallicity (dimensionless relative to solar); default to solar (Z=1) if absent
        if zprime_col is not None and zprime_col in t.colnames:
            Z: np.ndarray = np.asarray(t[zprime_col][mask], dtype=float)
        else:
            Z = np.ones(int(mask.sum()))

        # Solve
        wtot, H_cgs, se_cgs = get_self_consistent_solution(
            sg, ss, od, hs, sigma_eff=sigma_eff_model
        )

        # Back-convert to observer-friendly units
        P_weight[mask] = wtot / _kbol_cgs          # k_B K cm^-3
        H_gas[mask] = H_cgs / _pc_cgs              # pc
        sigma_sol[mask] = se_cgs / 1.0e5           # km/s

        # Predicted SFR
        sfr_cgs = get_sfr(wtot, Z=Z, Ytot=yield_model)
        Sigma_SFR_pred[mask] = sfr_cgs / _sfr_cgs  # M_sun/kpc²/yr

    # Attach columns with astropy units
    t["P_weight"] = P_weight * au.K / au.cm**3
    t["P_weight"].description = "Dynamical equilibrium pressure (P_DE / k_B)"

    t["H_gas"] = H_gas * au.pc
    t["H_gas"].description = "Self-consistent gas scale height"

    t["sigma_eff_sol"] = sigma_sol * au.km / au.s
    t["sigma_eff_sol"].description = "Self-consistent effective velocity dispersion"

    t["Sigma_SFR_pred"] = Sigma_SFR_pred * au.M_sun / au.kpc**2 / au.yr
    t["Sigma_SFR_pred"].description = "PRFM-predicted SFR surface density"

    return t
