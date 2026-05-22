"""
Shared plotting and data-preparation utilities for PHANGS megatable exploration.

Imported by scripts that produce exploratory figures; not part of the public
PRFM computation API.
"""

from collections.abc import Sequence
from typing import Any, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

# ---------------------------------------------------------------------------
# Aperture configuration
# ---------------------------------------------------------------------------

# Aperture type → column-name suffix used in the megatable
APERTURE_SUFFIX: dict[str, str] = {
    "annulus": "",
    "gauss": "_gauss",
    "hexagon": "_hexagon",
}

# Columns that never carry an aperture suffix (shared across all apertures)
NO_SUFFIX: set[str] = {
    "r_gal",
    "Zprime",
    "alpha_CO21_S20",
    "rho_star_mp",
    "V_circ_CO21_URC",
    "Omega_d",
    "H_star",
    "Sigma_gas",  # derived by compute_prfm_inputs
}

# ---------------------------------------------------------------------------
# Column label registry  (base name → LaTeX label with units)
# ---------------------------------------------------------------------------

COLUMN_LABELS: dict[str, str] = {
    "r_gal": r"$r_\mathrm{gal}$ [kpc]",
    "Sigma_mol": r"$\Sigma_\mathrm{mol}$ [$M_\odot\,\mathrm{pc}^{-2}$]",
    "Sigma_atom": r"$\Sigma_\mathrm{atom}$ [$M_\odot\,\mathrm{pc}^{-2}$]",
    "Sigma_gas": r"$\Sigma_\mathrm{gas}$ [$M_\odot\,\mathrm{pc}^{-2}$]",
    "Sigma_star": r"$\Sigma_\star$ [$M_\odot\,\mathrm{pc}^{-2}$]",
    "rho_star_mp": r"$\rho_\star$ [$M_\odot\,\mathrm{pc}^{-3}$]",
    "H_star": r"$H_\star$ [pc]",
    "V_circ_CO21_URC": r"$V_\mathrm{circ}$ [km s$^{-1}$]",
    "Omega_d": r"$\Omega_d$ [km s$^{-1}$ kpc$^{-1}$]",
    "Zprime": r"$Z'$ [$Z_\odot$]",
    "Sigma_SFR_HaW4recal": r"$\Sigma_\mathrm{SFR}^\mathrm{H\alpha+W4}$"
    r" [$M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}$]",
    "Sigma_SFR_FUVW4recal": r"$\Sigma_\mathrm{SFR}^\mathrm{FUV+W4}$"
    r" [$M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}$]",
    "Sigma_SFR_Hacorr": r"$\Sigma_\mathrm{SFR}^\mathrm{H\alpha,corr}$"
    r" [$M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}$]",
    "alpha_CO21_S20": r"$\alpha_\mathrm{CO}$ [$M_\odot\,\mathrm{pc}^{-2}\,(\mathrm{K\,km\,s}^{-1})^{-1}$]",
    "P_weight": r"$P_\mathrm{DE}/k_B$ [K cm$^{-3}$]",
    "H_gas": r"$H_\mathrm{gas}$ [pc]",
    "sigma_eff_sol": r"$\sigma_\mathrm{eff}$ [km s$^{-1}$]",
    "Sigma_SFR_pred": r"$\Sigma_\mathrm{SFR}^\mathrm{pred}$"
    r" [$M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}$]",
}


def col_label(col: str) -> str:
    """Return a LaTeX label for *col*, stripping aperture suffix if needed."""
    if col in COLUMN_LABELS:
        return COLUMN_LABELS[col]
    # Strip known aperture suffixes and try again
    for sfx in ("_gauss", "_hexagon"):
        if col.endswith(sfx):
            base = col[: -len(sfx)]
            if base in COLUMN_LABELS:
                return COLUMN_LABELS[base]
    return col


# ---------------------------------------------------------------------------
# Radial bin configuration
# ---------------------------------------------------------------------------

RGAL_BINS: list[tuple[float, float]] = [(0, 2), (2, 5), (5, 10), (10, 999)]
RGAL_LABELS: list[str] = [
    r"$r<2$ kpc",
    r"$2-5$ kpc",
    r"$5-10$ kpc",
    r"$r>10$ kpc",
]
RGAL_COLORS: list[str] = ["#c0392b", "#e67e22", "#27ae60", "#2980b9"]

# ---------------------------------------------------------------------------
# Column-resolution helper
# ---------------------------------------------------------------------------


def resolve_columns(
    base_fields: list[tuple[str, str, bool]],
    table: Table,
    aperture: str = "annulus",
) -> list[tuple[str, str, bool]]:
    """Resolve base column names to actual names present in *table*.

    For each ``(base_name, label, log_scale)`` entry, the suffix for the
    given *aperture* is appended unless the base name is in :data:`NO_SUFFIX`.
    If the suffixed name is absent, falls back to the unsuffixed name.
    Entries whose resolved name is absent from *table* are dropped.

    Parameters
    ----------
    base_fields:
        List of ``(base_name, axis_label, log_scale)`` tuples.
    table:
        Loaded PHANGS megatable (e.g. from :func:`prfm.phangs.load_all`).
    aperture:
        One of ``"annulus"``, ``"gauss"``, ``"hexagon"``.

    Returns
    -------
    list of (col_name, label, log_scale) with concrete column names.
    """
    sfx = APERTURE_SUFFIX[aperture]
    resolved: list[tuple[str, str, bool]] = []
    for base, label, log in base_fields:
        col = base if base in NO_SUFFIX else base + sfx
        if col in table.colnames:
            resolved.append((col, label, log))
        elif base in table.colnames:
            resolved.append((base, label, log))
    return resolved


# ---------------------------------------------------------------------------
# Array helpers
# ---------------------------------------------------------------------------


def to_log_or_raw(
    table: Table,
    col: str,
    log: bool,
) -> np.ndarray:
    """Return column values, optionally log10-transformed.

    Non-positive values are mapped to ``NaN`` before taking the logarithm.

    Parameters
    ----------
    table:
        Astropy table containing *col*.
    col:
        Column name.
    log:
        If ``True``, return ``log10(values)``; otherwise return raw values.

    Returns
    -------
    1-D float array of length ``len(table)``.
    """
    v = np.asarray(table[col], dtype=float)
    if log:
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.log10(np.where(v > 0, v, np.nan))
    return v


def build_matrix(
    table: Table,
    cols: list[tuple[str, str, bool]],
) -> tuple[np.ndarray, list[str]]:
    """Build a 2-D data array restricted to rows that are finite in all columns.

    Parameters
    ----------
    table:
        Astropy table.
    cols:
        List of ``(col_name, label, log_scale)`` tuples.

    Returns
    -------
    data : ndarray of shape ``(N_valid, len(cols))``
    labels : list of axis label strings
    """
    arrays = [to_log_or_raw(table, col, log) for col, _, log in cols]
    stacked = np.column_stack(arrays)
    valid = np.all(np.isfinite(stacked), axis=1)
    labels = [lbl for _, lbl, _ in cols]
    return stacked[valid], labels


def clean(vals: np.ndarray, log: bool) -> np.ndarray:
    """Return finite (and positive when *log* is ``True``) values.

    Parameters
    ----------
    vals:
        Input array (will be cast to ``float``).
    log:
        If ``True``, also require ``vals > 0``.

    Returns
    -------
    1-D float array of valid values.
    """
    v = np.asarray(vals, dtype=float)
    mask = np.isfinite(v)
    if log:
        mask &= v > 0
    return v[mask]


def bin_edges(
    vals: np.ndarray,
    log: bool,
    n: int = 40,
) -> np.ndarray:
    """Compute histogram bin edges for *vals*.

    Parameters
    ----------
    vals:
        Array of finite values (no NaNs expected).
    log:
        If ``True``, produce logarithmically spaced edges.
    n:
        Number of bins.

    Returns
    -------
    1-D array of ``n + 1`` bin edges.
    """
    if log:
        return np.logspace(np.log10(vals.min()), np.log10(vals.max()), n + 1)
    return np.linspace(vals.min(), vals.max(), n + 1)


# ---------------------------------------------------------------------------
# Scatter plot
# ---------------------------------------------------------------------------


def _asymmetric_errorbars(
    vals: np.ndarray,
    errs: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert symmetric errors to asymmetric (lo, hi) pairs for log-safe plotting.

    Clips the lower bar so it never reaches or crosses zero.

    Parameters
    ----------
    vals:
        Central values (positive, finite).
    errs:
        Symmetric 1-sigma uncertainties (same units as *vals*).

    Returns
    -------
    (err_lo, err_hi) each of the same shape as *vals*.
    """
    err_hi = errs.copy()
    err_lo = np.minimum(errs, vals * 0.9999)  # clip so lower bar stays positive
    return err_lo, err_hi


def _axis_limits(
    vals: np.ndarray,
    log: bool,
    plo: float = 1.0,
    phi: float = 99.0,
    margin: float = 0.5,
) -> tuple[float, float]:
    """Compute tight axis limits from the *plo*–*phi* percentile range of *vals*.

    Parameters
    ----------
    vals:
        Finite, positive data values.
    log:
        If ``True``, apply margin in log10 space.
    plo, phi:
        Lower and upper percentile bounds (default 1st–99th).
    margin:
        Fractional padding beyond the percentile range.
        In log space this is in dex; in linear space it is a fraction of the
        total range.

    Returns
    -------
    (lo, hi) axis limit pair.
    """
    v = vals[np.isfinite(vals) & (vals > 0)]
    if len(v) == 0:
        return (None, None)  # type: ignore[return-value]
    lo, hi = np.nanpercentile(v, [plo, phi])
    if log:
        log_lo = np.log10(lo) - margin
        log_hi = np.log10(hi) + margin
        return 10**log_lo, 10**log_hi
    span = hi - lo
    return lo - margin * span, hi + margin * span


def scatter_plot(
    table: Table,
    xcol: str,
    ycol: str,
    ax: Optional[plt.Axes] = None,
    log_x: bool = True,
    log_y: bool = True,
    errorbars: bool = True,
    bg_color: str = "0.6",
    sel_color: str = "tab:blue",
    slice_frac: float = 0.1,
    bg_alpha: float = 0.15,
    sel_alpha: float = 0.7,
    s: float = 4,
    **slice_kwargs: Union[float, Sequence[float]],
) -> plt.Axes:
    """Scatter plot of *ycol* vs *xcol* with optional error bars and slice highlighting.

    Parameters
    ----------
    table:
        PHANGS megatable (astropy Table).
    xcol:
        Column name for the x-axis.
    ycol:
        Column name for the y-axis.
    ax:
        Axes to draw on.  A new figure/axes is created if ``None``.
    log_x:
        Use log scale on x-axis.
    log_y:
        Use log scale on y-axis.
    errorbars:
        Draw error bars when ``e_<col>`` columns are present.
    bg_color:
        Colour for unselected (background) points.
    sel_color:
        Colour for points in the highlighted slice.
    slice_frac:
        Half-width of the slice as a fraction of the slice value.
        E.g. ``slice_frac=0.1`` with ``Sigma_gas=10`` selects rows
        where ``9 <= Sigma_gas <= 11``.
    bg_alpha:
        Alpha for background points.
    sel_alpha:
        Alpha for highlighted points.
    s:
        Marker size.
    **slice_kwargs:
        One optional keyword of the form ``column_name=value`` (scalar) or
        ``column_name=[v1, v2, ...]`` (multiple values) to highlight one or
        more narrow slices.  Example::

            scatter_plot(t, "Sigma_gas", "Sigma_SFR_HaW4recal", Sigma_gas=10)
            scatter_plot(t, "Sigma_gas", "Sigma_SFR_HaW4recal", Sigma_gas=[10, 30])

    Returns
    -------
    `~matplotlib.axes.Axes`
    """
    if ax is None:
        _, ax = plt.subplots()

    x = np.asarray(table[xcol], dtype=float)
    y = np.asarray(table[ycol], dtype=float)
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)

    slices = _parse_slices(table, valid, slice_kwargs, slice_frac, sel_color)
    any_sel = np.zeros(len(table), dtype=bool)
    for mask, _, _ in slices:
        any_sel |= mask
    bg_mask = valid & ~any_sel

    # Error arrays (None if column absent or errorbars disabled)
    def _errs(col: str) -> Optional[np.ndarray]:
        ecol = f"e_{col}"
        if errorbars and ecol in table.colnames:
            return np.asarray(table[ecol], dtype=float)
        return None

    ex = _errs(xcol)
    ey = _errs(ycol)

    def _plot_group(mask: np.ndarray, color: Any, alpha: float, zorder: int) -> None:
        if not mask.any():
            return
        xm, ym = x[mask], y[mask]
        if ex is not None or ey is not None:
            xerr = None if ex is None else get_symmetric_log_errorbars(xm, ex[mask])
            yerr = None if ey is None else get_symmetric_log_errorbars(ym, ey[mask])
            ax.errorbar(
                xm,
                ym,
                xerr=xerr,
                yerr=yerr,
                fmt="o",
                ms=np.sqrt(s),
                color=color,
                alpha=alpha,
                elinewidth=0.5,
                capsize=0,
                zorder=zorder,
                rasterized=True,
            )
        else:
            ax.scatter(
                xm,
                ym,
                s=s,
                color=color,
                alpha=alpha,
                zorder=zorder,
                rasterized=True,
                linewidths=0,
            )

    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")

    _plot_group(bg_mask, bg_color, bg_alpha, zorder=1)
    for i, (mask, color, _) in enumerate(slices):
        _plot_group(mask, color, sel_alpha, zorder=2 + i)

    # Tighten limits based on the data distribution, not the error bars
    xlim = _axis_limits(x[valid], log=log_x)
    ylim = _axis_limits(y[valid], log=log_y)
    if xlim[0] is not None:
        ax.set_xlim(xlim)
    if ylim[0] is not None:
        ax.set_ylim(ylim)

    ax.set_xlabel(col_label(xcol))
    ax.set_ylabel(col_label(ycol))

    return ax


def _pdf_step(
    vals: np.ndarray,
    edges: np.ndarray,
    total: int,
    log: bool = False,
) -> np.ndarray:
    """Compute PDF as count / total / delta_bin for each bin.

    Parameters
    ----------
    vals:
        Data values to histogram (finite, already filtered).
    edges:
        Bin edges (length n_bins + 1).
    total:
        Total count used for normalisation (typically the full-sample size).
    log:
        If ``True``, compute bin widths in log10 space so the PDF integrates
        to 1 when the x-axis is logarithmic.

    Returns
    -------
    pdf : ndarray of length n_bins
    """
    counts, _ = np.histogram(vals, bins=edges)
    if log:
        delta = np.log10(edges[1:]) - np.log10(edges[:-1])
    else:
        delta = edges[1:] - edges[:-1]
    return counts / total / delta


def hist_plot(
    table: Table,
    col: str,
    ax: Optional[plt.Axes] = None,
    log: bool = True,
    n_bins: int = 40,
    bg_color: str = "0.6",
    sel_color: str = "tab:blue",
    slice_frac: float = 0.1,
    bg_alpha: float = 0.6,
    sel_alpha: float = 0.85,
    **slice_kwargs: Union[float, Sequence[float]],
) -> plt.Axes:
    """Histogram of *col* with optional slice highlighting.

    The PDF is computed as ``count / N_total / delta_bin`` so the integral
    over all bins equals 1 for the background distribution.  Both
    distributions share the same bin edges and normalisation base (*N_total*
    of the full sample), making their shapes directly comparable.

    Parameters
    ----------
    table:
        PHANGS megatable (astropy Table).
    col:
        Column to histogram.
    ax:
        Axes to draw on.  A new figure/axes is created if ``None``.
    log:
        If ``True``, use logarithmic x-axis and log-spaced bins.
    n_bins:
        Number of histogram bins.
    bg_color:
        Colour for the full-distribution step.
    sel_color:
        Colour for the highlighted-slice step.
    slice_frac:
        Half-width of the slice as a fraction of the slice value
        (default ±10%).
    bg_alpha:
        Alpha for the background step.
    sel_alpha:
        Alpha for the slice step.
    **slice_kwargs:
        One optional keyword of the form ``column_name=value`` (scalar) or
        ``column_name=[v1, v2, ...]`` (multiple values) to highlight one or
        more narrow slices.  Example::

            hist_plot(t, "Sigma_SFR_HaW4recal", Sigma_gas=10)
            hist_plot(t, "Sigma_SFR_HaW4recal", Sigma_gas=[10, 30])

    Returns
    -------
    `~matplotlib.axes.Axes`
    """
    if ax is None:
        _, ax = plt.subplots()

    v = np.asarray(table[col], dtype=float)
    valid = np.isfinite(v) & (v > 0) if log else np.isfinite(v)
    v_all = v[valid]

    slices = _parse_slices(table, valid, slice_kwargs, slice_frac, sel_color)

    edges = bin_edges(v_all, log=log, n=n_bins)
    n_total = len(v_all)

    pdf_all = _pdf_step(v_all, edges, n_total, log=log)
    ax.step(
        edges[:-1],
        pdf_all,
        where="post",
        color=bg_color,
        alpha=bg_alpha,
        label="all",
        zorder=1,
    )

    for i, (mask, color, label) in enumerate(slices):
        if mask.any():
            pdf_sel = _pdf_step(v[mask], edges, n_total, log=log)
            ax.step(
                edges[:-1],
                pdf_sel,
                where="post",
                color=color,
                alpha=sel_alpha,
                zorder=2 + i,
                label=label,
            )

    if slices:
        ax.legend(fontsize="small")

    if log:
        ax.set_xscale("log")

    ax.set_xlabel(col_label(col))
    ax.set_ylabel("PDF")

    return ax


def _parse_slices(
    table: Table,
    valid: np.ndarray,
    slice_kwargs: dict[str, Any],
    slice_frac: float,
    sel_color: str,
) -> list[tuple[np.ndarray, str, str]]:
    """Resolve slice kwargs into a list of (mask, color, label) tuples.

    Accepts scalar or sequence values for each slice kwarg, e.g.::

        Sigma_gas=10           → one slice
        Sigma_gas=[10, 30]     → two slices, colors from tab10

    Parameters
    ----------
    table:
        Megatable.
    valid:
        Boolean array of rows that are finite/positive for the plotted columns.
    slice_kwargs:
        Raw ``**kwargs`` from the caller — values may be scalar or sequence.
    slice_frac:
        Fractional half-width of each slice window.
    sel_color:
        Colour used when there is exactly one slice value.

    Returns
    -------
    List of ``(mask, color, label)`` — one entry per slice value.
    """
    results: list[tuple[np.ndarray, str, str]] = []
    for scol, raw_val in slice_kwargs.items():
        if scol not in table.colnames:
            continue
        vals: list[float] = (
            list(raw_val)
            if isinstance(raw_val, Sequence) and not isinstance(raw_val, str)
            else [float(raw_val)]
        )
        n = len(vals)
        colors = [sel_color] if n == 1 else [plt.cm.tab10(i / 10) for i in range(n)]
        sym = COLUMN_LABELS.get(scol, scol).split("[")[0].strip()
        sv = np.asarray(table[scol], dtype=float)
        for val, color in zip(vals, colors):
            lo = val * (1.0 - slice_frac)
            hi = val * (1.0 + slice_frac)
            mask = valid & (sv >= lo) & (sv <= hi)
            label = rf"{sym} $\in$ [{lo:.3g}, {hi:.3g}]"
            results.append((mask, color, label))
        break  # one slice column supported
    return results


def get_symmetric_log_errorbars(
    linear_mean: np.ndarray,
    linear_std: np.ndarray,
) -> list[np.ndarray]:
    """Convert linear-space std to visually symmetric error bars on a log axis.

    Converts a standard deviation measured in linear space into a
    multiplicative factor (exp(σ/μ)) so that the resulting error bars
    appear symmetric when the axis is log-scaled.

    Parameters
    ----------
    linear_mean : array_like
        Mean values in linear space.
    linear_std : array_like
        Standard deviation values in linear space.

    Returns
    -------
    list of ndarray
        Two-element list ``[lower_errors, upper_errors]`` containing the
        absolute linear distances from the centre point, as expected by
        Matplotlib's ``yerr`` parameter.
    """
    linear_mean = np.asarray(linear_mean, dtype=float)
    linear_std = np.asarray(linear_std, dtype=float)

    # 1. Calculate the Coefficient of Variation (fractional error)
    cv = linear_std / linear_mean

    # 2. Calculate the multiplicative factor
    # Using np.exp ensures the math works perfectly regardless of the log base
    # Matplotlib uses for its axis scaling.
    multiplier = np.exp(cv)

    # 3. Calculate the new bounds based on the multiplier
    upper_bound = linear_mean * multiplier
    lower_bound = linear_mean / multiplier

    # 4. Calculate the distances from the center point for Matplotlib
    yerr_lower = linear_mean - lower_bound
    yerr_upper = upper_bound - linear_mean

    return [yerr_lower, yerr_upper]


def plot_weights(
    tbl: "Table",
    ax: "plt.Axes | None" = None,
    variation: dict | None = None,
) -> "plt.Figure":
    """Plot weight fractions (f_gas, f_star, f_DM) vs dynamical equilibrium pressure.

    Scatter plots each weight fraction against P_DE (from the ``P_weight``
    column) and overlays binned geometric-mean lines for each component.

    Parameters
    ----------
    tbl : Table
        PHANGS table with ``P_weight``, ``Sigma_gas``, ``Sigma_star``,
        ``Omega_d``, ``H_star``, and ``sigma_eff_sol`` columns.
    ax : Axes or None, optional
        Existing Matplotlib axes to draw on. A new figure is created if
        ``None``.
    variation : dict or None, optional
        Passed directly to :func:`prfm.phangs.get_weights`; see that
        function for supported keys.

    Returns
    -------
    Figure
        The Matplotlib figure containing the plot.
    """
    from prfm.phangs import get_weights

    f_gas, f_star, f_dm = get_weights(tbl, variation=variation)
    P_DE = np.asarray(tbl["P_weight"], dtype=float)  # k_B K cm^-3

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 3.5))
    else:
        fig = ax.figure
    kw = dict(s=2, alpha=0.5, rasterized=True, linewidths=0)
    ax.scatter(P_DE, f_gas, color="tab:blue", **kw)
    ax.scatter(P_DE, f_star, color="tab:orange", **kw)
    ax.scatter(P_DE, f_dm, color="tab:green", **kw)

    # Binned means
    P_edges = np.logspace(
        np.log10(np.nanpercentile(P_DE, 1)), np.log10(np.nanpercentile(P_DE, 99)), 20
    )
    P_centers = np.sqrt(P_edges[:-1] * P_edges[1:])  # geometric mid-points
    for f, color, label in [
        (f_gas, "tab:blue", r"$f_\mathrm{gas}$"),
        (f_star, "tab:orange", r"$f_\star$"),
        (f_dm, "tab:green", r"$f_\mathrm{DM}$"),
    ]:
        means = [
            np.nanmean(f[(P_DE >= lo) & (P_DE < hi)])
            for lo, hi in zip(P_edges[:-1], P_edges[1:])
        ]
        ax.plot(P_centers, means, color=color, lw=2, label=label)

    n_valid = np.isfinite(P_DE).sum()
    ax.set_xscale("log")
    ax.set_xlabel(r"$P_\mathrm{DE}/k_B\ [\mathrm{K\,cm^{-3}}]$")
    ax.set_ylabel(r"$\mathcal{W}_i\,/\,\mathcal{W}_\mathrm{tot}$")
    ax.set_ylim(0, 1)
    ax.legend(markerscale=5)
    ax.annotate(
        f"N={n_valid}",
        xy=(0.97, 0.97),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize="x-small",
        color="0.4",
    )
    fig.tight_layout()
    return fig
