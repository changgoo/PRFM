#!/usr/bin/env python3
"""Compare TIGRESS-NCR SFR and molecular fraction with the PHANGS parent band.

The completed core suite samples a finite-width Sigma_gas band.  At fixed
Sigma_gas its four varied environmental inputs are Sigma_star, H_star, Omega,
and qshear.  This script plots the two validation outcomes against those four
inputs and the derived quantities

    kappa = sqrt(2 (2 - qshear)) Omega
    rho_star = Sigma_star / (2 H_star).

The simulated molecular fraction is measured from the native history volume
integrals as 2 * scalar3 / mass.  ``scalar3`` is the H2 number abundance in the
five-scalar TIGRESS-NCR build used for this suite.  Every history sample from
400--600 Myr is plotted by default, with model time means overlaid.  A scatter
budget separates between-environment and within-model temporal variance in
log space and compares their quadrature sum with the PHANGS scatter.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
plt.style.use(SCRIPTS / "prfm.mplstyle")

from prfm.phangs import load_configured_phangs  # noqa: E402
from prfm.phangs_sampling import (  # noqa: E402
    PHANGSSamplingDesigner,
    SamplingConfig,
)

try:  # noqa: E402
    from pathena.hst_reader import read_hst
except ImportError as exc:  # pragma: no cover - environment diagnostic
    raise SystemExit(
        "pathena is required to read TIGRESS histories; use the pyathena "
        "environment or install tigress_ncr_tools"
    ) from exc


DEFAULT_SUITE = Path("/tigress/changgoo/anvil/TIGRESS-NCR-suite")
DEFAULT_SUMMARY = DEFAULT_SUITE / "prfm_diagnostics/prfm_model_summary.csv"
DEFAULT_CONFIG = ROOT / "config/phangs_prfm.yml"
DESIGN_FIELDS = ["Sigma_gas", "Sigma_star", "H_star", "Omega", "qshear"]
PREDICTORS = ["Sigma_star", "H_star", "Omega", "qshear", "kappa", "rho_star"]

LABELS = {
    "Sigma_star": r"$\Sigma_\star$ [$M_\odot\,{\rm pc}^{-2}$]",
    "H_star": r"$H_\star$ [pc]",
    "Omega": r"$\Omega$ [km s$^{-1}$ kpc$^{-1}$]",
    "qshear": r"$q$",
    "kappa": r"$\kappa$ [km s$^{-1}$ kpc$^{-1}$]",
    "rho_star": r"$\rho_\star=\Sigma_\star/(2H_\star)$ "
    r"[$M_\odot\,{\rm pc}^{-3}$]",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="Core-suite design CSV")
    parser.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--sigma-gas", type=float, default=None)
    parser.add_argument("--delta", type=float, default=0.3)
    parser.add_argument("--time-min", type=float, default=400.0)
    parser.add_argument("--time-max", type=float, default=600.0)
    parser.add_argument("--sfr-col", default="sfr40_mean")
    parser.add_argument("--obs-sfr-col", default="Sigma_SFR_HaW4recal")
    parser.add_argument("--h2-scalar", default="scalar3")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--figure-dir", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=180)
    return parser.parse_args()


def parse_sigma_gas(stem: str) -> float | None:
    match = re.search(r"Sgas([0-9.]+)", stem)
    return float(match.group(1)) if match else None


def row_indices(models: pd.Series) -> pd.Series:
    values = models.astype(str).str.extract(r"row(\d+)", expand=False)
    if values.isna().any():
        bad = models[values.isna()].tolist()
        raise ValueError(f"model names lack row suffixes: {bad}")
    return values.astype(int)


def add_environment_derived(frame: pd.DataFrame) -> pd.DataFrame:
    """Add kappa and the simulation-input stellar midplane density."""
    result = frame.copy()
    q = result["qshear"].to_numpy(dtype=float)
    omega = result["Omega"].to_numpy(dtype=float)
    epicycle = 2.0 * (2.0 - q)
    result["kappa"] = np.where(epicycle >= 0, np.sqrt(epicycle) * omega, np.nan)
    result["rho_star"] = result["Sigma_star"] / (2.0 * result["H_star"])
    return result


def primary_history(model: Path) -> Path:
    paths = [
        path
        for path in sorted((model / "hst").glob("*.hst"))
        if ".phase" not in path.name and ".whole" not in path.name
    ]
    if len(paths) != 1:
        raise FileNotFoundError(
            f"expected one primary history under {model / 'hst'}, found {len(paths)}"
        )
    return paths[0]


def _window_series(
    time: np.ndarray,
    values: np.ndarray,
    bounds: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """Return sorted samples with interpolated values at both time bounds."""
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(time) & np.isfinite(values)
    time, values = time[finite], values[finite]
    order = np.argsort(time, kind="stable")
    time, values = time[order], values[order]
    keep = np.r_[time[1:] > time[:-1], True]
    time, values = time[keep], values[keep]
    lower, upper = bounds
    if time.size < 2 or time[0] > lower or time[-1] < upper:
        covered = (time[0], time[-1]) if time.size else (np.nan, np.nan)
        raise ValueError(f"history covers {covered}, not [{lower}, {upper}]")
    inside = (time > lower) & (time < upper)
    sample_time = np.r_[lower, time[inside], upper]
    sample_values = np.r_[
        np.interp(lower, time, values), values[inside], np.interp(upper, time, values)
    ]
    return sample_time, sample_values


def history_diagnostics(
    history_path: Path,
    bounds: tuple[float, float],
    h2_scalar: str,
    sfr_history_col: str,
) -> tuple[dict[str, float | int], pd.DataFrame]:
    """Return model summary and every SFR/fmol sample in a time window."""
    history = read_hst(history_path)
    required = ["time", "mass", h2_scalar, sfr_history_col]
    missing = [name for name in required if name not in history]
    if missing:
        raise KeyError(f"{history_path} lacks history columns {missing}")
    mass = np.asarray(history["mass"], dtype=float)
    h2 = np.asarray(history[h2_scalar], dtype=float)
    fraction = np.divide(
        2.0 * h2,
        mass,
        out=np.full_like(mass, np.nan),
        where=np.isfinite(mass) & (mass > 0),
    )
    time, fraction = _window_series(history["time"], fraction, bounds)
    sfr_time, sfr = _window_series(history["time"], history[sfr_history_col], bounds)
    if not np.array_equal(time, sfr_time):
        raise ValueError(f"SFR and molecular history times differ in {history_path}")
    if np.any((fraction < -1.0e-8) | (fraction > 1.0 + 1.0e-8)):
        raise ValueError(f"unphysical molecular fractions in {history_path}")
    fraction = np.clip(fraction, 0.0, 1.0)
    fmol_mean = float(np.trapz(fraction, time) / (bounds[1] - bounds[0]))
    fmol_p16, fmol_p50, fmol_p84 = np.percentile(fraction, [16.0, 50.0, 84.0])
    sfr_mean = float(np.trapz(sfr, time) / (bounds[1] - bounds[0]))
    sfr_p16, sfr_p50, sfr_p84 = np.percentile(sfr, [16.0, 50.0, 84.0])
    summary = {
        "fmol_mean": fmol_mean,
        "fmol_p16": float(fmol_p16),
        "fmol_p50": float(fmol_p50),
        "fmol_p84": float(fmol_p84),
        "fmol_samples": int(len(fraction)),
        "sfr_history_mean": sfr_mean,
        "sfr_history_p16": float(sfr_p16),
        "sfr_history_p50": float(sfr_p50),
        "sfr_history_p84": float(sfr_p84),
        "sfr_history_samples": int(len(sfr)),
    }
    samples = pd.DataFrame({"time": time, "Sigma_SFR": sfr, "fmol": fraction})
    return summary, samples


def build_suite_table(
    design: pd.DataFrame,
    summary: pd.DataFrame,
    suite: Path,
    bounds: tuple[float, float],
    h2_scalar: str,
    sfr_col: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing = [field for field in DESIGN_FIELDS if field not in design]
    if missing:
        raise ValueError(f"design CSV lacks columns {missing}")
    if "model" not in summary or sfr_col not in summary:
        raise ValueError(f"summary must contain model and {sfr_col}")

    summary = summary.copy()
    summary["row"] = row_indices(summary["model"])
    if summary["row"].duplicated().any():
        raise ValueError("summary contains duplicate model row indices")
    joined = design.reset_index(drop=True).copy()
    joined["row"] = np.arange(len(joined))
    keep = [
        "row",
        "model",
        sfr_col,
        sfr_col.replace("_mean", "_p16"),
        sfr_col.replace("_mean", "_p50"),
        sfr_col.replace("_mean", "_p84"),
    ]
    keep = [column for column in keep if column in summary]
    joined = joined.merge(summary[keep], on="row", how="left", validate="one_to_one")
    if joined["model"].isna().any():
        rows = joined.loc[joined["model"].isna(), "row"].tolist()
        raise ValueError(f"design rows lack suite summaries: {rows}")

    sfr_history_col = sfr_col.removesuffix("_mean")
    model_summaries = []
    time_series = []
    for row, model_name in zip(joined["row"], joined["model"]):
        model = suite / str(model_name)
        stats, samples = history_diagnostics(
            primary_history(model), bounds, h2_scalar, sfr_history_col
        )
        model_summaries.append(stats)
        samples.insert(0, "row", int(row))
        samples.insert(1, "model", str(model_name))
        time_series.append(samples)
        print(f"  {model.name}: fmol={stats['fmol_mean']:.4f}")
    joined = pd.concat([joined, pd.DataFrame(model_summaries)], axis=1)
    joined["Sigma_SFR_summary_mean"] = joined[sfr_col]
    joined["Sigma_SFR"] = joined["sfr_history_mean"]
    joined["fmol"] = joined["fmol_mean"]
    joined = add_environment_derived(joined)
    temporal = pd.concat(time_series, ignore_index=True)
    temporal = temporal.merge(
        joined[["row", *DESIGN_FIELDS, "kappa", "rho_star"]],
        on="row",
        how="left",
        validate="many_to_one",
    )
    return joined, temporal


def _column_values(table, name: str) -> np.ndarray:
    return np.asarray(table[name], dtype=float)


def load_phangs_reference(
    config: Path,
    sigma_gas: float,
    delta: float,
    obs_sfr_col: str,
) -> pd.DataFrame:
    result = load_configured_phangs(str(config), base_dir=str(ROOT))
    designer = PHANGSSamplingDesigner(
        result["table"],
        config=SamplingConfig(
            target_sigma_gas=sigma_gas,
            delta_sigma_gas=delta,
        ),
    )
    reference = designer.select_reference_pixels()
    required = [*DESIGN_FIELDS, "Sigma_mol", obs_sfr_col]
    missing = [column for column in required if column not in reference.colnames]
    if missing:
        raise KeyError(f"PHANGS reference lacks columns {missing}")
    frame = pd.DataFrame({column: _column_values(reference, column) for column in required})
    frame["Sigma_SFR"] = frame[obs_sfr_col]
    frame["fmol"] = np.divide(
        frame["Sigma_mol"],
        frame["Sigma_gas"],
        out=np.full(len(frame), np.nan),
        where=frame["Sigma_gas"].to_numpy() > 0,
    )
    return add_environment_derived(frame)


def _error_columns(outcome: str, sfr_col: str) -> tuple[str, str] | None:
    if outcome == "fmol":
        return "fmol_p16", "fmol_p84"
    if outcome == "Sigma_SFR":
        return "sfr_history_p16", "sfr_history_p84"
    return None


def plot_outcome(
    reference: pd.DataFrame,
    suite: pd.DataFrame,
    temporal: pd.DataFrame,
    outcome: str,
    sfr_col: str,
    out: Path,
    dpi: int,
    title: str,
    budget_row: pd.Series,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4), sharey=True)
    axes = axes.ravel()
    error_columns = _error_columns(outcome, sfr_col)
    for index, (axis, predictor) in enumerate(zip(axes, PREDICTORS)):
        xref = reference[predictor].to_numpy(dtype=float)
        yref = reference[outcome].to_numpy(dtype=float)
        ref_good = np.isfinite(xref) & np.isfinite(yref) & (xref > 0) & (yref > 0)
        axis.scatter(
            xref[ref_good],
            yref[ref_good],
            s=10,
            color="0.72",
            alpha=0.48,
            linewidth=0,
            rasterized=True,
            zorder=1,
        )

        xtime = temporal[predictor].to_numpy(dtype=float)
        ytime = temporal[outcome].to_numpy(dtype=float)
        time_good = (
            np.isfinite(xtime) & np.isfinite(ytime) & (xtime > 0) & (ytime > 0)
        )
        axis.scatter(
            xtime[time_good],
            ytime[time_good],
            s=1.0,
            color="#3a86ff",
            alpha=0.035,
            linewidth=0,
            rasterized=True,
            zorder=2,
        )

        xsim = suite[predictor].to_numpy(dtype=float)
        ysim = suite[outcome].to_numpy(dtype=float)
        sim_good = np.isfinite(xsim) & np.isfinite(ysim) & (xsim > 0) & (ysim > 0)
        if error_columns and all(column in suite for column in error_columns):
            low = suite[error_columns[0]].to_numpy(dtype=float)
            high = suite[error_columns[1]].to_numpy(dtype=float)
            lower = np.maximum(ysim - low, 0.0)
            upper = np.maximum(high - ysim, 0.0)
            axis.errorbar(
                xsim[sim_good],
                ysim[sim_good],
                yerr=np.vstack([lower[sim_good], upper[sim_good]]),
                fmt="none",
                ecolor="#c52d3b",
                elinewidth=0.8,
                alpha=0.72,
                capsize=1.5,
                zorder=3,
            )
        axis.scatter(
            xsim[sim_good],
            ysim[sim_good],
            s=48,
            color="#e63946",
            edgecolor="black",
            linewidth=0.45,
            zorder=4,
        )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel(LABELS[predictor])
        if index % 3 == 0:
            ylabel = (
                r"$\Sigma_{\rm SFR}$ [$M_\odot\,{\rm kpc}^{-2}\,{\rm yr}^{-1}$]"
                if outcome == "Sigma_SFR"
                else r"$\Sigma_{\rm mol}/\Sigma_{\rm gas}$"
            )
            axis.set_ylabel(ylabel)

    obs_valid = int(np.sum(np.isfinite(reference[outcome]) & (reference[outcome] > 0)))
    sim_valid = int(np.sum(np.isfinite(suite[outcome]) & (suite[outcome] > 0)))
    handles = [
        Line2D(
            [0], [0], marker="o", linestyle="none", markerfacecolor="0.72",
            markeredgecolor="none", markersize=6,
            label=f"PHANGS band ({obs_valid:,} positive)",
        ),
        Line2D(
            [0], [0], marker=".", linestyle="none", color="#3a86ff",
            markersize=7, label="TIGRESS-NCR history samples",
        ),
        Line2D(
            [0], [0], marker="o", linestyle="none", markerfacecolor="#e63946",
            markeredgecolor="black", markeredgewidth=0.45, markersize=7,
            label=f"TIGRESS-NCR time means ({sim_valid} models; 16–84% bars)",
        ),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.965))
    scatter_note = (
        "log scatter [dex]: "
        f"PHANGS {budget_row['sigma_phangs_dex']:.2f}; "
        f"environment {budget_row['sigma_environment_dex']:.2f}; "
        f"temporal {budget_row['sigma_temporal_dex']:.2f}; "
        f"combined suite {budget_row['sigma_suite_total_dex']:.2f}"
    )
    fig.suptitle(f"{title}\n{scatter_note}", y=1.035, fontsize="medium")
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.09, top=0.88, wspace=0.08, hspace=0.24)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print(f"Wrote {out}")


def scatter_budget(
    reference: pd.DataFrame,
    temporal: pd.DataFrame,
    outcomes: tuple[str, ...] = ("Sigma_SFR", "fmol"),
) -> pd.DataFrame:
    """Decompose balanced model-time log variance and compare with PHANGS."""
    rows = []
    for outcome in outcomes:
        observed = reference[outcome].to_numpy(dtype=float)
        observed = np.log10(observed[np.isfinite(observed) & (observed > 0)])
        means = []
        within_variances = []
        sample_counts = []
        for _, group in temporal.groupby("model", sort=True):
            values = group[outcome].to_numpy(dtype=float)
            values = np.log10(values[np.isfinite(values) & (values > 0)])
            if values.size:
                means.append(float(np.mean(values)))
                within_variances.append(float(np.var(values, ddof=0)))
                sample_counts.append(int(values.size))
        environmental_variance = float(np.var(means, ddof=0))
        temporal_variance = float(np.mean(within_variances))
        total_variance = environmental_variance + temporal_variance
        observed_variance = float(np.var(observed, ddof=0))
        rows.append(
            {
                "outcome": outcome,
                "sigma_phangs_dex": np.sqrt(observed_variance),
                "sigma_environment_dex": np.sqrt(environmental_variance),
                "sigma_temporal_dex": np.sqrt(temporal_variance),
                "sigma_suite_total_dex": np.sqrt(total_variance),
                "suite_to_phangs_variance": (
                    total_variance / observed_variance if observed_variance > 0 else np.nan
                ),
                "environment_variance_fraction_of_suite": (
                    environmental_variance / total_variance if total_variance > 0 else np.nan
                ),
                "temporal_variance_fraction_of_suite": (
                    temporal_variance / total_variance if total_variance > 0 else np.nan
                ),
                "n_phangs_positive": int(observed.size),
                "n_models": int(len(means)),
                "n_time_samples": int(sum(sample_counts)),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    if args.time_max <= args.time_min:
        raise SystemExit("--time-max must exceed --time-min")
    sigma_gas = args.sigma_gas or parse_sigma_gas(args.csv.stem)
    if sigma_gas is None:
        raise SystemExit("give --sigma-gas or use Sgas<S> in the design filename")

    print("Reducing suite molecular histories ...")
    suite_table, temporal = build_suite_table(
        pd.read_csv(args.csv),
        pd.read_csv(args.summary),
        args.suite,
        (args.time_min, args.time_max),
        args.h2_scalar,
        args.sfr_col,
    )
    print("Loading PHANGS parent band ...")
    reference = load_phangs_reference(
        args.config, sigma_gas, args.delta, args.obs_sfr_col
    )

    output_dir = args.output_dir or args.csv.parent
    figure_dir = args.figure_dir or output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    joined_path = output_dir / f"{args.csv.stem}_validation_summary.csv"
    suite_table.to_csv(joined_path, index=False)
    print(f"Wrote {joined_path}")
    budget = scatter_budget(reference, temporal)
    budget_path = output_dir / f"{args.csv.stem}_scatter_budget.csv"
    budget.to_csv(budget_path, index=False)
    print(f"Wrote {budget_path}")

    band_title = (
        rf"PHANGS parent band and TIGRESS-NCR core suite: "
        rf"$\Sigma_{{\rm gas}}\approx{sigma_gas:g}$, $\Delta={args.delta:g}$ dex"
    )
    plot_outcome(
        reference,
        suite_table,
        temporal,
        "Sigma_SFR",
        args.sfr_col,
        figure_dir / f"{args.csv.stem}_sigma_sfr_vs_environment.png",
        args.dpi,
        band_title,
        budget.loc[budget["outcome"] == "Sigma_SFR"].iloc[0],
    )
    plot_outcome(
        reference,
        suite_table,
        temporal,
        "fmol",
        args.sfr_col,
        figure_dir / f"{args.csv.stem}_fmol_vs_environment.png",
        args.dpi,
        band_title,
        budget.loc[budget["outcome"] == "fmol"].iloc[0],
    )

    print(
        "Suite medians: "
        f"Sigma_SFR={suite_table['Sigma_SFR'].median():.4g}, "
        f"fmol={suite_table['fmol_mean'].median():.4g}; "
        "PHANGS positive medians: "
        f"Sigma_SFR={reference.loc[reference['Sigma_SFR'] > 0, 'Sigma_SFR'].median():.4g}, "
        f"fmol={reference.loc[reference['fmol'] > 0, 'fmol'].median():.4g}"
    )
    print("Log-scatter budget [dex]:")
    print(budget.to_string(index=False))


if __name__ == "__main__":
    main()
