"""Sampling utilities for PHANGS-conditioned PRFM parameter studies."""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field
from typing import Literal

import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.spatial import cKDTree
from scipy.stats import gaussian_kde, qmc

SynthesisMethod = Literal["kde_lhs", "observed_lhs"]


@dataclass
class SamplingConfig:
    """Configuration for a conditional PHANGS sampling experiment."""

    target_sigma_gas: float = 10.0
    delta_sigma_gas: float = 0.3
    fairness_tol_dex: float = 0.2
    random_seed: int = 27182
    n_repeats: int = 24
    sample_size_grid: list[int] = dc_field(
        default_factory=lambda: [8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]
    )
    primary_fields: list[str] = dc_field(
        default_factory=lambda: ["Sigma_star", "H_star", "Omega_d"]
    )
    gas_fields: list[str] = dc_field(
        default_factory=lambda: ["Sigma_gas", "Sigma_atom", "Sigma_mol"]
    )
    sfr_fields: list[str] = dc_field(
        default_factory=lambda: [
            "Sigma_SFR_HaW4recal",
            "Sigma_SFR_FUVW4recal",
            "Sigma_SFR_Hacorr",
        ]
    )
    comparison_sfr_fields: list[str] = dc_field(
        default_factory=lambda: [
            "Sigma_SFR_HaW4recal",
            "Sigma_SFR_FUVW4recal",
            "Sigma_SFR_Hacorr",
        ]
    )
    id_fields: list[str] = dc_field(
        default_factory=lambda: ["GALAXY", "ID", "RA", "DEC", "r_gal"]
    )
    quantiles: np.ndarray = dc_field(default_factory=lambda: np.arange(0.1, 1.0, 0.1))
    mol_floor: float = 1.0e-3
    fraction_floor: float = 1.0e-4
    kde_oversample_factor: int = 80
    kde_boundary_quantiles: tuple[float, float] = (0.005, 0.995)
    synthesis_method: SynthesisMethod = "kde_lhs"

    def __post_init__(self) -> None:
        if isinstance(self.sfr_fields, str):
            self.sfr_fields = [self.sfr_fields]
        if isinstance(self.comparison_sfr_fields, str):
            self.comparison_sfr_fields = [self.comparison_sfr_fields]


def _column_values(table: Table, col: str) -> np.ndarray:
    return np.asarray(table[col], dtype=float)


def _rank_coordinates(values: np.ndarray) -> np.ndarray:
    """Map each column to empirical rank coordinates in (0, 1)."""
    n, d = values.shape
    ranks = np.empty((n, d), dtype=float)
    for j in range(d):
        order = np.argsort(values[:, j], kind="mergesort")
        ranks[order, j] = (np.arange(n) + 0.5) / n
    return ranks


class PHANGSSamplingDesigner:
    """Design conditional PHANGS parameter samples.

    The class first selects rows within a target ``Sigma_gas`` band, then
    supports either empirical observed-pixel sampling or synthetic sampling
    from a smoothed log-space KDE.
    """

    def __init__(self, table: Table, config: SamplingConfig | None = None):
        self.table = table
        self.config = config or SamplingConfig()

    @property
    def available_sfr_fields(self) -> list[str]:
        return [
            field for field in self.config.sfr_fields if field in self.table.colnames
        ]

    @property
    def available_comparison_sfr_fields(self) -> list[str]:
        return [
            field
            for field in self.config.comparison_sfr_fields
            if field in self.table.colnames
        ]

    @property
    def primary_fit_fields(self) -> list[str]:
        return [*self.config.gas_fields, *self.config.primary_fields]

    @property
    def joint_fit_fields(self) -> list[str]:
        return [
            *self.config.gas_fields,
            *self.config.primary_fields,
            *self.available_sfr_fields,
        ]

    @property
    def visualization_fields(self) -> list[str]:
        fields = [
            *self.config.gas_fields,
            *self.config.primary_fields,
            *self.available_comparison_sfr_fields,
        ]
        return [field for field in fields if field in self.table.colnames]

    @property
    def matrix_base_fields(self) -> list[tuple[str, bool]]:
        """Default PHANGS correlation-matrix fields and log-scale flags."""
        return [
            ("Sigma_gas", True),
            ("Sigma_atom", True),
            ("Sigma_mol", True),
            ("Sigma_star", True),
            ("Sigma_SFR_HaW4recal", True),
            ("Sigma_SFR_FUVW4recal", True),
            ("Sigma_SFR_Hacorr", True),
            ("Omega_d", True),
            ("H_star", True),
        ]

    def with_config(self, **kwargs) -> "PHANGSSamplingDesigner":
        """Return a new designer with selected config fields changed."""
        values = self.config.__dict__.copy()
        values.update(kwargs)
        return type(self)(self.table, SamplingConfig(**values))

    def select_reference_pixels(self) -> Table:
        """Select observed pixels in the configured ``Sigma_gas`` band."""
        selected, _ = self.select_sigma_gas_pixels()
        return selected

    def _sigma_gas_mask(
        self,
        *,
        target_sigma_gas: float | None = None,
        delta_sigma_gas: float | None = None,
    ) -> np.ndarray:
        cfg = self.config
        target_sigma_gas = (
            cfg.target_sigma_gas if target_sigma_gas is None else target_sigma_gas
        )
        delta_sigma_gas = (
            cfg.delta_sigma_gas if delta_sigma_gas is None else delta_sigma_gas
        )
        if target_sigma_gas <= 0:
            raise ValueError("target_sigma_gas must be positive")
        if delta_sigma_gas < 0:
            raise ValueError("delta_sigma_gas must be non-negative")

        required = ["Sigma_gas", "Sigma_atom", "Omega_d", "H_star", "Sigma_star"]
        missing = [field for field in required if field not in self.table.colnames]
        if missing:
            raise KeyError(f"Missing required sampling columns: {missing}")

        valid = np.ones(len(self.table), dtype=bool)
        for field in required:
            values = _column_values(self.table, field)
            valid &= np.isfinite(values) & (values > 0)

        sigma_gas = _column_values(self.table, "Sigma_gas")
        log_gas = np.full(len(self.table), np.nan, dtype=float)
        positive = np.isfinite(sigma_gas) & (sigma_gas > 0)
        log_gas[positive] = np.log10(sigma_gas[positive])
        in_band = np.abs(log_gas - np.log10(target_sigma_gas)) <= delta_sigma_gas
        return valid & in_band

    def select_sigma_gas_pixels(
        self,
        *,
        target_sigma_gas: float | None = None,
        delta_sigma_gas: float | None = None,
        output_fields: list[str] | None = None,
    ) -> tuple[Table, np.ndarray]:
        """Return pixels in a log ``Sigma_gas`` band and their table mask."""
        mask = self._sigma_gas_mask(
            target_sigma_gas=target_sigma_gas,
            delta_sigma_gas=delta_sigma_gas,
        )
        selected = self.table[mask].copy()
        if "Sigma_mol" in selected.colnames:
            mol = _column_values(selected, "Sigma_mol")
            unit = getattr(selected["Sigma_mol"], "unit", None)
            mol = np.where(np.isfinite(mol), mol, 0.0)
            selected["Sigma_mol"] = mol if unit is None else mol * unit

        if output_fields is None:
            output_fields = [
                *self.config.id_fields,
                *self.config.gas_fields,
                *self.config.primary_fields,
                *self.available_sfr_fields,
                *self.available_comparison_sfr_fields,
            ]
        output_cols = []
        for field in output_fields:
            if field in selected.colnames and field not in output_cols:
                output_cols.append(field)
            err_field = f"e_{field}"
            if err_field in selected.colnames and err_field not in output_cols:
                output_cols.append(err_field)
        return selected[output_cols], mask

    def select_sigma_gas_targets(
        self,
        targets: list[float] | np.ndarray,
        *,
        delta_sigma_gas: float | None = None,
    ) -> dict[float, Table]:
        """Return selected-pixel tables for several target ``Sigma_gas`` values."""
        return {
            float(target): self.select_sigma_gas_pixels(
                target_sigma_gas=float(target),
                delta_sigma_gas=delta_sigma_gas,
            )[0]
            for target in targets
        }

    def resolve_matrix_columns(
        self,
        *,
        aperture: str = "hexagon",
        base_fields: list[tuple[str, bool]] | None = None,
    ) -> list[tuple[str, str, bool]]:
        """Resolve base matrix fields to concrete columns present in the table."""
        from prfm.phangs_plot import col_label, resolve_columns

        base_fields = base_fields or self.matrix_base_fields
        spec = [(field, col_label(field), log) for field, log in base_fields]
        return resolve_columns(spec, self.table, aperture=aperture)

    def plot_correlation_matrix(
        self,
        *,
        targets: list[float] | np.ndarray | None = None,
        delta_sigma_gas: float | None = None,
        cols: list[tuple[str, str, bool]] | None = None,
        aperture: str = "hexagon",
        table_label: str = "PHANGS",
        figsize_scale: float = 2.2,
    ):
        """Plot the PHANGS correlation matrix with optional Sigma-gas-band overlays."""
        import matplotlib.pyplot as plt
        from prfm.phangs_plot import col_label, hist_plot, scatter_plot

        cols = cols or self.resolve_matrix_columns(aperture=aperture)
        col_names = [col for col, _, _ in cols]
        col_logs = {col: log for col, _, log in cols}
        targets = (
            []
            if targets is None
            else list(np.atleast_1d(np.asarray(targets, dtype=float)))
        )
        colors = plt.cm.tab10(np.linspace(0, 1, max(len(targets), 1)))
        selections = [
            (
                target,
                self.select_sigma_gas_pixels(
                    target_sigma_gas=target,
                    delta_sigma_gas=delta_sigma_gas,
                )[1],
                color,
            )
            for target, color in zip(targets, colors)
        ]

        n = len(col_names)
        fig, axes = plt.subplots(n, n, figsize=(n * figsize_scale, n * figsize_scale))
        for i, ycol in enumerate(col_names):
            for j, xcol in enumerate(col_names):
                ax = axes[i, j]
                if i == j:
                    valid = np.isfinite(_column_values(self.table, ycol))
                    hist_plot(self.table, ycol, ax=ax, log=col_logs[ycol])
                    for target, mask, color in selections:
                        self._plot_selected_hist(
                            ax,
                            mask,
                            ycol,
                            log=col_logs[ycol],
                            color=color,
                            label=rf"$\Sigma_{{gas}}={target:g}$",
                        )
                    ax.set_ylabel("PDF")
                else:
                    valid = np.isfinite(_column_values(self.table, xcol)) & np.isfinite(
                        _column_values(self.table, ycol)
                    )
                    scatter_plot(
                        self.table,
                        xcol,
                        ycol,
                        ax=ax,
                        log_x=col_logs[xcol],
                        log_y=col_logs[ycol],
                        errorbars=True,
                        s=2,
                        bg_alpha=0.18 if selections else 0.3,
                    )
                    for target, mask, color in selections:
                        self._plot_selected_scatter(
                            ax,
                            mask,
                            xcol,
                            ycol,
                            color=color,
                            label=rf"$\Sigma_{{gas}}={target:g}$",
                        )
                    if ycol.startswith("Sigma_SFR") & xcol.startswith("Sigma_SFR"):
                        sfr = np.logspace(-5, 0, 100)
                        ax.plot(sfr, sfr, color="black", linestyle="--", linewidth=1)

                ax.annotate(
                    f"N={valid.sum()}",
                    xy=(0.97, 0.97),
                    xycoords="axes fraction",
                    ha="right",
                    va="top",
                    fontsize="xx-small",
                    color="0.4",
                )
                ax.tick_params(labelsize="x-small")
                if i < n - 1:
                    ax.set_xlabel("")
                    ax.tick_params(labelbottom=False)
                else:
                    ax.set_xlabel(col_label(xcol))
                if j > 0:
                    ax.set_ylabel("")
                    ax.tick_params(labelleft=False)

        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(
                handles,
                labels,
                loc="upper right",
                bbox_to_anchor=(1.0, 1.0),
                fontsize="small",
            )
        title = f"{table_label} {aperture} aperture"
        if targets:
            delta = (
                self.config.delta_sigma_gas
                if delta_sigma_gas is None
                else delta_sigma_gas
            )
            title += f" - Sigma_gas sampling (delta={delta:g} dex)"
        else:
            title += f" - correlation matrix ({len(self.table)} rows)"
        fig.suptitle(title, fontsize="large", y=1.002)
        fig.tight_layout(h_pad=0.3, w_pad=0.3)
        return fig, axes, selections

    def _plot_selected_hist(
        self,
        ax,
        mask: np.ndarray,
        col: str,
        *,
        log: bool = True,
        color="tab:blue",
        label: str | None = None,
        n_bins: int = 40,
    ) -> None:
        all_vals = _column_values(self.table, col)
        all_valid = (
            np.isfinite(all_vals) & (all_vals > 0) if log else np.isfinite(all_vals)
        )
        selected = all_vals[mask & all_valid]
        all_vals = all_vals[all_valid]
        if len(selected) == 0 or len(all_vals) == 0:
            return
        if log:
            edges = np.logspace(
                np.log10(all_vals.min()), np.log10(all_vals.max()), n_bins + 1
            )
            widths = np.diff(np.log10(edges))
        else:
            edges = np.linspace(all_vals.min(), all_vals.max(), n_bins + 1)
            widths = np.diff(edges)
        counts, _ = np.histogram(selected, bins=edges)
        pdf = counts / len(selected) / widths
        ax.step(
            edges[:-1],
            pdf,
            where="post",
            color=color,
            linewidth=1.2,
            label=label,
            zorder=3,
        )

    def _plot_selected_scatter(
        self,
        ax,
        mask: np.ndarray,
        xcol: str,
        ycol: str,
        *,
        color="tab:blue",
        label: str | None = None,
        s: float = 12,
    ) -> None:
        x = _column_values(self.table, xcol)
        y = _column_values(self.table, ycol)
        valid = mask & np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        if not valid.any():
            return
        ax.scatter(
            x[valid],
            y[valid],
            s=s,
            color=color,
            alpha=0.85,
            linewidths=0,
            label=label,
            rasterized=True,
            zorder=5,
        )

    def _log_dataframe(
        self,
        table: Table,
        fields: list[str],
        *,
        fill_mol_floor: bool = True,
    ) -> pd.DataFrame:
        columns = []
        valid = np.ones(len(table), dtype=bool)
        for field in fields:
            values = _column_values(table, field).copy()
            if fill_mol_floor and field == "Sigma_mol":
                values = np.where(
                    np.isfinite(values) & (values > 0), values, self.config.mol_floor
                )
                valid &= np.isfinite(values) & (values >= self.config.mol_floor)
            else:
                valid &= np.isfinite(values) & (values > 0)
            columns.append(values)
        values = np.column_stack(columns)
        return pd.DataFrame(np.log10(values[valid]), columns=fields)

    def _kde_fit_fields(self, fields: list[str]) -> list[str]:
        """Use total gas plus molecular fraction for constrained gas sampling."""
        if all(field in fields for field in self.config.gas_fields):
            return [
                *[
                    field
                    for field in fields
                    if field not in ("Sigma_atom", "Sigma_mol")
                ],
                "f_mol",
            ]
        return fields

    def _kde_log_dataframe(self, table: Table, fields: list[str]) -> pd.DataFrame:
        fit_fields = self._kde_fit_fields(fields)
        columns = {}
        valid = np.ones(len(table), dtype=bool)

        if "f_mol" in fit_fields:
            sigma_gas = _column_values(table, "Sigma_gas")
            sigma_mol = (
                _column_values(table, "Sigma_mol")
                if "Sigma_mol" in table.colnames
                else np.zeros(len(table))
            )
            sigma_mol = np.where(
                np.isfinite(sigma_mol) & (sigma_mol > 0), sigma_mol, 0.0
            )
            valid &= np.isfinite(sigma_gas) & (sigma_gas > 0)
            f_mol = np.divide(
                sigma_mol,
                sigma_gas,
                out=np.zeros_like(sigma_gas, dtype=float),
                where=np.isfinite(sigma_gas) & (sigma_gas > 0),
            )
            f_mol = np.clip(
                f_mol, self.config.fraction_floor, 1.0 - self.config.fraction_floor
            )
            columns["f_mol"] = np.log10(f_mol / (1.0 - f_mol))

        for field in fit_fields:
            if field == "f_mol":
                continue
            values = _column_values(table, field).copy()
            positive = np.isfinite(values) & (values > 0)
            valid &= positive
            log_values = np.full(len(table), np.nan, dtype=float)
            log_values[positive] = np.log10(values[positive])
            columns[field] = log_values

        return pd.DataFrame(
            {field: values[valid] for field, values in columns.items()},
            columns=fit_fields,
        )

    def _kde_candidate_pool(
        self,
        log_df: pd.DataFrame,
        n_candidates: int,
        seed: int,
    ) -> pd.DataFrame:
        rng = np.random.default_rng(seed)
        data = log_df.to_numpy(dtype=float)
        kde = gaussian_kde(data.T)
        lo = np.quantile(data, self.config.kde_boundary_quantiles[0], axis=0)
        hi = np.quantile(data, self.config.kde_boundary_quantiles[1], axis=0)

        chunks = []
        attempts = 0
        while sum(len(chunk) for chunk in chunks) < n_candidates and attempts < 50:
            needed = n_candidates - sum(len(chunk) for chunk in chunks)
            draw = kde.resample(max(n_candidates, 4 * needed), seed=rng).T
            chunks.append(draw[np.all((draw >= lo) & (draw <= hi), axis=1)])
            attempts += 1

        if not chunks or sum(len(chunk) for chunk in chunks) < n_candidates:
            raise RuntimeError("KDE candidate generation failed")
        return pd.DataFrame(np.vstack(chunks)[:n_candidates], columns=log_df.columns)

    def _lhs_select_log_candidates(
        self,
        candidates: pd.DataFrame,
        lhs_fields: list[str],
        n_samples: int,
        seed: int,
    ) -> pd.DataFrame:
        rank_coords = _rank_coordinates(candidates[lhs_fields].to_numpy(dtype=float))
        lhs = qmc.LatinHypercube(d=len(lhs_fields), seed=seed).random(n=n_samples)
        tree = cKDTree(rank_coords)
        k = min(len(candidates), max(16, 5 * n_samples))
        _, neighbor_positions = tree.query(lhs, k=k)
        if k == 1:
            neighbor_positions = neighbor_positions[:, None]

        chosen = []
        chosen_set = set()
        for row in neighbor_positions:
            for pos in np.atleast_1d(row):
                pos = int(pos)
                if pos not in chosen_set:
                    chosen.append(pos)
                    chosen_set.add(pos)
                    break
        if len(chosen) < n_samples:
            remaining = [i for i in range(len(candidates)) if i not in chosen_set]
            chosen.extend(remaining[: n_samples - len(chosen)])
        return candidates.iloc[chosen[:n_samples]].reset_index(drop=True)

    def synthesize_kde_lhs(
        self,
        reference: Table,
        n_samples: int,
        *,
        fit_fields: list[str] | None = None,
        lhs_fields: list[str] | None = None,
        seed: int | None = None,
    ) -> pd.DataFrame:
        """Draw a synthetic LHS sample from a log-space KDE."""
        cfg = self.config
        fit_fields = fit_fields or self.joint_fit_fields
        lhs_fields = lhs_fields or cfg.primary_fields
        seed = cfg.random_seed if seed is None else seed

        reference_log = self._kde_log_dataframe(reference, fit_fields)
        if len(reference_log) < max(20, 2 * len(fit_fields)):
            raise ValueError(
                f"Only {len(reference_log)} complete rows for {fit_fields}; KDE is underconstrained"
            )

        n_candidates = max(cfg.kde_oversample_factor * n_samples, 2000)
        candidates = self._kde_candidate_pool(reference_log, n_candidates, seed)
        sample_log = self._lhs_select_log_candidates(
            candidates, lhs_fields, n_samples, seed + 1
        )
        sample = self._linearize_kde_sample(sample_log)
        sample.attrs["log10_sample"] = sample_log
        sample.attrs["log10_reference"] = reference_log
        return sample

    def _linearize_kde_sample(self, sample_log: pd.DataFrame) -> pd.DataFrame:
        sample = pd.DataFrame(index=sample_log.index)
        for field in sample_log.columns:
            if field != "f_mol":
                sample[field] = 10.0 ** sample_log[field]

        if "f_mol" in sample_log.columns and "Sigma_gas" in sample.columns:
            odds = 10.0 ** sample_log["f_mol"]
            f_mol = odds / (1.0 + odds)
            f_mol = np.clip(f_mol, 0.0, 1.0)
            sample["Sigma_mol"] = sample["Sigma_gas"] * f_mol
            sample["Sigma_atom"] = sample["Sigma_gas"] - sample["Sigma_mol"]
            floor = self.config.mol_floor
            mol_floor_mask = sample["Sigma_mol"] <= 1.05 * floor
            sample.loc[mol_floor_mask, "Sigma_mol"] = 0.0
            sample.loc[mol_floor_mask, "Sigma_atom"] = sample.loc[
                mol_floor_mask, "Sigma_gas"
            ]

        ordered = []
        for field in [
            *self.config.gas_fields,
            *self.config.primary_fields,
            *self.available_sfr_fields,
        ]:
            if field in sample.columns and field not in ordered:
                ordered.append(field)
        ordered.extend([field for field in sample.columns if field not in ordered])
        return sample[ordered]

    def sample_observed_lhs(
        self,
        reference: Table,
        n_samples: int,
        *,
        lhs_fields: list[str] | None = None,
        seed: int | None = None,
    ) -> Table:
        """Select observed pixels closest to an LHS design in rank space."""
        cfg = self.config
        lhs_fields = lhs_fields or cfg.primary_fields
        seed = cfg.random_seed if seed is None else seed
        log_df = self._log_dataframe(reference, lhs_fields, fill_mol_floor=False)
        values = log_df.to_numpy(dtype=float)
        valid = np.ones(len(reference), dtype=bool)
        for field in lhs_fields:
            x = _column_values(reference, field)
            valid &= np.isfinite(x) & (x > 0)
        candidate_indices = np.flatnonzero(valid)
        n_available = len(candidate_indices)
        if n_samples > n_available:
            raise ValueError(
                f"Requested {n_samples} samples, but only {n_available} valid rows are available"
            )

        rank_coords = _rank_coordinates(values)
        lhs = qmc.LatinHypercube(d=len(lhs_fields), seed=seed).random(n=n_samples)
        tree = cKDTree(rank_coords)
        k = min(n_available, max(16, 4 * n_samples))
        _, neighbors = tree.query(lhs, k=k)
        if k == 1:
            neighbors = neighbors[:, None]

        chosen = []
        chosen_set = set()
        for row in neighbors:
            for pos in np.atleast_1d(row):
                pos = int(pos)
                if pos not in chosen_set:
                    chosen.append(pos)
                    chosen_set.add(pos)
                    break
        if len(chosen) < n_samples:
            remaining = [i for i in range(n_available) if i not in chosen_set]
            chosen.extend(remaining[: n_samples - len(chosen)])
        return reference[candidate_indices[np.asarray(chosen[:n_samples], dtype=int)]]

    def sample(
        self,
        reference: Table,
        n_samples: int,
        *,
        method: SynthesisMethod | None = None,
        fit_fields: list[str] | None = None,
        lhs_fields: list[str] | None = None,
        seed: int | None = None,
    ) -> pd.DataFrame | Table:
        """Sample with the configured or requested synthesis method."""
        method = method or self.config.synthesis_method
        if method == "kde_lhs":
            return self.synthesize_kde_lhs(
                reference,
                n_samples,
                fit_fields=fit_fields,
                lhs_fields=lhs_fields,
                seed=seed,
            )
        if method == "observed_lhs":
            return self.sample_observed_lhs(
                reference, n_samples, lhs_fields=lhs_fields, seed=seed
            )
        raise ValueError(f"Unknown synthesis method: {method}")

    def assign_sfr_from_neighbors(
        self,
        reference: Table,
        sample: pd.DataFrame,
        *,
        match_fields: list[str] | None = None,
        sfr_fields: list[str] | None = None,
        k_neighbors: int = 16,
        seed: int | None = None,
    ) -> pd.DataFrame:
        """Attach SFR values to a synthetic sample from nearby observed pixels.

        This is useful for primary-only synthetic samples. The sample is matched
        to observed pixels in log-space ``match_fields`` and SFR values are
        copied from one of the nearest observed neighbors. Neighbor choice is
        distance-weighted, so repeated calls can propagate local SFR scatter.
        """
        cfg = self.config
        match_fields = match_fields or self.primary_fit_fields
        sfr_fields = sfr_fields or self.available_sfr_fields
        seed = cfg.random_seed if seed is None else seed
        rng = np.random.default_rng(seed)

        missing_sample = [
            field for field in match_fields if field not in sample.columns
        ]
        if missing_sample:
            raise KeyError(f"Sample is missing match fields: {missing_sample}")
        missing_reference = [
            field
            for field in [*match_fields, *sfr_fields]
            if field not in reference.colnames
        ]
        if missing_reference:
            raise KeyError(f"Reference table is missing fields: {missing_reference}")

        ref_log = self._log_dataframe(reference, [*match_fields, *sfr_fields])
        # Recompute the mask explicitly so row indices in ref_log align to reference rows.
        valid_ref = np.ones(len(reference), dtype=bool)
        for field in [*match_fields, *sfr_fields]:
            values = _column_values(reference, field)
            if field == "Sigma_mol":
                values = np.where(
                    np.isfinite(values) & (values > 0), values, cfg.mol_floor
                )
                valid_ref &= np.isfinite(values) & (values >= cfg.mol_floor)
            else:
                valid_ref &= np.isfinite(values) & (values > 0)
        ref_indices = np.flatnonzero(valid_ref)

        sample_log = self._as_log_dataframe(sample, match_fields)
        if len(sample_log) != len(sample):
            raise ValueError(
                "All synthetic sample rows must be finite and positive in match_fields"
            )

        tree = cKDTree(ref_log[match_fields].to_numpy(dtype=float))
        k = min(k_neighbors, len(ref_log))
        distances, neighbor_pos = tree.query(
            sample_log[match_fields].to_numpy(dtype=float), k=k
        )
        if k == 1:
            distances = distances[:, None]
            neighbor_pos = neighbor_pos[:, None]

        out = sample.copy()
        chosen_reference_indices = []
        for i in range(len(out)):
            d = np.asarray(distances[i], dtype=float)
            pos = np.asarray(neighbor_pos[i], dtype=int)
            weights = 1.0 / np.maximum(d, 1.0e-6) ** 2
            weights = weights / weights.sum()
            chosen_pos = int(rng.choice(pos, p=weights))
            ref_idx = int(ref_indices[chosen_pos])
            chosen_reference_indices.append(ref_idx)
            for field in sfr_fields:
                out.loc[out.index[i], field] = float(reference[field][ref_idx])

        out.attrs.update(sample.attrs)
        out.attrs["sfr_assignment_method"] = "distance_weighted_nearest_neighbors"
        out.attrs["sfr_assignment_match_fields"] = match_fields
        out.attrs["sfr_assignment_k_neighbors"] = k
        out.attrs["sfr_assignment_reference_indices"] = np.asarray(
            chosen_reference_indices, dtype=int
        )
        return out

    def quantile_error_report(
        self,
        reference: Table | pd.DataFrame,
        sample: Table | pd.DataFrame,
        fields: list[str],
    ) -> tuple[float, pd.DataFrame]:
        """Return worst log10 quantile error and per-field report."""
        ref_log = self._as_log_dataframe(reference, fields)
        sample_log = self._as_log_dataframe(sample, fields)
        rows = []
        max_error = 0.0
        for field in fields:
            if field not in ref_log.columns or field not in sample_log.columns:
                continue
            ref = ref_log[field].to_numpy(dtype=float)
            draw = sample_log[field].to_numpy(dtype=float)
            if len(ref) == 0 or len(draw) == 0:
                error = np.nan
            else:
                error = float(
                    np.max(
                        np.abs(
                            np.quantile(draw, self.config.quantiles)
                            - np.quantile(ref, self.config.quantiles)
                        )
                    )
                )
                max_error = max(max_error, error)
            rows.append(
                {
                    "field": field,
                    "n_ref": len(ref),
                    "n_sample": len(draw),
                    "max_quantile_error_dex": error,
                    "fair": bool(
                        np.isfinite(error) and error <= self.config.fairness_tol_dex
                    ),
                }
            )
        return max_error, pd.DataFrame(rows)

    def quantile_error_report_against_sample_field(
        self,
        reference: Table | pd.DataFrame,
        sample: Table | pd.DataFrame,
        reference_fields: list[str],
        sample_field: str,
    ) -> tuple[float, pd.DataFrame]:
        """Compare several reference fields against one sampled field."""
        max_error = 0.0
        rows = []
        sample_log = self._as_log_dataframe(sample, [sample_field])
        if sample_field not in sample_log.columns:
            raise KeyError(f"Sample is missing field {sample_field}")
        draw = sample_log[sample_field].to_numpy(dtype=float)
        for field in reference_fields:
            ref_log = self._as_log_dataframe(reference, [field])
            if field not in ref_log.columns or len(draw) == 0:
                error = np.nan
                n_ref = 0
            else:
                ref = ref_log[field].to_numpy(dtype=float)
                n_ref = len(ref)
                error = float(
                    np.max(
                        np.abs(
                            np.quantile(draw, self.config.quantiles)
                            - np.quantile(ref, self.config.quantiles)
                        )
                    )
                )
                max_error = max(max_error, error)
            rows.append(
                {
                    "reference_field": field,
                    "sample_field": sample_field,
                    "n_ref": n_ref,
                    "n_sample": len(draw),
                    "max_quantile_error_dex": error,
                    "fair": bool(
                        np.isfinite(error) and error <= self.config.fairness_tol_dex
                    ),
                }
            )
        return max_error, pd.DataFrame(rows)

    def _as_log_dataframe(
        self, data: Table | pd.DataFrame, fields: list[str]
    ) -> pd.DataFrame:
        if isinstance(data, Table):
            table_fields = [field for field in fields if field in data.colnames]
            if not table_fields:
                return pd.DataFrame()
            return self._log_dataframe(data, table_fields)

        columns = {}
        for field in fields:
            if field not in data.columns:
                continue
            values = np.asarray(data[field], dtype=float)
            if field == "Sigma_mol":
                values = np.where(
                    np.isfinite(values) & (values > 0), values, self.config.mol_floor
                )
            valid = np.isfinite(values) & (values > 0)
            columns[field] = np.log10(values[valid])
        return pd.DataFrame(columns)

    def estimate_required_samples(
        self,
        reference: Table,
        *,
        method: SynthesisMethod | None = None,
        fit_fields: list[str] | None = None,
        lhs_fields: list[str] | None = None,
        fairness_fields: list[str] | None = None,
    ) -> tuple[int | None, pd.DataFrame]:
        """Scan sample sizes and return the first meeting the fairness criterion."""
        cfg = self.config
        method = method or cfg.synthesis_method
        lhs_fields = lhs_fields or cfg.primary_fields
        fairness_fields = fairness_fields or lhs_fields
        rows = []
        required = None
        sizes = list(cfg.sample_size_grid)
        if method == "observed_lhs":
            sizes = sorted({n for n in sizes if n <= len(reference)} | {len(reference)})

        for n_samples in sizes:
            errors = []
            for repeat in range(cfg.n_repeats):
                draw = self.sample(
                    reference,
                    n_samples,
                    method=method,
                    fit_fields=fit_fields,
                    lhs_fields=lhs_fields,
                    seed=cfg.random_seed + 1000 * n_samples + repeat,
                )
                error, _ = self.quantile_error_report(reference, draw, fairness_fields)
                errors.append(error)
            errors = np.asarray(errors, dtype=float)
            p90_error = float(np.nanquantile(errors, 0.9))
            median_error = float(np.nanmedian(errors))
            fair = p90_error <= cfg.fairness_tol_dex
            rows.append(
                {
                    "n_samples": n_samples,
                    "median_max_quantile_error_dex": median_error,
                    "p90_max_quantile_error_dex": p90_error,
                    "fair": fair,
                }
            )
            if fair and required is None:
                required = n_samples
        return required, pd.DataFrame(rows)

    def plot_distribution_overlay(
        self,
        reference: Table,
        sample: Table | pd.DataFrame,
        *,
        fields: list[str] | None = None,
        bins: int = 36,
        bin_source: Literal["full", "reference", "combined"] = "full",
        show_kde: bool = False,
        kde_n_draws: int = 10000,
        kde_seed: int | None = None,
        sample_sfr_field: str | None = None,
        comparison_sfr_fields: list[str] | None = None,
    ):
        """Plot 1-D observed-vs-sample distribution overlays."""
        import matplotlib.pyplot as plt
        from prfm.phangs_plot import col_label

        fields = fields or self.visualization_fields
        comparison_sfr_fields = comparison_sfr_fields or []
        fields = self._available_plot_fields(
            reference,
            sample,
            fields,
            sample_sfr_field=sample_sfr_field,
            comparison_sfr_fields=comparison_sfr_fields,
        )
        kde_draw = (
            self._sample_kde_overlay(sample, n_draws=kde_n_draws, seed=kde_seed)
            if show_kde
            else None
        )
        ncols = 3
        nrows = int(np.ceil(len(fields) / ncols))
        fig, axes = plt.subplots(
            nrows, ncols, figsize=(4.0 * ncols, 3.0 * nrows), squeeze=False
        )

        for ax, field in zip(axes.ravel(), fields):
            sample_field = self._plot_sample_field(
                field,
                sample,
                sample_sfr_field=sample_sfr_field,
                comparison_sfr_fields=comparison_sfr_fields,
            )
            ref = self._positive_plot_values(reference, field)
            draw = self._positive_plot_values(sample, sample_field)
            if len(ref) == 0 or len(draw) == 0:
                ax.set_axis_off()
                continue
            edges = self._hist_edges(field, ref, draw, bins=bins, bin_source=bin_source)
            ax.hist(
                ref,
                bins=edges,
                histtype="stepfilled",
                density=True,
                alpha=0.22,
                color="0.55",
                label="observed",
            )
            sample_label = (
                "sample" if sample_field == field else f"sample {sample_field}"
            )
            ax.hist(
                draw,
                bins=edges,
                histtype="step",
                density=True,
                linewidth=1.8,
                color="tab:red",
                label=sample_label,
            )
            if kde_draw is not None and sample_field in kde_draw.columns:
                kde_vals = self._positive_plot_values(kde_draw, sample_field)
                if len(kde_vals):
                    ax.hist(
                        kde_vals,
                        bins=edges,
                        histtype="step",
                        density=True,
                        linewidth=1.4,
                        color="tab:blue",
                        linestyle="--",
                        label="KDE model",
                    )
            ax.set_xscale("log")
            ax.set_xlabel(col_label(field))
            ax.set_ylabel("PDF")
            title = f"{field}: Nobs={len(ref)}, Nsamp={len(draw)}"
            if sample_field != field:
                title += f"\n(sample={sample_field})"
            ax.set_title(title, fontsize="small")

        for ax in axes.ravel()[len(fields) :]:
            ax.set_axis_off()
        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper right")
        fig.suptitle("Observed conditional distribution vs sample", y=1.01)
        fig.tight_layout()
        return fig, axes

    def _sample_kde_overlay(
        self,
        sample: Table | pd.DataFrame,
        *,
        n_draws: int,
        seed: int | None,
    ) -> pd.DataFrame | None:
        if not isinstance(sample, pd.DataFrame):
            return None
        reference_log = sample.attrs.get("log10_reference")
        if reference_log is None:
            return None
        seed = self.config.random_seed + 99173 if seed is None else seed
        candidates = self._kde_candidate_pool(reference_log, max(n_draws, 2000), seed)
        return self._linearize_kde_sample(
            candidates.iloc[:n_draws].reset_index(drop=True)
        )

    def plot_pair_overlay(
        self,
        reference: Table,
        sample: Table | pd.DataFrame,
        *,
        fields: list[str] | None = None,
    ):
        """Plot pairwise observed-vs-sample overlays."""
        import matplotlib.pyplot as plt
        from prfm.phangs_plot import col_label

        fields = fields or self.visualization_fields
        fields = self._available_plot_fields(reference, sample, fields)
        n = len(fields)
        fig, axes = plt.subplots(n, n, figsize=(1.75 * n, 1.75 * n), squeeze=False)
        for i, yfield in enumerate(fields):
            for j, xfield in enumerate(fields):
                ax = axes[i, j]
                if i == j:
                    ref = self._positive_plot_values(reference, xfield)
                    draw = self._positive_plot_values(sample, xfield)
                    if len(ref) and len(draw):
                        edges = self._hist_edges(
                            xfield, ref, draw, bins=24, bin_source="full"
                        )
                        ax.hist(
                            ref,
                            bins=edges,
                            histtype="stepfilled",
                            density=True,
                            alpha=0.18,
                            color="0.55",
                        )
                        ax.hist(
                            draw,
                            bins=edges,
                            histtype="step",
                            density=True,
                            linewidth=1.0,
                            color="tab:red",
                        )
                        ax.set_xscale("log")
                else:
                    x_ref = self._raw_values(reference, xfield)
                    y_ref = self._raw_values(reference, yfield)
                    x_draw = self._raw_values(sample, xfield)
                    y_draw = self._raw_values(sample, yfield)
                    ref_valid = self._plot_pair_mask(x_ref, y_ref, xfield, yfield)
                    draw_valid = self._plot_pair_mask(x_draw, y_draw, xfield, yfield)
                    ax.scatter(
                        x_ref[ref_valid],
                        y_ref[ref_valid],
                        s=2,
                        color="0.65",
                        alpha=0.18,
                        linewidths=0,
                        rasterized=True,
                    )
                    ax.scatter(
                        x_draw[draw_valid],
                        y_draw[draw_valid],
                        s=16,
                        color="tab:red",
                        alpha=0.85,
                        linewidths=0,
                        rasterized=True,
                    )
                    ax.set_xscale("log")
                    ax.set_yscale("log")

                if i < n - 1:
                    ax.tick_params(axis="x", which="both", labelbottom=False)
                    ax.set_xlabel("")
                else:
                    ax.set_xlabel(col_label(xfield), fontsize="x-small")
                if j > 0:
                    ax.tick_params(axis="y", which="both", labelleft=False)
                    ax.set_ylabel("")
                else:
                    ax.set_ylabel(col_label(yfield), fontsize="x-small")
                ax.tick_params(labelsize="xx-small")
        fig.suptitle("Observed vs sample pair distributions", y=1.002)
        fig.tight_layout(h_pad=0.2, w_pad=0.2)
        return fig, axes

    def _available_plot_fields(
        self,
        reference: Table,
        sample: Table | pd.DataFrame,
        fields: list[str],
        *,
        sample_sfr_field: str | None = None,
        comparison_sfr_fields: list[str] | None = None,
    ) -> list[str]:
        sample_cols = sample.colnames if isinstance(sample, Table) else sample.columns
        comparison_sfr_fields = comparison_sfr_fields or []
        available = []
        for field in fields:
            if field not in reference.colnames:
                continue
            sample_field = self._plot_sample_field(
                field,
                sample,
                sample_sfr_field=sample_sfr_field,
                comparison_sfr_fields=comparison_sfr_fields,
            )
            if sample_field in sample_cols:
                available.append(field)
        return available

    def _plot_sample_field(
        self,
        field: str,
        sample: Table | pd.DataFrame,
        *,
        sample_sfr_field: str | None,
        comparison_sfr_fields: list[str],
    ) -> str:
        sample_cols = sample.colnames if isinstance(sample, Table) else sample.columns
        if (
            sample_sfr_field is not None
            and field in comparison_sfr_fields
            and sample_sfr_field in sample_cols
        ):
            return sample_sfr_field
        return field

    def _raw_values(self, data: Table | pd.DataFrame, field: str) -> np.ndarray:
        return (
            _column_values(data, field)
            if isinstance(data, Table)
            else np.asarray(data[field], dtype=float)
        )

    def _positive_plot_values(
        self, data: Table | pd.DataFrame, field: str
    ) -> np.ndarray:
        values = self._raw_values(data, field)
        threshold = self.config.mol_floor * 1.05 if field == "Sigma_mol" else 0.0
        return values[np.isfinite(values) & (values > threshold)]

    def _hist_edges(
        self,
        field: str,
        reference_values: np.ndarray,
        sample_values: np.ndarray,
        *,
        bins: int,
        bin_source: Literal["full", "reference", "combined"],
    ) -> np.ndarray:
        if bin_source == "full" and field in self.table.colnames:
            values = self._positive_plot_values(self.table, field)
        elif bin_source == "reference":
            values = reference_values
        elif bin_source == "combined":
            values = np.concatenate([reference_values, sample_values])
        else:
            raise ValueError("bin_source must be 'full', 'reference', or 'combined'")

        values = values[np.isfinite(values) & (values > 0)]
        if len(values) == 0:
            values = np.concatenate([reference_values, sample_values])
            values = values[np.isfinite(values) & (values > 0)]
        if len(values) == 0:
            raise ValueError(f"No positive values available to bin {field}")
        return np.logspace(np.log10(values.min()), np.log10(values.max()), bins + 1)

    def _plot_pair_mask(
        self,
        x: np.ndarray,
        y: np.ndarray,
        xfield: str,
        yfield: str,
    ) -> np.ndarray:
        mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        if xfield == "Sigma_mol":
            mask &= x > self.config.mol_floor * 1.05
        if yfield == "Sigma_mol":
            mask &= y > self.config.mol_floor * 1.05
        return mask
