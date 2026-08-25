from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "project" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import plot_validation_diagnostics as diagnostics  # noqa: E402


def test_add_environment_derived():
    frame = pd.DataFrame(
        {
            "Sigma_star": [42.0],
            "H_star": [245.0],
            "Omega": [28.0],
            "qshear": [1.0],
        }
    )
    result = diagnostics.add_environment_derived(frame)
    assert result.loc[0, "kappa"] == pytest.approx(np.sqrt(2.0) * 28.0)
    assert result.loc[0, "rho_star"] == pytest.approx(42.0 / 490.0)


def test_window_series_interpolates_bounds_and_keeps_last_duplicate():
    time, values = diagnostics._window_series(
        np.array([0.0, 1.0, 1.0, 2.0]),
        np.array([0.0, 1.0, 3.0, 4.0]),
        (0.5, 1.5),
    )
    np.testing.assert_allclose(time, [0.5, 1.0, 1.5])
    np.testing.assert_allclose(values, [1.5, 3.0, 3.5])


def test_history_diagnostics_reduces_full_window(monkeypatch, tmp_path):
    history = {
        "time": np.array([0.0, 1.0, 2.0]),
        "mass": np.ones(3),
        "scalar3": np.array([0.1, 0.2, 0.3]),
        "sfr40": np.array([1.0, 2.0, 3.0]),
    }
    monkeypatch.setattr(diagnostics, "read_hst", lambda _: history)
    summary, samples = diagnostics.history_diagnostics(
        tmp_path / "model.hst", (0.5, 1.5), "scalar3", "sfr40"
    )
    assert summary["fmol_mean"] == pytest.approx(0.4)
    assert summary["sfr_history_mean"] == pytest.approx(2.0)
    np.testing.assert_allclose(samples["fmol"], [0.3, 0.4, 0.5])
    np.testing.assert_allclose(samples["Sigma_SFR"], [1.5, 2.0, 2.5])


def test_scatter_budget_obeys_variance_decomposition():
    root_two = np.sqrt(2.0)
    reference = pd.DataFrame({"Sigma_SFR": 10.0 ** np.array([-root_two, root_two])})
    temporal = pd.DataFrame(
        {
            "model": ["a", "a", "b", "b"],
            "Sigma_SFR": 10.0 ** np.array([0.0, 2.0, 2.0, 4.0]),
        }
    )
    row = diagnostics.scatter_budget(
        reference, temporal, outcomes=("Sigma_SFR",)
    ).iloc[0]
    assert row["sigma_phangs_dex"] == pytest.approx(root_two)
    assert row["sigma_environment_dex"] == pytest.approx(1.0)
    assert row["sigma_temporal_dex"] == pytest.approx(1.0)
    assert row["sigma_suite_total_dex"] == pytest.approx(root_two)
    assert row["suite_to_phangs_variance"] == pytest.approx(1.0)
    assert row["environment_variance_fraction_of_suite"] == pytest.approx(0.5)
    assert row["temporal_variance_fraction_of_suite"] == pytest.approx(0.5)
