"""Unit tests for project/scripts/make_ablation_grid.py.

Network- and data-free: builds a tiny physical design in memory and checks the
OAT ablation grid it produces (anchor spanning, one-at-a-time structure,
Z_dust = f_dtm * Z_gas, suffix format), plus that the output feeds cleanly into
csv_to_slurm_yaml.build_config.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]


def _load(name: str):
    spec = importlib.util.spec_from_file_location(
        name, _ROOT / "project" / "scripts" / f"{name}.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mag = _load("make_ablation_grid")
csy = _load("csv_to_slurm_yaml")


def _design(n: int = 6) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    return pd.DataFrame(
        {
            "Sigma_gas": np.full(n, 10.0),
            "Sigma_star": np.linspace(20.0, 200.0, n),  # monotone: easy to span
            "H_star": rng.uniform(200, 600, n),
            "Omega": rng.uniform(15, 40, n),
            "qshear": rng.uniform(0.5, 1.4, n),
        }
    )


def test_select_anchor_indices_spans_min_median_max():
    df = _design(6)  # Sigma_star increasing with index
    idx = mag.select_anchor_indices(df, 3, by="Sigma_star")
    assert idx == [0, 2, 5]  # min, median-ish, max positions (linspace + round)


def test_grid_shape_and_columns():
    df = _design(6)
    grid = mag.build_ablation_grid(df, anchor_indices=[0, 5])
    # 2 anchors * 3 params * 3 sub-fiducial levels (4 levels minus base 1.0)
    assert len(grid) == 2 * 3 * 3
    for col in ("Z_gas", "Z_dust", "xi_CR_amp", "suffix", "sweep_param", "anchor_row"):
        assert col in grid.columns
    for f in mag.PHYSICAL_FIELDS:
        assert f in grid.columns


def test_one_at_a_time_structure():
    df = _design(6)
    grid = mag.build_ablation_grid(df, anchor_indices=[0])
    for _, r in grid.iterrows():
        off_base = [r["Z_gas"] != 1.0, r["f_dtm"] != 1.0, r["xi_CR_amp"] != 1.0]
        assert sum(off_base) == 1  # exactly one knob perturbed
    # The (1,1,1) base corner is never emitted (those runs already exist).
    assert not (
        (grid["Z_gas"] == 1.0) & (grid["f_dtm"] == 1.0) & (grid["xi_CR_amp"] == 1.0)
    ).any()


def test_z_dust_is_fdtm_times_zgas():
    df = _design(6)
    grid = mag.build_ablation_grid(df, anchor_indices=[0, 3])
    np.testing.assert_allclose(grid["Z_dust"], grid["f_dtm"] * grid["Z_gas"])


def test_suffix_encodes_anchor_param_level():
    df = _design(6)
    grid = mag.build_ablation_grid(df, anchor_indices=[5])
    suffixes = set(grid["suffix"])
    assert "a0005_Z0p10" in suffixes  # Z_gas sweep, level 0.1
    assert "a0005_fd0p46" in suffixes  # f_dtm sweep, level 0.46
    assert "a0005_cr0p22" in suffixes  # xi_CR_amp sweep, level 0.22


def test_physical_params_match_anchor():
    df = _design(6)
    grid = mag.build_ablation_grid(df, anchor_indices=[3])
    for f in mag.PHYSICAL_FIELDS:
        np.testing.assert_allclose(grid[f].to_numpy(), df.at[3, f])


def test_grid_feeds_csv_to_slurm_yaml(tmp_path):
    df = _design(6)
    grid = mag.build_ablation_grid(df, anchor_indices=[0, 5])
    csv = tmp_path / "ablation.csv"
    grid.to_csv(csv, index=False)

    cfg = csy.build_config(csv, base="R8_8pc", decimals=4, machine_yaml=None)
    # Physics knobs vary per row -> dropped from the fixed block.
    assert "problem/Z_gas" not in cfg["fixed_overrides"]
    assert "problem/xi_CR_amp" not in cfg["fixed_overrides"]
    # Descriptive suffixes flow through from the CSV.
    assert all(m["suffix"].startswith("a") for m in cfg["models"])
    m0 = cfg["models"][0]["extra_overrides"]
    assert "problem/Z_gas" in m0 and "problem/xi_CR_amp" in m0
