"""Smoke test for project/scripts/plot_augmentation_coverage.py.

Verifies the figure builder returns a well-formed 2x3 marginals figure from
in-memory DataFrames (no PHANGS data, no display).
"""

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

_ROOT = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location(
        "plot_augmentation_coverage",
        _ROOT / "project" / "scripts" / "plot_augmentation_coverage.py",
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pac = _load()


def _design(n, seed):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "Sigma_gas":  10 ** rng.normal(1.0, 0.1, n),
        "Sigma_star": 10 ** rng.normal(1.7, 0.3, n),
        "H_star":     10 ** rng.normal(2.5, 0.2, n),
        "Omega":      10 ** rng.normal(1.4, 0.2, n),
        "qshear":     10 ** rng.normal(-0.05, 0.1, n),
    })


def test_build_coverage_figure_shape():
    fig = pac.build_coverage_figure(_design(300, 0), _design(32, 1), _design(32, 2))
    assert isinstance(fig, Figure)
    # 2x3 grid: five field panels + one legend panel.
    assert len(fig.axes) == 6


def test_augmented_fields_labeled():
    fig = pac.build_coverage_figure(
        _design(300, 0), _design(32, 1), _design(32, 2),
        augmented=["H_star", "Omega"],
    )
    xlabels = [ax.get_xlabel() for ax in fig.axes]
    augmented_labeled = [x for x in xlabels if "augmented" in x]
    assert len(augmented_labeled) == 2
