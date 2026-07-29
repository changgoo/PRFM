"""Unit tests for project/scripts/csv_to_slurm_yaml.py.

Network- and data-free: builds a tiny design CSV in a temp dir and checks the
suite-YAML config it produces (pdf2d bin derivation, restart_overrides, Omega
unit conversion, per-row model structure).
"""

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_SCRIPT = _ROOT / "project" / "scripts" / "csv_to_slurm_yaml.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("csv_to_slurm_yaml", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


csy = _load_module()


def test_pdf2d_overrides_derived_from_box():
    """Bins = box/dx along the projected axes; ranges span +/- box/2."""
    domain = {"dx": 8, "box_size": [1024, 1024, 4096], "grid_size": [32, 32, 32]}
    ov = csy.build_pdf2d_overrides(domain)

    # output8 = x-y face-on
    assert ov["output8/Nbinx"] == 128
    assert ov["output8/Nbiny"] == 128
    assert (ov["output8/binx_min"], ov["output8/binx_max"]) == (-512, 512)
    assert (ov["output8/biny_min"], ov["output8/biny_max"]) == (-512, 512)

    # output9 = x-z edge-on (vertical axis spans the full 4096 pc box)
    assert ov["output9/Nbinx"] == 128
    assert ov["output9/Nbiny"] == 512
    assert (ov["output9/binx_min"], ov["output9/binx_max"]) == (-512, 512)
    assert (ov["output9/biny_min"], ov["output9/biny_max"]) == (-2048, 2048)


def test_pdf2d_bins_track_resolution():
    """Halving dx doubles the bin counts; ranges are unchanged."""
    domain = {"dx": 4, "box_size": [1024, 1024, 4096], "grid_size": [64, 64, 64]}
    ov = csy.build_pdf2d_overrides(domain)
    assert ov["output8/Nbinx"] == 256
    assert ov["output9/Nbiny"] == 1024
    assert ov["output9/biny_max"] == 2048  # range still spans the box


def _write_csv(tmp_path: Path) -> Path:
    df = pd.DataFrame(
        {
            "Sigma_gas": [11.0192, 14.2106],
            "Sigma_star": [148.4237, 33.2941],
            "H_star": [556.4261, 199.6901],
            "Omega": [17.943, 31.6103],  # km/s/kpc (PHANGS convention)
            "qshear": [0.8072, 1.0916],
        }
    )
    csv = tmp_path / "design_test.csv"
    df.to_csv(csv, index=False)
    return csv


def test_build_config_blocks_present(tmp_path):
    csv = _write_csv(tmp_path)
    cfg = csy.build_config(csv, base="R8_8pc", decimals=4, machine_yaml=None)

    assert set(cfg) == {"domain", "fixed_overrides", "restart_overrides", "models"}

    fixed = cfg["fixed_overrides"]
    # pdf2d keys come first, then the scalar physics/output defaults.
    assert next(iter(fixed)) == "output8/Nbinx"
    assert fixed["problem/rhodm"] == 0.0064
    assert fixed["problem/beta"] == 10

    # restart_overrides carries exactly the derived pdf2d geometry.
    assert cfg["restart_overrides"] == csy.build_pdf2d_overrides(cfg["domain"])


def test_build_config_models_and_omega(tmp_path):
    csv = _write_csv(tmp_path)
    cfg = csy.build_config(csv, base="R8_8pc", decimals=4, machine_yaml=None)

    models = cfg["models"]
    assert [m["suffix"] for m in models] == ["row0000", "row0001"]

    m0 = models[0]["extra_overrides"]
    assert m0["problem/surf"] == 11.0192
    assert m0["problem/SurfS"] == 148.4237
    assert m0["problem/zstar"] == 556.4261
    assert m0["problem/qshear"] == 0.8072
    # Omega converted km/s/kpc -> km/s/pc (divide by 1000).
    assert m0["problem/Omega"] == pytest.approx(0.017943)


def test_unknown_base_omits_domain_and_pdf2d(tmp_path):
    """A base without domain defaults yields no domain/restart/pdf2d blocks."""
    csv = _write_csv(tmp_path)
    cfg = csy.build_config(csv, base="unknown_base", decimals=4, machine_yaml=None)
    assert "domain" not in cfg
    assert "restart_overrides" not in cfg
    assert not any(k.startswith("output8/") for k in cfg["fixed_overrides"])


def test_five_column_csv_keeps_fixed_physics(tmp_path):
    """Backward compat: a physical-only CSV keeps Z_gas/Z_dust fixed at 1.0."""
    csv = _write_csv(tmp_path)
    cfg = csy.build_config(csv, base="R8_8pc", decimals=4, machine_yaml=None)
    fixed = cfg["fixed_overrides"]
    assert fixed["problem/Z_gas"] == 1.0
    assert fixed["problem/Z_dust"] == 1.0
    # xi_CR_amp is never emitted unless it varies (athinput default applies).
    assert "problem/xi_CR_amp" not in fixed
    for m in cfg["models"]:
        assert "problem/Z_gas" not in m["extra_overrides"]


def _write_physics_csv(tmp_path: Path) -> Path:
    df = pd.DataFrame(
        {
            "Sigma_gas": [11.0192, 11.0192],
            "Sigma_star": [148.4237, 148.4237],
            "H_star": [556.4261, 556.4261],
            "Omega": [17.943, 17.943],
            "qshear": [0.8072, 0.8072],
            "Z_gas": [1.0, 0.1],
            "Z_dust": [1.0, 0.046],  # f_dtm * Z_gas
            "xi_CR_amp": [1.0, 0.22],
            "suffix": ["a07_base", "a07_Z0p10"],
        }
    )
    csv = tmp_path / "ablation_test.csv"
    df.to_csv(csv, index=False)
    return csv


def test_optional_physics_columns_vary_per_row(tmp_path):
    csv = _write_physics_csv(tmp_path)
    cfg = csy.build_config(csv, base="R8_8pc", decimals=4, machine_yaml=None)

    fixed = cfg["fixed_overrides"]
    # Now that they vary, the fixed defaults must be gone (no conflict).
    assert "problem/Z_gas" not in fixed
    assert "problem/Z_dust" not in fixed
    assert "problem/xi_CR_amp" not in fixed
    # Non-physics fixed params remain.
    assert fixed["problem/rhodm"] == 0.0064

    m1 = cfg["models"][1]["extra_overrides"]
    assert m1["problem/Z_gas"] == 0.1
    assert m1["problem/Z_dust"] == 0.046
    assert m1["problem/xi_CR_amp"] == 0.22


def test_suffix_column_overrides_row_index(tmp_path):
    csv = _write_physics_csv(tmp_path)
    cfg = csy.build_config(csv, base="R8_8pc", decimals=4, machine_yaml=None)
    assert [m["suffix"] for m in cfg["models"]] == ["a07_base", "a07_Z0p10"]
