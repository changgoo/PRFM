"""Unit tests for the pure helpers of project/scripts/augment_design.py.

The data-dependent generation (select_reference / generate_block) requires the
PHANGS megatable and is exercised end-to-end elsewhere; here we cover the
argument parsing and block-assembly logic that need no data.
"""

import importlib.util
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]


def _load(name: str):
    spec = importlib.util.spec_from_file_location(
        name, _ROOT / "project" / "scripts" / f"{name}.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ad = _load("augment_design")


def test_parse_augment():
    assert ad.parse_augment("H_star=0.3,Omega=0.3") == {"H_star": 0.3, "Omega": 0.3}
    assert ad.parse_augment("H_star=0.25") == {"H_star": 0.25}
    assert ad.parse_augment("") == {}


def test_parse_sigma_gas_from_stem():
    assert ad.parse_sigma_gas_from_stem("design_Sgas10.0_n0032") == 10.0
    assert ad.parse_sigma_gas_from_stem("design_Sgas5.0_n0064") == 5.0
    assert ad.parse_sigma_gas_from_stem("nothing_here") is None


def test_add_block_suffixes_continues_from_offset():
    block = pd.DataFrame({"Sigma_gas": [1.0, 2.0, 3.0]})
    out = ad.add_block_suffixes(block, offset=32)
    assert list(out["suffix"]) == ["row0032", "row0033", "row0034"]
    # original columns preserved
    assert list(out["Sigma_gas"]) == [1.0, 2.0, 3.0]


def test_default_augment_is_hstar_omega():
    assert ad.DEFAULT_AUGMENT == {"H_star": 0.3, "Omega": 0.3}


def test_block_suffixes_do_not_collide_with_existing():
    existing_n = 32
    block = pd.DataFrame({"Sigma_gas": [1.0] * 4})
    out = ad.add_block_suffixes(block, offset=existing_n)
    existing_suffixes = {f"row{i:04d}" for i in range(existing_n)}
    assert existing_suffixes.isdisjoint(set(out["suffix"]))
