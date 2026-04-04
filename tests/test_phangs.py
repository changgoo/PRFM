"""
TDD tests for prfm.phangs — PHANGS megatable data handling.

Tests are split into two groups:
  - Unit tests (no network/disk required): use minimal synthetic ECSV data
    to verify parsing, derived-quantity formulas, and filtering logic.
  - Integration tests (require downloaded data): marked with
    @pytest.mark.integration and skipped automatically when
    data/phangs_megatable/ is absent.

Run only unit tests (fast):
    pytest tests/test_phangs.py -m "not integration"

Run all including integration:
    pytest tests/test_phangs.py
"""

import textwrap
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

DATA_DIR = Path(__file__).parent.parent / "data" / "phangs_megatable"
HAS_DATA = DATA_DIR.exists() and any(DATA_DIR.glob("*.ecsv"))

integration = pytest.mark.skipif(not HAS_DATA, reason="PHANGS data not downloaded")

# ---------------------------------------------------------------------------
# Minimal synthetic ECSV fixture
# ---------------------------------------------------------------------------

SYNTHETIC_ECSV = textwrap.dedent("""\
    # %ECSV 1.0
    # ---
    # datatype:
    # - {name: r_gal, unit: kpc, datatype: float64}
    # - {name: V_circ_CO21_URC, unit: km / s, datatype: float64}
    # - {name: Sigma_mol, unit: solMass / pc2, datatype: float64}
    # - {name: Sigma_atom, unit: solMass / pc2, datatype: float64}
    # - {name: Sigma_star, unit: solMass / pc2, datatype: float64}
    # - {name: rho_star_mp, unit: solMass / pc3, datatype: float64}
    # - {name: Sigma_SFR_HaW4recal, unit: solMass / (kpc2 yr), datatype: float64}
    # - {name: Zprime, datatype: float64}
    # meta:
    #   GALAXY: TESTGAL
    #   DIST_MPC: 10.0
    #   R25_KPC: 5.0
    # schema: astropy-2.0
    r_gal V_circ_CO21_URC Sigma_mol Sigma_atom Sigma_star rho_star_mp Sigma_SFR_HaW4recal Zprime
    1.0 150.0 10.0 5.0 200.0 0.1 0.01 1.0
    2.0 180.0 8.0 4.0 150.0 0.08 0.008 0.9
    3.0 200.0 5.0 3.0 100.0 0.05 0.005 0.8
    4.0 210.0 3.0 2.0 50.0 0.02 0.002 0.7
    5.0 220.0 1.0 1.5 20.0 0.008 0.001 0.6
""")


@pytest.fixture
def synthetic_table(synthetic_file):
    return phangs.load(synthetic_file)


@pytest.fixture
def synthetic_file(tmp_path):
    f = tmp_path / "TESTGAL_annulus_0p5kpc.ecsv"
    f.write_text(SYNTHETIC_ECSV)
    return f


# ---------------------------------------------------------------------------
# Import the module under test
# ---------------------------------------------------------------------------

from prfm import phangs  # noqa: E402


# ---------------------------------------------------------------------------
# File listing and URL construction
# ---------------------------------------------------------------------------


class TestFileListAndURL:
    def test_list_files_returns_list(self):
        files = phangs.list_files()
        assert isinstance(files, list)
        assert len(files) > 0

    def test_list_files_ecsv_extension(self):
        for f in phangs.list_files():
            assert f.endswith(".ecsv")

    def test_list_files_filter_aperture(self):
        annulus = phangs.list_files(aperture="annulus")
        assert all("annulus" in f for f in annulus)
        hexagon = phangs.list_files(aperture="hexagon")
        assert all("hexagon" in f for f in hexagon)

    def test_list_files_all_three_apertures(self):
        for ap in ("annulus", "gauss", "hexagon"):
            files = phangs.list_files(aperture=ap)
            assert len(files) > 0, f"No files found for aperture={ap}"

    def test_build_url_contains_filename(self):
        url = phangs.build_url("NGC0628_annulus_0p5kpc.ecsv")
        assert "NGC0628_annulus_0p5kpc.ecsv" in url
        assert url.startswith("http")

    def test_galaxies_list(self):
        galaxies = phangs.list_galaxies()
        assert "NGC0628" in galaxies
        assert "IC1954" in galaxies
        assert len(galaxies) >= 10

    def test_filename_for_galaxy(self):
        fname = phangs.filename("NGC0628", aperture="annulus")
        assert fname == "NGC0628_annulus_0p5kpc.ecsv"
        fname = phangs.filename("NGC0628", aperture="hexagon")
        assert fname == "NGC0628_hexagon_1p5kpc.ecsv"


# ---------------------------------------------------------------------------
# Loading a single table
# ---------------------------------------------------------------------------


class TestLoadTable:
    def test_load_returns_astropy_table(self, synthetic_file):
        t = phangs.load(synthetic_file)
        assert isinstance(t, Table)

    def test_load_has_required_columns(self, synthetic_file):
        t = phangs.load(synthetic_file)
        for col in ("r_gal", "Sigma_mol", "Sigma_atom", "Sigma_star",
                    "rho_star_mp", "V_circ_CO21_URC", "Sigma_SFR_HaW4recal"):
            assert col in t.colnames, f"Missing column: {col}"

    def test_load_preserves_units(self, synthetic_file):
        t = phangs.load(synthetic_file)
        assert t["r_gal"].unit is not None
        assert t["Sigma_mol"].unit is not None

    def test_load_attaches_galaxy_name(self, synthetic_file):
        t = phangs.load(synthetic_file)
        assert t.meta.get("GALAXY") == "TESTGAL"

    def test_load_nonexistent_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            phangs.load(tmp_path / "does_not_exist.ecsv")


# ---------------------------------------------------------------------------
# Derived PRFM quantities
# ---------------------------------------------------------------------------


class TestComputePRFMInputs:
    def test_returns_table_with_new_columns(self, synthetic_table):
        out = phangs.compute_prfm_inputs(synthetic_table)
        for col in ("Sigma_gas", "Omega_d", "H_star"):
            assert col in out.colnames, f"Missing derived column: {col}"

    def test_sigma_gas_equals_mol_plus_atom(self, synthetic_table):
        out = phangs.compute_prfm_inputs(synthetic_table)
        expected = synthetic_table["Sigma_mol"] + synthetic_table["Sigma_atom"]
        np.testing.assert_allclose(out["Sigma_gas"].value, expected.value)

    def test_omega_d_formula(self, synthetic_table):
        """Omega_d = V_circ / r_gal  (km/s/kpc)"""
        out = phangs.compute_prfm_inputs(synthetic_table)
        V = synthetic_table["V_circ_CO21_URC"].value   # km/s
        r = synthetic_table["r_gal"].value              # kpc
        expected = V / r
        np.testing.assert_allclose(out["Omega_d"].value, expected, rtol=1e-10)

    def test_h_star_formula(self, synthetic_table):
        """H_star = Sigma_star / (2 * rho_star_mp)  [pc]"""
        out = phangs.compute_prfm_inputs(synthetic_table)
        Ss = synthetic_table["Sigma_star"].value   # M_sun/pc^2
        rs = synthetic_table["rho_star_mp"].value  # M_sun/pc^3
        expected = Ss / (2.0 * rs)                 # pc
        np.testing.assert_allclose(out["H_star"].value, expected, rtol=1e-10)

    def test_units_attached(self, synthetic_table):
        out = phangs.compute_prfm_inputs(synthetic_table)
        assert out["Sigma_gas"].unit is not None
        assert out["Omega_d"].unit is not None
        assert out["H_star"].unit is not None

    def test_sigma_gas_positive(self, synthetic_table):
        out = phangs.compute_prfm_inputs(synthetic_table)
        assert np.all(out["Sigma_gas"].value > 0)

    def test_h_star_positive(self, synthetic_table):
        out = phangs.compute_prfm_inputs(synthetic_table)
        assert np.all(out["H_star"].value > 0)

    def test_omega_d_positive(self, synthetic_table):
        out = phangs.compute_prfm_inputs(synthetic_table)
        assert np.all(out["Omega_d"].value > 0)


# ---------------------------------------------------------------------------
# Filtering helpers
# ---------------------------------------------------------------------------


class TestFiltering:
    def test_valid_rows_removes_nans(self):
        t = Table({
            "Sigma_gas": [1.0, float("nan"), 3.0],
            "Omega_d": [1.0, 2.0, float("nan")],
            "H_star": [100.0, 200.0, 300.0],
        })
        mask = phangs.valid_rows(t, cols=["Sigma_gas", "Omega_d", "H_star"])
        assert mask.sum() == 1
        assert mask[0] is np.bool_(True)

    def test_valid_rows_removes_non_positive(self):
        t = Table({
            "Sigma_gas": [1.0, 0.0, -1.0, 5.0],
            "Omega_d": [1.0, 2.0, 3.0, 4.0],
        })
        mask = phangs.valid_rows(t, cols=["Sigma_gas", "Omega_d"])
        assert mask.sum() == 2  # rows 0 and 3

    def test_valid_rows_all_good(self, synthetic_table):
        out = phangs.compute_prfm_inputs(synthetic_table)
        mask = phangs.valid_rows(out, cols=["Sigma_gas", "Omega_d", "H_star"])
        assert mask.sum() == len(out)


# ---------------------------------------------------------------------------
# run_prfm
# ---------------------------------------------------------------------------


class TestRunPRFM:
    def test_returns_table(self, synthetic_table):
        out = phangs.run_prfm(synthetic_table)
        assert isinstance(out, Table)

    def test_output_columns_present(self, synthetic_table):
        out = phangs.run_prfm(synthetic_table)
        for col in ("P_weight", "H_gas", "sigma_eff_sol", "Sigma_SFR_pred"):
            assert col in out.colnames, f"Missing column: {col}"

    def test_same_length_as_input(self, synthetic_table):
        out = phangs.run_prfm(synthetic_table)
        assert len(out) == len(synthetic_table)

    def test_output_units_attached(self, synthetic_table):
        out = phangs.run_prfm(synthetic_table)
        assert out["P_weight"].unit is not None
        assert out["H_gas"].unit is not None
        assert out["sigma_eff_sol"].unit is not None
        assert out["Sigma_SFR_pred"].unit is not None

    def test_all_rows_finite_for_clean_input(self, synthetic_table):
        out = phangs.run_prfm(synthetic_table)
        for col in ("P_weight", "H_gas", "sigma_eff_sol", "Sigma_SFR_pred"):
            assert np.all(np.isfinite(out[col].value)), f"Non-finite values in {col}"

    def test_invalid_rows_give_nan(self, synthetic_file, tmp_path):
        """Rows with NaN inputs should produce NaN outputs."""
        import textwrap
        ecsv_nan = textwrap.dedent("""\
            # %ECSV 1.0
            # ---
            # datatype:
            # - {name: r_gal, unit: kpc, datatype: float64}
            # - {name: V_circ_CO21_URC, unit: km / s, datatype: float64}
            # - {name: Sigma_mol, unit: solMass / pc2, datatype: float64}
            # - {name: Sigma_atom, unit: solMass / pc2, datatype: float64}
            # - {name: Sigma_star, unit: solMass / pc2, datatype: float64}
            # - {name: rho_star_mp, unit: solMass / pc3, datatype: float64}
            # - {name: Sigma_SFR_HaW4recal, unit: solMass / (kpc2 yr), datatype: float64}
            # - {name: Zprime, datatype: float64}
            # meta:
            #   GALAXY: TESTNAN
            # schema: astropy-2.0
            r_gal V_circ_CO21_URC Sigma_mol Sigma_atom Sigma_star rho_star_mp Sigma_SFR_HaW4recal Zprime
            1.0 150.0 10.0 5.0 200.0 0.1 0.01 1.0
            2.0 180.0 nan 4.0 150.0 0.08 0.008 0.9
            3.0 200.0 nan nan 100.0 0.05 0.005 0.8
        """)
        f = tmp_path / "TESTNAN_annulus_0p5kpc.ecsv"
        f.write_text(ecsv_nan)
        t = phangs.load(f)
        import warnings
        with warnings.catch_warnings(record=True):
            out = phangs.run_prfm(t)
        assert np.isfinite(out["Sigma_SFR_pred"].value[0])
        # mol=nan, atom=valid → Sigma_gas filled with atom only → valid output
        assert np.isfinite(out["Sigma_SFR_pred"].value[1])
        # both mol and atom nan → Sigma_gas=nan → output must be nan
        assert np.isnan(out["Sigma_SFR_pred"].value[2])

    def test_sfr_pred_positive_for_valid_rows(self, synthetic_table):
        out = phangs.run_prfm(synthetic_table)
        assert np.all(out["Sigma_SFR_pred"].value > 0)

    def test_p_weight_positive_for_valid_rows(self, synthetic_table):
        out = phangs.run_prfm(synthetic_table)
        assert np.all(out["P_weight"].value > 0)

    def test_h_gas_positive_for_valid_rows(self, synthetic_table):
        out = phangs.run_prfm(synthetic_table)
        assert np.all(out["H_gas"].value > 0)

    def test_accepts_precomputed_inputs(self, synthetic_table):
        """run_prfm should not recompute if Sigma_gas already present."""
        t = phangs.compute_prfm_inputs(synthetic_table)
        out = phangs.run_prfm(t)
        assert "P_weight" in out.colnames

    def test_zprime_none_skips_metallicity(self, synthetic_table):
        """zprime_col=None should still produce valid outputs."""
        out = phangs.run_prfm(synthetic_table, zprime_col=None)
        assert np.all(np.isfinite(out["Sigma_SFR_pred"].value))

    def test_input_table_not_mutated(self, synthetic_table):
        """run_prfm must return a copy; the original table is unchanged."""
        cols_before = list(synthetic_table.colnames)
        phangs.run_prfm(synthetic_table)
        assert list(synthetic_table.colnames) == cols_before

    def test_variation_omega_d_none(self, synthetic_table):
        """variation={'Omega_d': None} disables rotation term; output still finite."""
        out = phangs.run_prfm(synthetic_table, variation={"Omega_d": None})
        assert np.all(np.isfinite(out["Sigma_SFR_pred"].value))

    def test_variation_sigma_star_scaling(self, synthetic_table):
        """Scaling Sigma_star via variation changes P_weight relative to default."""
        out_default = phangs.run_prfm(synthetic_table)
        out_scaled = phangs.run_prfm(synthetic_table, variation={"Sigma_star": 2.0})
        assert np.all(out_scaled["P_weight"].value >= out_default["P_weight"].value)

    def test_variation_h_star_scaling(self, synthetic_table):
        """Scaling H_star via variation changes P_weight relative to default."""
        out_default = phangs.run_prfm(synthetic_table)
        out_scaled = phangs.run_prfm(synthetic_table, variation={"H_star": 2.0})
        assert not np.allclose(
            out_scaled["P_weight"].value, out_default["P_weight"].value
        )

    def test_prfm_cols_without_omega_d(self, synthetic_table):
        """prfm_cols without Omega_d (gravity-only) produces finite outputs."""
        out = phangs.run_prfm(
            synthetic_table,
            prfm_cols=["Sigma_gas", "Sigma_star", "H_star"],
        )
        assert np.all(np.isfinite(out["Sigma_SFR_pred"].value))

    def test_invalid_rows_emit_warning(self, synthetic_table):
        """A warning is emitted when any input row is invalid."""
        import warnings

        t = phangs.compute_prfm_inputs(synthetic_table.copy())
        t["Sigma_gas"][0] = np.nan
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            phangs.run_prfm(t)
        assert any("invalid" in str(w.message).lower() for w in caught)

    def test_alternative_models_produce_finite_output(self, synthetic_table):
        """Non-default sigma_eff and yield models should still produce finite outputs."""
        out = phangs.run_prfm(
            synthetic_table,
            sigma_eff_model="tigress-classic-mid",
            yield_model="tigress-classic",
        )
        assert np.all(np.isfinite(out["Sigma_SFR_pred"].value))


# ---------------------------------------------------------------------------
# Integration tests (require downloaded data)
# ---------------------------------------------------------------------------


@integration
class TestIntegrationLoad:
    def test_one_file_readable(self):
        files = sorted(DATA_DIR.glob("*_annulus_0p5kpc.ecsv"))
        assert len(files) > 0
        t = phangs.load(files[0])
        assert len(t) > 0

    def test_galaxy_name_in_meta(self):
        files = sorted(DATA_DIR.glob("*_annulus_0p5kpc.ecsv"))
        t = phangs.load(files[0])
        assert "GALAXY" in t.meta

    def test_prfm_inputs_computable(self):
        files = sorted(DATA_DIR.glob("*_annulus_0p5kpc.ecsv"))
        t = phangs.load(files[0])
        out = phangs.compute_prfm_inputs(t)
        assert "Sigma_gas" in out.colnames


@integration
class TestIntegrationLoadAll:
    def test_load_all_returns_stacked_table(self):
        t = phangs.load_all(DATA_DIR, aperture="annulus")
        assert isinstance(t, dict)
        t = phangs.vstack_tables(t)
        assert isinstance(t, Table)
        assert len(t) > 0
        assert "GALAXY" in t.colnames

    def test_load_all_contains_multiple_galaxies(self):
        t = phangs.load_all(DATA_DIR, aperture="annulus")
        t = phangs.vstack_tables(t)
        assert len(set(t["GALAXY"])) > 1
