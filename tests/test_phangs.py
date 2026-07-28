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
    # - {name: beta_CO21_URC, datatype: float64}
    # - {name: e_beta_CO21_URC, datatype: float64}
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
    r_gal V_circ_CO21_URC beta_CO21_URC e_beta_CO21_URC Sigma_mol Sigma_atom Sigma_star rho_star_mp Sigma_SFR_HaW4recal Zprime
    1.0 150.0 0.0 0.05 10.0 5.0 200.0 0.1 0.01 1.0
    2.0 180.0 0.1 0.05 8.0 4.0 150.0 0.08 0.008 0.9
    3.0 200.0 0.2 0.05 5.0 3.0 100.0 0.05 0.005 0.8
    4.0 210.0 0.3 0.05 3.0 2.0 50.0 0.02 0.002 0.7
    5.0 220.0 0.4 0.05 1.0 1.5 20.0 0.008 0.001 0.6
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
        for col in (
            "r_gal",
            "Sigma_mol",
            "Sigma_atom",
            "Sigma_star",
            "rho_star_mp",
            "V_circ_CO21_URC",
            "Sigma_SFR_HaW4recal",
        ):
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
        for col in ("Sigma_gas", "Omega", "Omega_d", "qshear", "H_star"):
            assert col in out.colnames, f"Missing derived column: {col}"

    def test_sigma_gas_equals_mol_plus_atom(self, synthetic_table):
        out = phangs.compute_prfm_inputs(synthetic_table)
        expected = synthetic_table["Sigma_mol"] + synthetic_table["Sigma_atom"]
        np.testing.assert_allclose(out["Sigma_gas"].value, expected.value)

    def test_compute_inputs_accepts_aperture_specific_gas_and_sfr_columns(
        self, synthetic_table
    ):
        t = synthetic_table.copy()
        t["Sigma_mol_gauss"] = 2.0 * t["Sigma_mol"]
        t["Sigma_atom_gauss"] = 3.0 * t["Sigma_atom"]
        t["e_Sigma_mol_gauss"] = 0.1 * t["Sigma_mol_gauss"]
        t["e_Sigma_atom_gauss"] = 0.1 * t["Sigma_atom_gauss"]
        t["Sigma_SFR_HaW4recal_gauss"] = 4.0 * t["Sigma_SFR_HaW4recal"]

        out = phangs.compute_prfm_inputs(
            t,
            sigma_mol_col="Sigma_mol_gauss",
            sigma_atom_col="Sigma_atom_gauss",
            sfr_suffix="_gauss",
        )

        expected = t["Sigma_mol_gauss"] + t["Sigma_atom_gauss"]
        np.testing.assert_allclose(out["Sigma_gas"].value, expected.value)
        np.testing.assert_allclose(out["Sigma_mol"].value, t["Sigma_mol_gauss"].value)
        np.testing.assert_allclose(out["Sigma_atom"].value, t["Sigma_atom_gauss"].value)
        np.testing.assert_allclose(
            out["Sigma_SFR_HaW4recal"].value,
            t["Sigma_SFR_HaW4recal_gauss"].value,
        )

    def test_omega_formula(self, synthetic_table):
        """Omega = V_circ / r_gal  (km/s/kpc)."""
        out = phangs.compute_prfm_inputs(synthetic_table)
        V = synthetic_table["V_circ_CO21_URC"].value  # km/s
        r = synthetic_table["r_gal"].value  # kpc
        expected = V / r
        np.testing.assert_allclose(out["Omega"].value, expected, rtol=1e-10)
        np.testing.assert_allclose(out["Omega_d"].value, expected, rtol=1e-10)

    def test_qshear_formula(self, synthetic_table):
        """qshear = 1 - beta_CO21_URC."""
        out = phangs.compute_prfm_inputs(synthetic_table)
        expected = 1.0 - synthetic_table["beta_CO21_URC"]
        np.testing.assert_allclose(out["qshear"], expected, rtol=1e-10)

    def test_h_star_formula(self, synthetic_table):
        """H_star = Sigma_star / (4 * rho_star_mp)  [pc].

        The factor of 4 corresponds to H_star being the sech^2 scale height
        (rho ~ sech^2(z / (2 H_star)), so Sigma_star = 4 rho_star_mp H_star),
        matching the TIGRESS `zstar` convention that H_star maps onto.
        """
        out = phangs.compute_prfm_inputs(synthetic_table)
        Ss = synthetic_table["Sigma_star"].value  # M_sun/pc^2
        rs = synthetic_table["rho_star_mp"].value  # M_sun/pc^3
        expected = Ss / (4.0 * rs)  # pc
        np.testing.assert_allclose(out["H_star"].value, expected, rtol=1e-10)

    def test_units_attached(self, synthetic_table):
        out = phangs.compute_prfm_inputs(synthetic_table)
        assert out["Sigma_gas"].unit is not None
        assert out["Omega"].unit is not None
        assert out["Omega_d"].unit is not None
        assert out["H_star"].unit is not None

    def test_sigma_gas_positive(self, synthetic_table):
        out = phangs.compute_prfm_inputs(synthetic_table)
        assert np.all(out["Sigma_gas"].value > 0)

    def test_h_star_positive(self, synthetic_table):
        out = phangs.compute_prfm_inputs(synthetic_table)
        assert np.all(out["H_star"].value > 0)

    def test_omega_positive(self, synthetic_table):
        out = phangs.compute_prfm_inputs(synthetic_table)
        assert np.all(out["Omega"].value > 0)


# ---------------------------------------------------------------------------
# Filtering helpers
# ---------------------------------------------------------------------------


class TestFiltering:
    def test_valid_rows_removes_nans(self):
        t = Table(
            {
                "Sigma_gas": [1.0, float("nan"), 3.0],
                "Omega": [1.0, 2.0, float("nan")],
                "H_star": [100.0, 200.0, 300.0],
            }
        )
        mask = phangs.valid_rows(t, cols=["Sigma_gas", "Omega", "H_star"])
        assert mask.sum() == 1
        assert mask[0] is np.bool_(True)

    def test_valid_rows_removes_non_positive(self):
        t = Table(
            {
                "Sigma_gas": [1.0, 0.0, -1.0, 5.0],
                "Omega": [1.0, 2.0, 3.0, 4.0],
            }
        )
        mask = phangs.valid_rows(t, cols=["Sigma_gas", "Omega"])
        assert mask.sum() == 2  # rows 0 and 3

    def test_valid_rows_all_good(self, synthetic_table):
        out = phangs.compute_prfm_inputs(synthetic_table)
        mask = phangs.valid_rows(out, cols=["Sigma_gas", "Omega", "H_star"])
        assert mask.sum() == len(out)


class TestConfiguredLoading:
    def test_load_configured_phangs_joins_apertures(self, tmp_path):
        data_dir = tmp_path / "phangs_megatable"
        data_dir.mkdir()
        hex_file = data_dir / "TESTGAL_hexagon_1p5kpc.ecsv"
        hex_file.write_text(
            textwrap.dedent("""\
                # %ECSV 1.0
                # ---
                # datatype:
                # - {name: ID, datatype: int64}
                # - {name: r_gal, unit: kpc, datatype: float64}
                # - {name: RA, unit: deg, datatype: float64}
                # - {name: DEC, unit: deg, datatype: float64}
                # - {name: V_circ_CO21_URC, unit: km / s, datatype: float64}
                # - {name: beta_CO21_URC, datatype: float64}
                # - {name: Sigma_mol, unit: solMass / pc2, datatype: float64}
                # - {name: Sigma_atom, unit: solMass / pc2, datatype: float64}
                # - {name: Sigma_star, unit: solMass / pc2, datatype: float64}
                # - {name: rho_star_mp, unit: solMass / pc3, datatype: float64}
                # meta:
                #   GALAXY: TESTGAL
                # schema: astropy-2.0
                ID r_gal RA DEC V_circ_CO21_URC beta_CO21_URC Sigma_mol Sigma_atom Sigma_star rho_star_mp
                1 1.0 10.0 20.0 150.0 0.2 10.0 5.0 200.0 0.1
            """)
        )
        gauss_file = data_dir / "TESTGAL_gauss_1p5kpc.ecsv"
        gauss_file.write_text(
            textwrap.dedent("""\
                # %ECSV 1.0
                # ---
                # datatype:
                # - {name: ID, datatype: int64}
                # - {name: r_gal, unit: kpc, datatype: float64}
                # - {name: Sigma_mol_gauss, unit: solMass / pc2, datatype: float64}
                # - {name: e_Sigma_mol_gauss, unit: solMass / pc2, datatype: float64}
                # - {name: Sigma_atom_gauss, unit: solMass / pc2, datatype: float64}
                # - {name: e_Sigma_atom_gauss, unit: solMass / pc2, datatype: float64}
                # - {name: Sigma_SFR_HaW4recal_gauss, unit: solMass / (kpc2 yr), datatype: float64}
                # meta:
                #   GALAXY: TESTGAL
                # schema: astropy-2.0
                ID r_gal Sigma_mol_gauss e_Sigma_mol_gauss Sigma_atom_gauss e_Sigma_atom_gauss Sigma_SFR_HaW4recal_gauss
                1 1.0 20.0 2.0 7.0 0.7 0.03
            """)
        )
        config_file = tmp_path / "phangs.yml"
        config_file.write_text(
            textwrap.dedent("""\
                data_dir: phangs_megatable
                apertures:
                  context: hexagon
                  canonical: gauss
                join:
                  keys: [GALAXY, ID]
                  join_type: inner
                  table_names: [hexagon, gauss]
                  uniq_col_name: "{col_name}_{table_name}"
                canonical:
                  gas:
                    sigma_mol_col: Sigma_mol_gauss
                    sigma_atom_col: Sigma_atom_gauss
                  sfr_suffix: _gauss
                preserve_context_fields: [Sigma_mol, Sigma_atom]
                geometry_fields: [RA, DEC, r_gal]
                plot_aperture: gauss
            """)
        )

        loaded = phangs.load_configured_phangs(config_file, base_dir=tmp_path)
        out = loaded["table"]

        assert loaded["apertures"] == ("hexagon", "gauss")
        assert loaded["plot_aperture"] == "gauss"
        assert len(out) == 1
        np.testing.assert_allclose(out["Sigma_gas"].value, [27.0])
        np.testing.assert_allclose(out["Sigma_mol"].value, [20.0])
        np.testing.assert_allclose(out["Sigma_atom"].value, [7.0])
        np.testing.assert_allclose(out["Sigma_mol_hexagon"].value, [10.0])
        np.testing.assert_allclose(out["Sigma_SFR_HaW4recal"].value, [0.03])
        np.testing.assert_allclose(out["Omega"].value, [150.0])
        np.testing.assert_allclose(out["qshear"], [0.8])


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
            # - {name: beta_CO21_URC, datatype: float64}
            # - {name: Sigma_mol, unit: solMass / pc2, datatype: float64}
            # - {name: Sigma_atom, unit: solMass / pc2, datatype: float64}
            # - {name: Sigma_star, unit: solMass / pc2, datatype: float64}
            # - {name: rho_star_mp, unit: solMass / pc3, datatype: float64}
            # - {name: Sigma_SFR_HaW4recal, unit: solMass / (kpc2 yr), datatype: float64}
            # - {name: Zprime, datatype: float64}
            # meta:
            #   GALAXY: TESTNAN
            # schema: astropy-2.0
            r_gal V_circ_CO21_URC beta_CO21_URC Sigma_mol Sigma_atom Sigma_star rho_star_mp Sigma_SFR_HaW4recal Zprime
            1.0 150.0 0.0 10.0 5.0 200.0 0.1 0.01 1.0
            2.0 180.0 0.1 nan 4.0 150.0 0.08 0.008 0.9
            3.0 200.0 0.2 nan nan 100.0 0.05 0.005 0.8
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

    def test_variation_omega_none(self, synthetic_table):
        """variation={'Omega': None} disables rotation term; output still finite."""
        out = phangs.run_prfm(synthetic_table, variation={"Omega": None})
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

    def test_prfm_cols_without_omega(self, synthetic_table):
        """prfm_cols without Omega (gravity-only) produces finite outputs."""
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


class TestSamplingConfig:
    def test_default_synthesis_method_is_kde_sobol(self):
        from prfm.phangs_sampling import SamplingConfig
        cfg = SamplingConfig()
        assert cfg.synthesis_method == "kde_sobol"

    def test_design_fields_are_five(self):
        from prfm.phangs_sampling import SamplingConfig
        cfg = SamplingConfig()
        assert cfg.design_fields == [
            "Sigma_gas", "Sigma_star", "H_star", "Omega", "qshear"
        ]

    def test_qshear_max_default(self):
        from prfm.phangs_sampling import SamplingConfig
        assert SamplingConfig().qshear_max == 1.5

    def test_sobol_seed_default(self):
        from prfm.phangs_sampling import SamplingConfig
        assert SamplingConfig().sobol_seed == 42

    def test_kde_aux_sample_size_default(self):
        from prfm.phangs_sampling import SamplingConfig
        assert SamplingConfig().kde_aux_sample_size == 100_000


class TestMarginalQuantiles:
    """Unit test for _compute_marginal_quantiles — uses synthetic DataFrame."""

    def _make_designer(self):
        from prfm.phangs_sampling import PHANGSSamplingDesigner, SamplingConfig
        from astropy.table import Table
        import numpy as np

        rng = np.random.default_rng(0)
        n = 200
        t = Table({
            "Sigma_gas":  rng.lognormal(np.log(10.0), 0.3, n),
            "Sigma_star": rng.lognormal(np.log(50.0), 0.5, n),
            "H_star":     rng.lognormal(np.log(300.0), 0.4, n),
            "Omega":      rng.lognormal(np.log(30.0), 0.4, n),
            "qshear":     rng.uniform(0.5, 1.4, n),
        })
        cfg = SamplingConfig(kde_aux_sample_size=5_000)
        return PHANGSSamplingDesigner(t, config=cfg)

    def test_returns_dict_with_five_callables(self):
        d = self._make_designer()
        ref = d.table
        fields = d.config.design_fields
        qfns = d._compute_marginal_quantiles(ref, fields)
        assert set(qfns.keys()) == set(fields)
        for fn in qfns.values():
            assert callable(fn)

    def test_quantile_fn_maps_zero_to_min(self):
        import numpy as np
        d = self._make_designer()
        qfns = d._compute_marginal_quantiles(d.table, d.config.design_fields)
        q0 = qfns["Sigma_gas"](np.array([0.01]))
        q1 = qfns["Sigma_gas"](np.array([0.99]))
        assert q0 < q1

    def test_quantile_fn_output_in_log_space(self):
        """Q_j maps [0,1] -> log10 values; log10(Sigma_gas~10) should be near 1."""
        import numpy as np
        d = self._make_designer()
        qfns = d._compute_marginal_quantiles(d.table, d.config.design_fields)
        median_log = qfns["Sigma_gas"](np.array([0.5]))[0]
        assert 0.5 < median_log < 1.5  # log10(10) = 1


class TestKdeSobol:
    """Unit tests for synthesize_kde_sobol — no real PHANGS data needed."""

    @pytest.fixture
    def designer_and_reference(self):
        from prfm.phangs_sampling import PHANGSSamplingDesigner, SamplingConfig
        from astropy.table import Table
        import numpy as np

        rng = np.random.default_rng(7)
        n = 500
        qshear_vals = rng.uniform(0.3, 1.4, n)
        t = Table({
            "Sigma_gas":  rng.lognormal(np.log(10.0), 0.3, n),
            "Sigma_star": rng.lognormal(np.log(50.0), 0.5, n),
            "H_star":     rng.lognormal(np.log(300.0), 0.4, n),
            "Omega":      rng.lognormal(np.log(30.0), 0.4, n),
            "qshear":     qshear_vals,
            # validation fields — not used in KDE fit
            "Sigma_mol":  rng.lognormal(np.log(5.0), 0.5, n),
            "Sigma_atom": rng.lognormal(np.log(5.0), 0.5, n),
        })
        cfg = SamplingConfig(kde_aux_sample_size=5_000, sobol_seed=42)
        return PHANGSSamplingDesigner(t, config=cfg), t

    def test_returns_dataframe_with_correct_shape(self, designer_and_reference):
        import pandas as pd
        d, ref = designer_and_reference
        result = d.synthesize_kde_sobol(ref, n_samples=64)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 64

    def test_design_fields_present(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.synthesize_kde_sobol(ref, n_samples=64)
        for f in ["Sigma_gas", "Sigma_star", "H_star", "Omega", "qshear"]:
            assert f in result.columns

    def test_all_values_positive(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.synthesize_kde_sobol(ref, n_samples=64)
        for f in d.config.design_fields:
            assert (result[f] > 0).all(), f"{f} has non-positive values"

    def test_qshear_bounded(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.synthesize_kde_sobol(ref, n_samples=64)
        assert (result["qshear"] <= 1.5).all()

    def test_nesting_property_64_subset_of_128(self, designer_and_reference):
        """First 64 rows of n=128 must equal the n=64 design."""
        import pandas as pd
        d, ref = designer_and_reference
        s64 = d.synthesize_kde_sobol(ref, n_samples=64)
        s128 = d.synthesize_kde_sobol(ref, n_samples=128)
        pd.testing.assert_frame_equal(
            s64.reset_index(drop=True),
            s128.iloc[:64].reset_index(drop=True),
            check_like=False,
        )

    def test_validation_fields_not_in_output(self, designer_and_reference):
        """f_mol and SFR should not appear in the Sobol sample output."""
        d, ref = designer_and_reference
        result = d.synthesize_kde_sobol(ref, n_samples=64)
        for f in ["Sigma_mol", "Sigma_atom", "Sigma_SFR_HaW4recal"]:
            assert f not in result.columns

    def test_attrs_stores_n_extra(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.synthesize_kde_sobol(ref, n_samples=64)
        assert "n_extra" in result.attrs
        assert result.attrs["n_extra"] >= 0


class TestExpandedKdeSobol:
    @pytest.fixture
    def designer_and_reference(self):
        from prfm.phangs_sampling import PHANGSSamplingDesigner, SamplingConfig
        from astropy.table import Table
        import numpy as np

        rng = np.random.default_rng(99)
        n = 400
        t = Table({
            "Sigma_gas":  rng.lognormal(np.log(10.0), 0.3, n),
            "Sigma_star": rng.lognormal(np.log(50.0), 0.5, n),
            "H_star":     rng.lognormal(np.log(300.0), 0.4, n),
            "Omega":      rng.lognormal(np.log(30.0), 0.4, n),
            "qshear":     rng.uniform(0.4, 1.3, n),
        })
        cfg = SamplingConfig(kde_aux_sample_size=5_000)
        return PHANGSSamplingDesigner(t, config=cfg), t

    def test_expanded_sobol_has_correct_shape(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.synthesize_expanded_kde_sobol(ref, n_samples=32)
        assert len(result) == 32
        assert set(d.config.design_fields).issubset(result.columns)

    def test_sample_dispatch_kde_sobol(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.sample(ref, n_samples=32, method="kde_sobol")
        assert len(result) == 32

    def test_sample_dispatch_expanded(self, designer_and_reference):
        d, ref = designer_and_reference
        result = d.sample(ref, n_samples=32, method="expanded_kde_sobol")
        assert len(result) == 32

    def test_sample_uses_default_method(self, designer_and_reference):
        """Default method is now kde_sobol."""
        import pandas as pd
        d, ref = designer_and_reference
        result = d.sample(ref, n_samples=32)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 32


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
