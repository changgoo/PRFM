"""
Unit tests for prfm.prfm based on the toy galaxy model from book/prfm.ipynb.

The reference model is a radial profile (R=2-16 kpc) with:
  - Sigma_gas = 10 / (R/R0)  [M_sun/pc^2]  (linear in R)
  - Sigma_star = 42 * exp(R0/Re - R/Re)    (exponential disk)
  - Omega_d = 28 / (R/R0)   [km/s/kpc]     (flat rotation curve)
  - H_star = 245             [pc]            (no flaring)
  - sigma_eff = 15           [km/s]          (constant)
"""

import numpy as np
import pytest

import prfm
from prfm.prfm import (
    _Gconst_cgs,
    _kbol_cgs,
    _pc_cgs,
    _sfr_cgs,
    _surf_cgs,
    get_feedback_yield,
    get_feedback_yield_comp,
    get_pressure,
    get_sfr,
    get_sigma_eff,
    get_weight_dm,
    get_weight_gas,
    get_weight_star,
)

# ---------------------------------------------------------------------------
# Shared toy-galaxy fixture
# ---------------------------------------------------------------------------

_kms_kpc_cgs = 1.0e5 / (3.085677581e21)  # km/s/kpc -> 1/s


@pytest.fixture(scope="module")
def toy_model():
    """Toy galaxy radial profile matching the notebook example."""
    R0, Re = 8.0, 3.5
    R = np.linspace(2, 16, 100)
    r0 = R / R0
    re = R / Re

    Sigma_gas = 10.0 / r0
    Sigma_star = 42.0 * np.exp(R0 / Re - re)
    Omega_d = 28.0 / r0
    H_star = 245.0
    sigma_eff = 15.0

    return dict(
        Sigma_gas=Sigma_gas,
        Sigma_star=Sigma_star,
        Omega_d=Omega_d,
        H_star=H_star,
        sigma_eff=sigma_eff,
    )


@pytest.fixture(scope="module")
def toy_model_cgs():
    """Same toy galaxy in CGS."""
    R0, Re = 8.0, 3.5
    R = np.linspace(2, 16, 100)
    r0 = R / R0
    re = R / Re

    return dict(
        Sigma_gas=10.0 / r0 * _surf_cgs,
        Sigma_star=42.0 * np.exp(R0 / Re - re) * _surf_cgs,
        Omega_d=28.0 / r0 * _kms_kpc_cgs,
        H_star=245.0 * _pc_cgs,
        sigma_eff=15.0 * 1.0e5,
    )


# ---------------------------------------------------------------------------
# Low-level weight & pressure functions
# ---------------------------------------------------------------------------


class TestWeightFunctions:
    """Physics checks on the primitive CGS weight functions."""

    # Solar-circle CGS reference values
    Sg = 10.0 * _surf_cgs
    Ss = 42.0 * _surf_cgs
    Hs = 245.0 * _pc_cgs
    Hg = 100.0 * _pc_cgs
    Od = 28.0 * _kms_kpc_cgs
    se = 15.0 * 1.0e5

    def test_weight_gas_formula(self):
        """W_gas = 0.5 * pi * G * Sigma_gas^2"""
        expected = 0.5 * np.pi * _Gconst_cgs * self.Sg**2
        assert np.isclose(get_weight_gas(self.Sg), expected)

    def test_weight_gas_scales_as_sigma_squared(self):
        w1 = get_weight_gas(self.Sg)
        w2 = get_weight_gas(2.0 * self.Sg)
        assert np.isclose(w2 / w1, 4.0)

    def test_weight_star_formula(self):
        """W_star = pi * G * Sigma_gas * Sigma_star * H_gas / (H_gas + H_star)"""
        expected = np.pi * _Gconst_cgs * self.Sg * self.Ss * self.Hg / (self.Hg + self.Hs)
        assert np.isclose(get_weight_star(self.Sg, self.Hg, self.Ss, self.Hs), expected)

    def test_weight_star_approaches_thin_limit(self):
        """For H_gas >> H_star, W_star -> pi*G*Sigma_gas*Sigma_star (thin-disk limit)."""
        W = get_weight_star(self.Sg, 1e10 * self.Hg, self.Ss, self.Hs)
        thin = np.pi * _Gconst_cgs * self.Sg * self.Ss
        assert np.isclose(W, thin, rtol=1e-4)

    def test_weight_dm_formula(self):
        """W_dm = zeta_d * Sigma_gas * H_gas * Omega_d^2"""
        zeta_d = 1.0 / 3.0
        expected = zeta_d * self.Sg * self.Hg * self.Od**2
        assert np.isclose(get_weight_dm(self.Sg, self.Hg, self.Od), expected)

    def test_pressure_formula(self):
        """P = 0.5 * Sigma_gas / H_gas * sigma_eff^2"""
        expected = 0.5 * self.Sg / self.Hg * self.se**2
        assert np.isclose(get_pressure(self.Sg, self.Hg, self.se), expected)

    def test_pressure_weight_balance_at_equilibrium(self, toy_model_cgs):
        """At equilibrium H, total pressure should equal total weight."""
        d = toy_model_cgs
        # Use PRFM class to find equilibrium H, then check P = W_tot
        m = prfm.PRFM(**{k: v / (1e5 if k == "sigma_eff" else (_surf_cgs if "Sigma" in k else (_pc_cgs if "H" in k else _kms_kpc_cgs))) for k, v in d.items()}, astro_units=True)
        m.calc_weights()
        P = get_pressure(d["Sigma_gas"], m._H_gas, d["sigma_eff"])
        np.testing.assert_allclose(P, m._Wtot, rtol=1e-5)


# ---------------------------------------------------------------------------
# Scale height: analytic vs. numerical consistency
# ---------------------------------------------------------------------------


WEIGHT_COMBOS = [
    (1, 0, 0),  # gas only
    (0, 1, 0),  # star only
    (0, 0, 1),  # dm only
    (1, 1, 0),  # gas + star
    (1, 0, 1),  # gas + dm
    (1, 1, 1),  # full
]


@pytest.mark.parametrize("wgas,wstar,wdm", WEIGHT_COMBOS)
def test_analytic_numerical_scale_height_agreement(toy_model, wgas, wstar, wdm):
    """Analytic and numerical scale heights must agree to within 1% across the radial profile."""
    m = prfm.PRFM(**toy_model, astro_units=True)
    H_ana = m.get_scale_height(method="analytic", wgas=wgas, wstar=wstar, wdm=wdm)
    H_num = m.get_scale_height(method="numerical", wgas=wgas, wstar=wstar, wdm=wdm)
    np.testing.assert_allclose(H_ana, H_num, rtol=0.01, err_msg=f"wgas={wgas} wstar={wstar} wdm={wdm}")


# ---------------------------------------------------------------------------
# PRFM class: weights and weight contributions
# ---------------------------------------------------------------------------


class TestPRFMWeights:
    def test_calc_weights_sets_attributes(self, toy_model):
        m = prfm.PRFM(**toy_model, astro_units=True)
        m.calc_weights()
        for attr in ["H_gas", "Wgas", "Wstar", "Wdm", "Wtot"]:
            assert hasattr(m, attr), f"missing attribute {attr}"
            assert np.all(np.isfinite(getattr(m, attr))), f"{attr} contains non-finite values"

    def test_weight_contributions_sum_to_one(self, toy_model):
        m = prfm.PRFM(**toy_model, astro_units=True)
        m.calc_weights()
        fgas, fstar, fdm = m.get_weight_contribution()
        np.testing.assert_allclose(fgas + fstar + fdm, 1.0, rtol=1e-10)

    def test_weight_contributions_non_negative(self, toy_model):
        m = prfm.PRFM(**toy_model, astro_units=True)
        m.calc_weights()
        fgas, fstar, fdm = m.get_weight_contribution()
        assert np.all(fgas >= 0) and np.all(fstar >= 0) and np.all(fdm >= 0)

    def test_wtot_equals_sum_of_components(self, toy_model):
        m = prfm.PRFM(**toy_model, astro_units=True)
        m.calc_weights()
        np.testing.assert_allclose(m.Wgas + m.Wstar + m.Wdm, m.Wtot, rtol=1e-10)


# ---------------------------------------------------------------------------
# Self-consistent solutions
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "sigma_eff_model",
    [15.0, "tigress-classic-mid", "tigress-classic-avg"],
    ids=["constant", "tigress-classic-mid", "tigress-classic-avg"],
)
@pytest.mark.parametrize("method", ["analytic", "numerical"])
def test_self_consistent_solution(toy_model, sigma_eff_model, method):
    """Self-consistent solution converges and produces positive, finite outputs."""
    params = {**toy_model, "sigma_eff": sigma_eff_model}
    m = prfm.PRFM(**params, astro_units=True)
    H = m.calc_self_consistent_solution(method=method)
    assert np.all(np.isfinite(H)), "H_gas contains non-finite values"
    assert np.all(H > 0), "H_gas contains non-positive values"
    assert np.all(np.isfinite(m.Wtot))
    assert np.all(m.Wtot > 0)


def test_model_sigma_eff_consistent_with_pressure(toy_model):
    """For a sigma_eff model, sigma_eff stored in the PRFM object must match
    what get_sigma_eff returns for the converged pressure."""
    params = {**toy_model, "sigma_eff": "tigress-classic-mid"}
    m = prfm.PRFM(**params, astro_units=True)
    m.calc_self_consistent_solution()
    expected_sigma = get_sigma_eff(m._Wtot, model="tigress-classic-mid")
    np.testing.assert_allclose(m._sigma_eff, expected_sigma, rtol=1e-6)


# ---------------------------------------------------------------------------
# Feedback yield models
# ---------------------------------------------------------------------------


class TestFeedbackYield:
    # Reference pressure ~ 1e4 k_B cm^-3 K in CGS
    P_ref = 1.0e4 * _kbol_cgs

    def test_classic_yield_positive(self):
        Y = get_feedback_yield(self.P_ref, model="tigress-classic")
        assert Y > 0

    def test_decomp_yield_equals_sum_of_components(self):
        P = self.P_ref
        Ytot = get_feedback_yield(P, model="tigress-classic-decomp")
        Yth = get_feedback_yield_comp(P, comp="th", model="tigress-classic-decomp")
        Ytrb = get_feedback_yield_comp(P, comp="trb", model="tigress-classic-decomp")
        assert np.isclose(Ytot, Yth + Ytrb)

    def test_ncr_decomp_all_yield_positive(self):
        Y = get_feedback_yield(self.P_ref, Z=1.0, model="tigress-ncr-decomp-all")
        assert Y > 0

    def test_yield_decreases_with_pressure_classic(self):
        """Classic model has negative slope — higher P gives lower yield."""
        Y_lo = get_feedback_yield(self.P_ref, model="tigress-classic")
        Y_hi = get_feedback_yield(10 * self.P_ref, model="tigress-classic")
        assert Y_hi < Y_lo

    @pytest.mark.parametrize(
        "model",
        ["tigress-classic", "tigress-classic-decomp", "tigress-ncr", "tigress-ncr-decomp"],
    )
    def test_yield_model_returns_positive_array(self, model):
        P = self.P_ref * np.logspace(-1, 2, 20)
        kwargs = {"Z": 1.0} if "ncr" in model else {}
        Y = get_feedback_yield(P, model=model, **kwargs)
        assert np.all(Y > 0)


# ---------------------------------------------------------------------------
# SFR calculation
# ---------------------------------------------------------------------------


class TestSFRCalculation:
    def test_sfr_positive_constant_yield(self, toy_model):
        m = prfm.PRFM(**toy_model, astro_units=True)
        m.calc_self_consistent_solution()
        Sigma_SFR = m.calc_sfr(Ytot=1.0e3)
        assert np.all(Sigma_SFR > 0)

    @pytest.mark.parametrize(
        "ytot_model", ["tigress-classic", "tigress-classic-decomp"]
    )
    def test_sfr_positive_yield_model(self, toy_model, ytot_model):
        m = prfm.PRFM(**toy_model, astro_units=True)
        m.calc_self_consistent_solution()
        Sigma_SFR = m.calc_sfr(Ytot=ytot_model)
        assert np.all(Sigma_SFR > 0)
        assert np.all(np.isfinite(Sigma_SFR))

    def test_sfr_inversely_proportional_to_yield(self, toy_model):
        """Doubling Ytot should halve Sigma_SFR."""
        m = prfm.PRFM(**toy_model, astro_units=True)
        m.calc_self_consistent_solution()
        sfr1 = m.calc_sfr(Ytot=500.0).copy()
        sfr2 = m.calc_sfr(Ytot=1000.0).copy()
        np.testing.assert_allclose(sfr1 / sfr2, 2.0, rtol=1e-10)

    def test_sfr_function_agrees_with_class(self, toy_model):
        """get_sfr functional API should match PRFM.calc_sfr."""
        m = prfm.PRFM(**toy_model, astro_units=True)
        m.calc_self_consistent_solution()
        Sigma_SFR_class = m.calc_sfr(Ytot=1.0e3)
        Sigma_SFR_func = get_sfr(m._Wtot, Ytot=1.0e3) / _sfr_cgs
        np.testing.assert_allclose(Sigma_SFR_class, Sigma_SFR_func, rtol=1e-10)
