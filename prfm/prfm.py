from argparse import ArgumentError
import astropy.constants as ac
import astropy.units as au
import numpy as np
from scipy.optimize import brentq

# setting up some convenient conversion factors
_Gconst_cgs = ac.G.cgs.value
_pc_cgs = ac.pc.cgs.value
_kbol_cgs = ac.k_B.cgs.value

_surf_cgs = (ac.M_sun / ac.pc**2).cgs.value
_rho_cgs = (ac.M_sun / ac.pc**3).cgs.value
_sfr_cgs = (ac.M_sun / ac.kpc**2 / au.yr).cgs.value
_kms_cgs = (1 * au.km / au.s).cgs.value

_sigma_eff_models = dict()
_sigma_eff_models["tigress-classic-mid"] = dict(sigma_0=9.8, expo=0.15, sigma_min=5)
_sigma_eff_models["tigress-classic-avg"] = dict(sigma_0=12, expo=0.22, sigma_min=5)
_sigma_eff_models["tigress-ncr-mid"] = dict(
    sigma_0=8.9, expo=0.08, expo_Z=-0.005, sigma_min=5
)
_sigma_eff_models["tigress-ncr-avg"] = dict(
    sigma_0=11.7, expo=0.12, expo_Z=-0.03, sigma_min=5
)

_yield_models = dict()
# TIGRESS-NCR: Ostriker & Kim 2022
# Equations (25)
# renormalized for P/10^4
# turbulent is multiplied by 1.5 to accomodate magnetic
_yield_models["tigress-classic"] = dict(Y0=10 ** (3.86 - 0.212 * 4), expo=-0.212)
_yield_models["tigress-classic-decomp"] = dict(
    Yth0=10 ** (4.45 - 0.506 * 4),
    expo_th=-0.506,
    Ytrb0=1.5 * 10 ** (2.81 - 0.06 * 4),
    expo_trb=-0.06,
)
# TIGRESS-NCR: Kim et al. 2024
# Equations 14, 15, 16
_yield_models["tigress-ncr-decomp-all"] = dict(
    Yth0=390,
    expo_th=-0.46,
    expo_th_Z=-0.53,
    Ytrb0=561,
    expo_trb=-0.21,
    expo_trb_Z=-0.04,
    Ymag0=578,
    expo_mag=-0.40,
    expo_mag_Z=-0.44,
)
# Equations 14 and 17
_yield_models["tigress-ncr-decomp"] = dict(
    Yth0=390,
    expo_th=-0.46,
    expo_th_Z=-0.53,
    Ytrb0=1.17e3,
    expo_trb=-0.22,
    expo_trb_Z=-0.18,
)
# Equation 20
_yield_models["tigress-ncr"] = dict(Y0=1.65e3, expo=-0.29, expo_Z=-0.27)


def get_weight_gas(Sigma_gas):
    """weight due to gas self-gravity"""
    return 0.5 * np.pi * _Gconst_cgs * Sigma_gas**2


def get_weight_star(Sigma_gas, H_gas, Sigma_star, H_star):
    """weight due to stellar gravity"""
    return np.pi * _Gconst_cgs * Sigma_gas * Sigma_star * H_gas / (H_gas + H_star)


def get_weight_star_gaussian(Sigma_gas, H_gas, Sigma_star, H_star):
    """weight due to stellar gravity for Gaussian profiles"""
    return 2 * _Gconst_cgs * Sigma_gas * Sigma_star * np.arctan(H_gas / H_star)


def get_weight_dm(Sigma_gas, H_gas, Omega_d, zeta_d=1 / 3.0):
    """weight due to dark matter gravity"""
    return zeta_d * Sigma_gas * H_gas * Omega_d**2


def get_pressure(Sigma_gas, H_gas, sigma_eff):
    """total pressure"""
    return 0.5 * Sigma_gas / H_gas * sigma_eff**2


def get_weights(
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff, zeta_d=1 / 3.0, method="analytic"
):
    """calculate gas scale height and then each weight term

    Parameters
    ----------
    Sigma_gas: float or array-like
        gas surface density
    Sigma_star: float or array-like
        stellar surface density
    Omega_d: float or array-like
        galactic rotation speed
    H_star: float or array-like
        stellar scale height
    sigma_eff: str or float or array-like
        velocity dispersion
    zeta_d: float [1/3.]
        parameter for dm and gas distribution
    method: str [analytic or numerical]
        calculate scale height either analytically or numerically

    Returns
    -------
    H_gas: float or array-like
        gas scale height
    W_gas: float or array-like
        weight by gas
    W_star: float or array-like
        weight by star
    W_dm: float or array-like
        weight by dark matter

    Example
    -------
    >>> import prfm
    >>> H, W_gas, W_star, W_dm = prfm.get_weights(Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff)
    """
    H_gas = get_scale_height(
        Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff, zeta_d=zeta_d, method=method
    )

    W_gas = get_weight_gas(Sigma_gas)
    W_star = get_weight_star(Sigma_gas, H_gas, Sigma_star, H_star)
    W_dm = get_weight_dm(Sigma_gas, H_gas, Omega_d, zeta_d=zeta_d)

    return H_gas, W_gas, W_star, W_dm


def get_weight_star_thick(Sigma_gas, rho_star, sigma_eff):
    W_star = Sigma_gas * np.sqrt(2 * _Gconst_cgs * rho_star) * sigma_eff
    return W_star


def get_weights_thick(Sigma_gas, rho_star, sigma_eff):
    """calculate gas scale height and then each weight term for H_*>>H_gas

    Parameters
    ----------
    Sigma_gas: float or array-like
        gas surface density
    rho_star: float or array-like
        stellar volume density
    sigma_eff: str or float or array-like
        velocity dispersion

    Returns
    -------
    H_gas: float or array-like
        gas scale height
    W_gas: float or array-like
        weight by gas
    W_star: float or array-like
        weight by star
    W_dm: float or array-like
        weight by dark matter

    Note
    ----
    Used approximate formula Equation (7) in Ostriker and Kim (2022)
    """
    H_gas = get_scale_height_thick(Sigma_gas, rho_star, sigma_eff)

    W_gas = get_weight_gas(Sigma_gas)
    W_star = get_weight_star_thick(Sigma_gas, rho_star, sigma_eff)
    W_dm = 0.0 * H_gas / H_gas

    return H_gas, W_gas, W_star, W_dm


def get_weight_star_thin(Sigma_gas, Sigma_star, sigma_eff):
    W_star = np.pi * _Gconst_cgs * Sigma_gas * Sigma_star
    return W_star


def get_weights_thin(Sigma_gas, Sigma_star, sigma_eff):
    """calculate gas scale height and then each weight term for H_*>>H_gas

    Parameters
    ----------
    Sigma_gas: float or array-like
        gas surface density
    Sigma_star: float or array-like
        stellar surface density
    sigma_eff: str or float or array-like
        velocity dispersion

    Returns
    -------
    H_gas: float or array-like
        gas scale height
    W_gas: float or array-like
        weight by gas
    W_star: float or array-like
        weight by star
    W_dm: float or array-like
        weight by dark matter

    Note
    ----
    Used approximate formula Equation (7) in Ostriker and Kim (2022)
    """
    H_gas = get_scale_height_thin(Sigma_gas, Sigma_star, sigma_eff)

    W_gas = get_weight_gas(Sigma_gas)
    W_star = get_weight_star_thin(Sigma_gas, Sigma_star, sigma_eff)
    W_dm = 0.0 * H_gas / H_gas

    return H_gas, W_gas, W_star, W_dm


def get_weight_contribution(
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff, zeta_d=1 / 3.0, method="analytic"
):
    """Wrapper function to calculate the ratio of each weight term to the total weight

    Parameters
    ----------
    Sigma_gas: float or array-like
        gas surface density
    Sigma_star: float or array-like
        stellar surface density
    Omega_d: float or array-like
        galactic rotation speed
    H_star: float or array-like
        stellar scale height
    sigma_eff: str or float or array-like
        velocity dispersion
    zeta_d: float [1/3.]
        parameter for dm and gas distribution
    method: str [analytic or numerical]
        calculate scale height either analytically or numerically

    Returns
    -------
    f_gas: float or array-like
        weight by gas/total weight
    f_star: float or array-like
        weight by star/total weight
    f_dm: float or array-like
        weight by dark matter/total weight
    """
    H, wgas, wstar, wdm = get_weights(
        Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff, zeta_d=zeta_d, method=method
    )
    wtot = wgas + wstar + wdm

    return wgas / wtot, wstar / wtot, wdm / wtot


def get_scale_height_gas_only(*args, **kwargs):
    """analytic solution for gas only case"""
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff = args
    return sigma_eff**2 / (np.pi * _Gconst_cgs * Sigma_gas)


def get_scale_height_star_only(*args, **kwargs):
    """analytic solution for star only case"""
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff = args
    h = sigma_eff**2 / (4 * np.pi * _Gconst_cgs * Sigma_star)
    h_star = H_star / h
    return h * (1 + np.sqrt(1 + 2 * h_star))


def get_scale_height_dm_only(*args, **kwargs):
    """analytic solution for dark matter only case"""
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff = args
    zeta_d = kwargs["zeta_d"]
    return sigma_eff / np.sqrt(2 * zeta_d) / Omega_d


def get_scale_height_star_gas(*args, **kwargs):
    """analytic solution neglecting dark matter"""
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff = args
    H_gas_only = get_scale_height_gas_only(*args, **kwargs)
    eta_star = H_star / H_gas_only
    s_star = Sigma_star / Sigma_gas
    # h = ((1-eta_star)+np.sqrt((eta_star-1)**2+4*eta_star*(1+2.*s_star)))/(2*(1+2*s_star))
    h = (
        2
        * eta_star
        / ((eta_star - 1) + np.sqrt((eta_star + 1) ** 2 + 8 * s_star * eta_star))
    )
    return h * H_gas_only


def get_scale_height_star_gas_approx(*args, **kwargs):
    """analytic solution neglecting dark matter"""
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff = args
    H_gas_only = get_scale_height_gas_only(*args, **kwargs)
    eta_star = H_star / H_gas_only
    s_star = Sigma_star / Sigma_gas
    h = (1 + np.sqrt(1 + 8 * s_star * eta_star)) / (4 * s_star)
    return h * H_gas_only


def get_scale_height_dm_gas(*args, **kwargs):
    """analytic solution neglectic star"""
    H_gas_only = get_scale_height_gas_only(*args, **kwargs)
    H_dm_only = get_scale_height_dm_only(*args, **kwargs)
    eta_D = H_dm_only / H_gas_only
    h = (-(eta_D**2) + np.sqrt(eta_D**4 + 4 * eta_D**2)) * 0.5
    return h * H_gas_only


def get_sigma_eff(P, Z=None, model="tigress-classic-mid"):
    """pressure-dependent velocity dispersion model"""
    sigma_eff_model = _sigma_eff_models[model]
    sigma_0 = sigma_eff_model["sigma_0"]
    expo = sigma_eff_model["expo"]
    sigma_min = sigma_eff_model["sigma_min"]

    veld = sigma_0 * (1.0e-4 * P / _kbol_cgs) ** expo
    if (Z is not None) and ("expo_Z" in sigma_eff_model):
        veld *= Z ** sigma_eff_model["expo_Z"]
    return np.clip(veld, sigma_min, None) * 1.0e5


def get_feedback_yield(P, Z=None, model="tigress-classic"):
    """total feedback yield as a function of weight

    Paramerters
    -----------
    P : float
        midplane pressure or weight
    model : str
        tigress-classic for OK22
        tigress-NCR for K24

    Returns
    -------
    Ytot : float
        feedback yield in km/s
    """
    yield_model = _yield_models[model]
    if "Y0" in yield_model:
        Y0 = yield_model["Y0"]
        slope = yield_model["expo"]
        Ytot = Y0 * (1.0e-4 * P / _kbol_cgs) ** slope
        if (Z is not None) and ("expo_Z" in yield_model):
            Ytot *= Z ** yield_model["expo_Z"]
        return Ytot
    elif "Yth0" in yield_model:
        Yth0 = yield_model["Yth0"]
        expo_th = yield_model["expo_th"]
        Ytrb0 = yield_model["Ytrb0"]
        expo_trb = yield_model["expo_trb"]
        Yth = Yth0 * (1.0e-4 * P / _kbol_cgs) ** expo_th
        if (Z is not None) and ("expo_th_Z" in yield_model):
            Yth *= Z ** yield_model["expo_th_Z"]
        Ytrb = Ytrb0 * (1.0e-4 * P / _kbol_cgs) ** expo_trb
        if (Z is not None) and ("expo_trb_Z" in yield_model):
            Ytrb *= Z ** yield_model["expo_trb_Z"]
        if "Ymag0" in yield_model:
            Ymag0 = yield_model["Ymag0"]
            expo_mag = yield_model["expo_mag"]
            expo_mag_Z = yield_model["expo_mag_Z"]
            # add magnetic contribution to turbulent for fully decomposed yield
            Ytrb += Ymag0 * (1.0e-4 * P / _kbol_cgs) ** expo_mag * Z**expo_mag_Z
        return Yth + Ytrb


def get_feedback_yield_comp(P, Z=None, comp="th", model="tigress-classic-decomp"):
    """feedback yield of each component as a function of weight

    Paramerters
    -----------
    P : float
        midplane pressure or weight
    comp : str
        component name ['th','trb']
    model : str
        tigress-classic for OK22
        tigress-NCR for K23

    Returns
    -------
    Y : float
        feedback yield in km/s
    """
    yield_model = _yield_models[model]
    Y0 = yield_model[f"Y{comp}0"]
    expo = yield_model[f"expo_{comp}"]
    Y = Y0 * (1.0e-4 * P / _kbol_cgs) ** expo
    if (Z is not None) and (f"expo_{comp}_Z" in yield_model):
        Y *= Z ** yield_model[f"expo_{comp}_Z"]
    return Y


def get_sfr(P, Z=None, Ytot="tigress-classic"):
    """calculate SFR surface density using feedback yield

    Paramerters
    -----------
    P : float
        midplane pressure or weight
    Ytot : float, str
        constant value in km/s or
        P-dependent feedback yield model e.g., `tigress-classic`

    Returns
    -------
    sfr : float
        SFR surface density in g/cm^2/yr
    """
    Ytot = get_feedback_yield(P, Z=Z, model=Ytot) if isinstance(Ytot, str) else Ytot
    return P / (Ytot * 1.0e5)


def get_scale_height(
    Sigma_gas,
    Sigma_star,
    Omega_d,
    H_star,
    sigma_eff,
    zeta_d=1 / 3.0,
    method="analytic",
    wgas=1,
    wstar=1,
    wdm=1,
):
    """wrapper function to calculate gas scale height either numerically or analytically

    Parameters
    ----------
    Sigma_gas : float or array-like
        gas surface density
    Sigma_star : float or array-like
        stellar surface density
    Omega_d : float or array-like
        galactic rotation speed
    H_star : float or array-like
        stellar scale height
    sigma_eff : str or float or array-like
        velocity dispersion
    zeta_d : float [1/3.]
        parameter for dm and gas distribution
    method : str [analytic or numerical]
        calculate scale height either analytically or numerically

    wgas : int [0 or 1]
        toggle weight term from gas
    wstar : int [0 or 1]
        toggle weight term from stars
    wdm : int [0 or 1]
        toggle weight term from dark matter

    Returns
    -------
    H : float or array-like
        scale height in cm

    Example
    -------
    >>> import prfm
    >>> H=prfm.get_scale_height(Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff)
    """
    args = (Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff)
    kwargs = dict(zeta_d=zeta_d, wgas=wgas, wstar=wstar, wdm=wdm)
    w = (wgas << 2) + (wstar << 1) + (wdm << 0)

    if method == "numerical":
        return get_scale_height_numerical(*args, **kwargs)
    else:
        if w == 7:  # 111
            return get_scale_height_analytic(*args, zeta_d=zeta_d)
        elif w == 6:  # 110
            return get_scale_height_star_gas(*args, zeta_d=zeta_d)
        elif w == 5:  # 101
            return get_scale_height_dm_gas(*args, zeta_d=zeta_d)
        elif w == 4:  # 100
            return get_scale_height_gas_only(*args, zeta_d=zeta_d)
        elif w == 3:  # 011
            # no analytic solution provided
            return
        elif w == 2:  # 010
            return get_scale_height_star_only(*args, zeta_d=zeta_d)
        elif w == 1:  # 001
            return get_scale_height_dm_only(*args, zeta_d=zeta_d)
        else:
            # no such gas is expected
            raise ("at least one term must be considered")


def get_scale_height_thick(Sigma_gas, rho_star, sigma_eff, zeta_d=1 / np.pi):
    """approximate solution for H_* >> H_gas without explicitly defining H_*

    Parameters
    ----------
    Sigma_gas : float or array-like
        gas surface density
    rho_star : float or array-like
        stellar(+dark matter) volume density
    sigma_eff : str or float or array-like
        velocity dispersion

    Returns
    -------
    H : float or array-like
        scale height in cm

    Notes
    -----
    Ostriker & Kim (2022) Equation (5)
    """
    g1 = np.pi * _Gconst_cgs * Sigma_gas
    g2 = 32 * np.pi * zeta_d * _Gconst_cgs * rho_star * sigma_eff**2
    return 2 * sigma_eff**2 / (g1 + np.sqrt(g1**2 + g2))


def get_scale_height_thin(Sigma_gas, Sigma_star, sigma_eff):
    """approximate solution for H_* << H_gas without explicitly defining H_*

    Parameters
    ----------
    Sigma_gas : float or array-like
        gas surface density
    Sigma_star : float or array-like
        stellar surface density
    sigma_eff : str or float or array-like
        velocity dispersion

    Returns
    -------
    H : float or array-like
        scale height in cm

    Notes
    -----
    Ostriker & Kim (2022) Equation (5)
    """
    g1 = np.pi * _Gconst_cgs * Sigma_gas
    g2 = 2 * np.pi * _Gconst_cgs * Sigma_star
    return sigma_eff**2 / (g1 + g2)


# @np.vectorize
def get_scale_height_analytic(
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff, zeta_d=1 / 3.0
):
    """Analytic solution of the cubic equation for the vertical dynamical equilibirum
    including all three weight terms.

    All inputs must be in c.g.s. units
    """
    args = (Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff)
    kwargs = dict(zeta_d=zeta_d)
    H_gas_only = get_scale_height_gas_only(*args, **kwargs)
    H_dm_only = get_scale_height_dm_only(*args, **kwargs)

    s_star = Sigma_star / Sigma_gas
    eta_star = H_star / H_gas_only
    eta_dm_sq = (H_dm_only / H_gas_only) ** 2

    a0 = -eta_star * eta_dm_sq
    a1 = (eta_star - 1) * eta_dm_sq
    a2 = (1 + 2 * s_star) * eta_dm_sq + eta_star

    Q = (3 * a1 - a2**2) / 9.0
    R = (9 * a2 * a1 - 27 * a0 - 2 * a2**3) / 54.0

    theta = np.arccos(R / np.sqrt(-(Q**3)))
    h = 2 * np.sqrt(-Q) * np.cos(theta / 3) - a2 / 3.0

    return h * H_gas_only


@np.vectorize
def get_scale_height_numerical(
    Sigma_gas,
    Sigma_star,
    Omega_d,
    H_star,
    sigma_eff,
    zeta_d=1 / 3.0,
    wgas=1,
    wstar=1,
    wdm=1,
    fgas=get_weight_gas,
    fstar=get_weight_star,
    fdm=get_weight_dm,
):
    """Numerical solution of the vertical dynamical equilibrium equation
    including all weight terms.

    All inputs must be in c.g.s. units

    Parameters
    ----------
    fgas : func
    fstar : func
    fdm : func
    """

    def fun(x):
        return (
            get_pressure(Sigma_gas, x, sigma_eff)
            - wgas * fgas(Sigma_gas)
            - wstar * fstar(Sigma_gas, x, Sigma_star, H_star)
            - wdm * fdm(Sigma_gas, x, Omega_d, zeta_d=zeta_d)
        )

    # soln=root(fun,H_gas_init)
    # return soln.x
    return brentq(fun, 1.0e-3 * _pc_cgs, 1.0e6 * _pc_cgs)


def get_self_consistent_solution(
    Sigma_gas,
    Sigma_star,
    Omega_d,
    H_star,
    sigma_eff="tigress-claasic-mid",
    zeta_d=1 / 3.0,
    method="analytic",
    tol=1.0e-5,
    niter=16,
):
    """Iteratively solve for sigma_eff(P) to get converged weight, H, and sigma_eff

    Parameters
    ----------
    Sigma_gas: float or array-like
        gas surface density
    Sigma_star: float or array-like
        stellar surface density
    Omega_d: float or array-like
        galactic rotation speed
    H_star: float or array-like
        stellar scale height
    sigma_eff : float or str
        effective vertical velocity dispersion
    zeta_d: float [1/3.]
        parameter for dm and gas distribution
    method: str [analytic or numerical]
        calculate scale height either analytically or numerically

    Returns
    -------
    Wtot : float or array
    H : float or array
    sigma_eff: float or array
    """
    # make initial guess
    if not isinstance(sigma_eff, str):
        sigma0 = sigma_eff
    else:
        if sigma_eff not in _sigma_eff_models:
            raise IndexError("sigma_eff models: ", list(_sigma_eff_models.keys()))
        sigma0 = 15.0e5
    args = (Sigma_gas, Sigma_star, Omega_d, H_star, sigma0)
    kwargs = dict(zeta_d=zeta_d, method=method)
    Hprev, wgas, wstar, wdm = get_weights(*args, **kwargs)
    wtot_prev = wgas + wstar + wdm

    # return if velocity diseprsion is a constant
    if not isinstance(sigma_eff, str):
        return wtot_prev, Hprev, sigma0 * np.ones_like(Hprev)

    # iterative solve
    for i in range(niter):
        new_sigma_eff = get_sigma_eff(wtot_prev, model=sigma_eff)
        args = (Sigma_gas, Sigma_star, Omega_d, H_star, new_sigma_eff)
        Hnext, wgas, wstar, wdm = get_weights(*args, **kwargs)
        wtot_next = wgas + wstar + wdm
        # L1_norm = np.sum(np.abs(wtot_next/wtot_prev - 1))
        L1_norm = np.sum(np.abs(Hnext / Hprev - 1))
        if L1_norm < tol:
            break
        wtot_prev = np.copy(wtot_next)
        Hprev = np.copy(Hnext)

    return wtot_next, Hnext, get_sigma_eff(wtot_next, model=sigma_eff)


class PRFM(object):
    """a wrapper class for scale height, weight, SFR calculations

    Parameters
    ----------
    Sigma_gas : float, array-like
        gas surface density
    Sigma_star : float, array-like
        stellar surface density
    Omega_d : float, array-like
        galactic rotation speed
    H_star : float, array-like
        stellar scale height
    sigma_eff : str, float, array-like
        velocity dispersion
    zeta_d : float [1/3.]
        parameter for dm and gas distribution
    method : str [analytic or numerical]
        calculate scale height either analytically or numerically
    """

    def __init__(
        self,
        Sigma_gas=None,
        Sigma_star=None,
        H_star=None,
        rho_star=None,
        Omega_d=None,
        rho_dm=None,
        sigma_eff="tigress-classic-mid",
        Ytot="tigress-classic",
        astro_units=True,
    ):
        # initialize method attribute
        self._method = "analytic"

        # set model name based on sigma_eff
        self.name = sigma_eff if isinstance(sigma_eff, str) else "constant"

        # initialize unit conversion factors
        self.units = self._set_units()

        # setting up gas surface density
        if Sigma_gas is None:
            # gas density must be given
            raise ArgumentError("Sigma_gas must be provided")
        self._Sigma_gas = Sigma_gas

        # setting up stellar disk
        if (Sigma_star is None) and (rho_star is None):
            self._Sigma_star = 0
            self._rho_star = 0
            self._wstar = 0
        else:
            self._wstar = 1
            if H_star is None:
                if rho_star is None:
                    # this is limit for H_star << H_gas
                    self._stellar_disk = "thick"
                    self._Sigma_star = Sigma_star
                elif Sigma_star is None:
                    # this is limit for H_star >> H_gas
                    self._stellar_disk = "thick"
                    self._rho_star = rho_star
                else:
                    # both given
                    self._stellar_disk = "general"
                    self._H_star = Sigma_star / (2 * rho_star)
            else:
                if rho_star is None:
                    # Sigma_star is given
                    rho_star = Sigma_star / (2 * H_star)
                elif Sigma_star is None:
                    # rho_star is given
                    Sigma_star = rho_star * 2 * H_star
                else:
                    # both given
                    raise ArgumentError(
                        "cannot provide all three: Sigma_star, rho_star, H_star"
                    )

                self._stellar_disk = "general"
                self._Sigma_star = Sigma_star
                self._rho_star = rho_star
                self._H_star = H_star
        # unit conversion done before DM part
        if astro_units:
            self._astro_to_cgs()

        # setting up dark matter
        if (Omega_d is None) and (rho_dm is None):
            self._Omega_d = 0
            self._rho_dm = 0
            self._wdm = 0
        else:
            self._wdm = 1
            # this conversion assumes flat rotation
            _fourpiG = 4 * np.pi * ac.G.cgs.value
            if rho_dm is None:
                # Omega_d is given
                if astro_units:
                    Omega_d = Omega_d * self.units["Omega_d"].cgs.value
                rho_dm = Omega_d**2 / _fourpiG
            elif Omega_d is None:
                # rho_dm is given
                if astro_units:
                    rho_dm = rho_dm * self.units["rho_dm"].cgs.value
                Omega_d = np.sqrt(_fourpiG * rho_dm)
            else:
                # both given
                raise ArgumentError("cannot provide both Omega_d and rho_dm")

            self._Omega_d = Omega_d
            self._rho_dm = rho_dm

        self.set_feedback_yield(Ytot)
        self.set_sigma_eff(sigma_eff, astro_units=astro_units)

        # store parameters in astro-friendly units
        self._cgs_to_astro()

        # set weight functions
        self._weight_functions = [get_weight_gas, get_weight_star, get_weight_dm]

    def set_sigma_eff(self, sigma_eff, astro_units=True):
        """For plug and play"""
        if not isinstance(sigma_eff, str):
            if astro_units:
                u = self.units["sigma_eff"].cgs.value
                sigma_eff = sigma_eff * u
            self._sigma_eff = sigma_eff
            self._sigma_eff_model = "constant"
        else:
            if sigma_eff not in _sigma_eff_models:
                print(
                    f"{sigma_eff} effective EoS is not supported. "
                    f"Try one of {list(_sigma_eff_models.keys())}"
                )
            self._sigma_eff_model = sigma_eff
        self.updated = dict(SFR=False, H=False)
        # set arguments that will be passed to the functions
        self.reset_arg_list(sigma_eff)

    def set_feedback_yield(self, Ytot):
        """For plug and play"""
        if not isinstance(Ytot, str):
            self._yield_model = "constant"
        else:
            self._yield_model = Ytot
        self.updated = dict(SFR=False, H=False)
        self._Ytot = Ytot

    def set_weight_function_gas(self, func):
        """Update weight function for gas"""
        self._weight_functions[0] = func
        self._method = "numerical"

    def set_weight_function_star(self, func):
        """Update weight function for star"""
        self._weight_functions[1] = func
        self._method = "numerical"

    def set_weight_function_dm(self, func):
        """Update weight function for star"""
        self._weight_functions[2] = func
        self._method = "numerical"

    def reset_arg_list(self, sigma_eff):
        # reset arguments with new velocity dispersion
        if self._stellar_disk == "general":
            self._args = (
                self._Sigma_gas,
                self._Sigma_star,
                self._Omega_d,
                self._H_star,
                sigma_eff,
            )
        elif self._stellar_disk == "thick":
            self._args = (self._Sigma_gas, self._rho_star, sigma_eff)
        elif self._stellar_disk == "thin":
            self._args = (self._Sigma_gas, self._Sigma_star, sigma_eff)

    def _astro_to_cgs(self):
        for var in [
            "Sigma_gas",
            "Sigma_star",
            "H_star",
            "rho_star",
            "rho_dm",
            "Omega_d",
            "sigma_eff",
        ]:
            if hasattr(self, "_" + var):
                u_cgs = self.units[var].cgs.value
                v = getattr(self, "_" + var) * u_cgs
                setattr(self, "_" + var, v)

    def _cgs_to_astro(self):
        for var in [
            "Sigma_gas",
            "Sigma_star",
            "H_star",
            "rho_star",
            "rho_dm",
            "Omega_d",
            "sigma_eff",
        ]:
            if hasattr(self, "_" + var):
                u_cgs = self.units[var].cgs.value
                v = getattr(self, "_" + var) / u_cgs
                setattr(self, var, v)

    def _set_units(self):
        units = dict()
        units["Sigma_gas"] = units["Sigma_star"] = ac.M_sun / ac.pc**2
        units["rho_star"] = units["rho_dm"] = ac.M_sun / ac.pc**3
        units["Omega_d"] = au.km / au.s / ac.kpc
        units["H_star"] = units["H_gas"] = ac.pc
        units["Wtot"] = ac.k_B / au.cm**3 * au.K
        units["Wgas"] = units["Wstar"] = units["Wdm"] = units["Wtot"]
        units["Ytot"] = 1.0 * au.km / au.s
        units["sigma_eff"] = 1.0 * au.km / au.s
        units["Sigma_SFR"] = ac.M_sun / ac.kpc**2 / au.yr

        return units

    def __repr__(self) -> str:
        args = "PRFM calculator is prepared for\n"
        for var in [
            "Sigma_gas",
            "Sigma_star",
            "rho_star",
            "Omega_d",
            "rho_dm",
            "H_star",
            "sigma_eff",
        ]:
            if not hasattr(self, "_" + var):
                continue
            v = np.atleast_1d(getattr(self, var))
            args += "  {}[0]: {}, N={}\n".format(var, v[0], len(v))
        if hasattr(self, "_sigma_eff_model"):
            args += "  sigma_eff_model: {}\n".format(self._sigma_eff_model)
        if hasattr(self, "_yield_model"):
            args += "  feedback_yield_model: {}\n".format(self._yield_model)
        return args

    def get_scale_height(self, method=None, wgas=1, wstar=1, wdm=1):
        """Wrapper function to get the gas scale height

        Parameters
        ----------
        method : str
            override method if provided
        """

        if method is None:
            method = self._method

        if self._stellar_disk == "general":
            H_gas = get_scale_height(*self._args, wgas=wgas, wstar=wstar, wdm=wdm)
        elif self._stellar_disk == "thick":
            H_gas = get_scale_height_thick(*self._args)
        elif self._stellar_disk == "thin":
            H_gas = get_scale_height_thin(*self._args)

        return H_gas / self.units["H_gas"].cgs.value

    def get_scale_height_numerical(self):
        """Wrapper function to get the gas scale height for generic weight functions"""

        H_gas = get_scale_height_numerical(
            *self._args,
            fgas=self._weight_functions[0],
            fstar=self._weight_functions[1],
            fdm=self._weight_functions[2],
        )
        return H_gas / self.units["H_gas"].cgs.value

    def get_weight_contribution(self):
        """Wrapper function to get weight contribution"""
        if not hasattr(self, "Wtot"):
            self.calc_weights()
        return self.Wgas / self.Wtot, self.Wstar / self.Wtot, self.Wdm / self.Wtot

    def calc_weights(self, method=None):
        """Wrapper function to calculate scale height and all weights

        Parameters
        ----------
        method : str
            override method if provided
        """

        if method is None:
            method = self._method

        if method == "analytic":
            H_gas = self.get_scale_height(wstar=self._wstar, wdm=self._wdm)
        elif method == "numerical":
            H_gas = self.get_scale_height_numerical()

        _H_gas = H_gas * self.units["H_gas"].cgs.value

        self.H_gas = H_gas
        self._H_gas = _H_gas

        W_gas = self._weight_functions[0](self._Sigma_gas)
        if self._stellar_disk == "thick":
            W_star = get_weight_star_thick(*self._args)
        elif self._stellar_disk == "thin":
            W_star = get_weight_star_thin(*self._args)
        else:
            W_star = self._weight_functions[1](
                self._Sigma_gas, _H_gas, self._Sigma_star, self._H_star
            )
        W_dm = self._weight_functions[2](self._Sigma_gas, _H_gas, self._Omega_d)

        for var, v in zip(
            ["H_gas", "Wgas", "Wstar", "Wdm"], [_H_gas, W_gas, W_star, W_dm]
        ):
            u_cgs = self.units[var].cgs.value
            setattr(self, "_" + var, v)  # results in cgs
            setattr(self, var, v / u_cgs)  # results in astro units

        # total weight
        Wtot = W_gas + W_star + W_dm
        u_cgs = self.units["Wtot"].cgs.value
        self._Wtot = Wtot
        self.Wtot = Wtot / u_cgs

        if self._sigma_eff_model != "constant":
            self._sigma_eff = get_sigma_eff(Wtot, model=self._sigma_eff_model)
            self.sigma_eff = self._sigma_eff / self.units["sigma_eff"].cgs.value

    def calc_self_consistent_solution(
        self, method=None, niter=16, tol=1.0e-6, verbose=False
    ):
        """Wrapper function to calculate self-consistent solutions

        Update H_gas, W_tot (and all components), sigma_eff
        """
        if method is None:
            method = self._method

        print(
            f"Calculating {method} H and weight solutions "
            f"using the {self._sigma_eff_model} effective EoS"
        )

        # for model sigma, need to give an initial guess
        if self._sigma_eff_model != "constant":
            self.reset_arg_list(15.0e5)
        self.calc_weights(method=method)

        # that's it for constant sigma_eff
        # iterative solve for model sigma
        if self._sigma_eff_model != "constant":
            # iterative solve
            for i in range(niter):
                Hprev = np.copy(self._H_gas)
                self.reset_arg_list(self._sigma_eff)
                self.calc_weights(method=method)
                L1_norm = np.sum(np.abs(self._H_gas / Hprev - 1))
                if verbose:
                    print(i, L1_norm)
                if L1_norm < tol:
                    break
        self.updated["H"] = True
        return self.H_gas

    def calc_sfr(self, Z=None, Ytot=None):
        """Wrapper function to calculate SFR surface density

        Update Sigma_SFR (astro units) and _Sigma_SFR (cgs)
        """
        if not self.updated["H"]:
            self.calc_self_consistent_solution()

        if Ytot is None:
            Ytot = self._Ytot

        print(
            f"Calculating Sigma_SFR with the {self._yield_model} feedback yield model"
        )

        if Z is not None:
            print(" also considering metallicity dependence!")
        # results in c.g.s.
        Sigma_SFR = get_sfr(self._Wtot, Z=Z, Ytot=Ytot)
        u_cgs = self.units["Sigma_SFR"].cgs.value
        self._Sigma_SFR = Sigma_SFR
        self.Sigma_SFR = Sigma_SFR / u_cgs

        self.updated["SFR"] = True
        return self.Sigma_SFR

    def calc_times(self, Z=None, Ytot=None):
        if not self.updated["SFR"]:
            if Ytot is None:
                Ytot = self._Ytot
            self.calc_sfr(Ytot, Z=Z)

        # calculate tdep in c.g.s.
        self._tdep = self._Sigma_gas / self._Sigma_SFR
        self._tdyn = 2 * self._H_gas / self._sigma_eff

        self.tdep = (self._tdep * au.s).to("Gyr").value
        self.tdyn = (self._tdyn * au.s).to("Gyr").value

        return self.tdyn, self.tdep

    def check_solutions(self):
        self.calc_self_consistent_solution(method="analytic")
        sola = np.copy(self.Wtot)
        self.calc_self_consistent_solution(method="numerical")
        soln = np.copy(self.Wtot)

        L1 = np.mean(np.abs(sola - soln))
        print("difference between analytic and numerical solutions: {}".format(L1))
        return sola, soln
