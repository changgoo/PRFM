import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

import numpy as np
import astropy.units as au
import astropy.constants as ac
import pandas as pd
import xarray as xr
import cmasher as cmr

import os
from prfm import get_sfr, get_feedback_yield_comp

dirpath = os.path.dirname(__file__)

# setting up some convenient conversion factors
_yield_conv = (
    ((ac.k_B * au.K / au.cm**3) / (ac.M_sun / ac.kpc**2 / au.yr)).to("km/s").value
)
_kbol_cgs = ac.k_B.cgs.value
_surf_cgs = (ac.M_sun / ac.pc**2).cgs.value
_sfr_cgs = (ac.M_sun / ac.kpc**2 / au.yr).cgs.value
_kms_cgs = (1 * au.km / au.s).cgs.value

# =====================================================================================
# data container class
# =====================================================================================


class PRFM_data(object):
    """Simulation data container class."""

    def __init__(self, paper):
        self.paper = paper
        self.field_list = ["SFR", "Pturb", "Pth", "oPimag", "dPimag"]

    def convert_log_linear(self):
        """convert data in log10 into linear"""
        for f in self.field_list:
            logf = f"log_{f}"
            ferr = f"{f}_std"
            logferr = f"{logf}_std"
            if hasattr(self, logf):
                v = getattr(self, logf)
                setattr(self, f, 10**v)
                if hasattr(self, logferr):
                    verr = getattr(self, logferr)
                    setattr(self, ferr, 10**v * np.log(10) * verr)

    def convert_linear_log(self):
        """convert data in linear into log10"""
        for f in self.field_list:
            logf = f"log_{f}"
            ferr = f"{f}_std"
            logferr = f"{logf}_std"
            if hasattr(self, f):
                v = getattr(self, f)
                setattr(self, logf, np.log10(v))
                if hasattr(self, ferr):
                    verr = getattr(self, ferr)
                    setattr(self, logferr, verr / v / np.log(10))

    def get_Ptotal(self):
        """calculate *Ptot*, *Pnonth*, *Pimag* given `Pturb`, `Pth`, `oPimag`, `dPimag`"""
        self.Ptot = np.zeros_like(self.Pturb)
        self.Ptot_std = np.zeros_like(self.Pturb)
        self.Pnonth = np.zeros_like(self.Pturb)
        self.Pnonth_std = np.zeros_like(self.Pturb)

        # add thermal and turbulent
        for f in ["Pturb", "Pth"]:
            ferr = f"{f}_std"
            self.Ptot += getattr(self, f)
            self.Ptot_std += getattr(self, ferr) ** 2
            if f != "Pth":
                self.Pnonth += getattr(self, f)
                self.Pnonth_std += getattr(self, ferr) ** 2
        # get total mag first
        if hasattr(self, "oPimag"):
            self.Pimag = self.oPimag + self.dPimag
            self.Pimag_std = np.sqrt(self.oPimag_std**2 + self.dPimag_std**2)

            v = self.Pimag
            verr = self.Pimag_std
        # add magnetic term
        if hasattr(self, "Pimag"):
            self.Ptot += np.nan_to_num(getattr(self, "Pimag"))
            self.Ptot_std += np.nan_to_num(getattr(self, "Pimag_std") ** 2)
            self.Pnonth += np.nan_to_num(getattr(self, "Pimag"))
            self.Pnonth_std += np.nan_to_num(getattr(self, "Pimag_std") ** 2)
        self.Ptot_std = np.sqrt(self.Ptot_std)
        self.Pnonth_std = np.sqrt(self.Pnonth_std)

        # get log
        v = self.Ptot
        verr = self.Ptot_std
        self.log_Ptot = np.log10(v)
        self.log_Ptot_std = verr / v / np.log(10)
        self.field_list += ["Ptot"]

        # non thermal
        v = self.Pnonth
        verr = self.Pnonth_std
        self.log_Pnonth = np.log10(v)
        self.log_Pnonth_std = verr / v / np.log(10)
        self.field_list += ["Pnonth"]

        if hasattr(self, "Pimag"):
            v = self.Pimag
            verr = self.Pimag_std
            self.log_Pimag = np.log10(v)
            self.log_Pimag_std = verr / v / np.log(10)
            if "Pimag" not in self.field_list:
                self.field_list += ["Pimag"]

    def get_yield(self):
        """calculate feedback yield P/SFR in km/s"""
        SFR = self.SFR
        SFRerr = self.SFR_std
        for f in self.field_list:
            if f in ["SFR"]:
                continue
            if not hasattr(self, f):
                continue
            P = getattr(self, f)
            Perr = getattr(self, f"{f}_std")
            Y = P / SFR * _yield_conv
            Yerr = np.sqrt((Perr / P) ** 2 + (SFRerr / SFR) ** 2) * Y * _yield_conv
            setattr(self, f"Y_{f}", Y)
            setattr(self, f"Y_{f}_std", Yerr)
            setattr(self, f"log_Y_{f}", np.log10(Y))
            setattr(self, f"log_Y_{f}_std", Yerr / Y / np.log(10))


# =====================================================================================
# plotting utilities
# =====================================================================================


def add_one_sim(data, xf="Ptot", yf="SFR", ms=5):
    c = data.color
    m = "*" if data.paper == "TIGRESS-GC" else "o"
    x = getattr(data, "log_{}".format(xf))
    y = getattr(data, "log_{}".format(yf))
    xerr = getattr(data, "log_{}_std".format(xf))
    yerr = getattr(data, "log_{}_std".format(yf))
    p = plt.errorbar(
        x,
        y,
        yerr=yerr,
        xerr=xerr,
        marker=m,
        ls="",
        mew=1,
        mec="k",
        lw=1,
        color=c,
        ms=ms * (1.5 if m == "*" else 1),
        label=data.paper,
    )  # ,clip_on=False)
    return p


def get_ncr_color(Z, cmap=cmr.guppy, Zmin=10 ** (-1.3), Zmax=10 ** (0.5)):
    norm = LogNorm(vmin=Zmin, vmax=Zmax)
    return cmap(norm(Z))


def add_ncr_sim(data, xf="Ptot", yf="SFR", ms=5, legend=4):
    m = "o"
    x = getattr(data, "log_{}".format(xf))
    y = getattr(data, "log_{}".format(yf))
    xerr = getattr(data, "log_{}_std".format(xf))
    yerr = getattr(data, "log_{}_std".format(yf))
    colors = get_ncr_color(data.Zdust)
    for x_, y_, xerr_, yerr_, c_ in zip(x, y, xerr, yerr, colors):
        p = plt.errorbar(
            x_,
            y_,
            yerr=yerr_,
            xerr=xerr_,
            marker=m,
            ls="",
            mew=1,
            mec="k",
            lw=1,
            color=c_,
            ms=ms * (1.5 if m == "*" else 1),
        )

    # legend like K24
    if legend:
        ec = plt.rcParams["axes.labelcolor"]
        pkwargs = dict(ls="", marker="o", markeredgecolor=ec)

        Zmodels = [0.025, 0.1, 0.3, 1, 3]
        colors = get_ncr_color(Zmodels)
        custom_lines = []
        for c in colors:
            custom_lines.append(Line2D([0], [0], color=c, **pkwargs))
        _ = plt.legend(
            custom_lines,
            Zmodels,
            title=r"$Z_{\rm d}^\prime$",
            fontsize="xx-small",
            title_fontsize="x-small",
            loc=legend,
        )

    return p


def add_PSFR_model_line(
    Wmin=2, Wmax=8, Z=None, labels=True, model="tigress-classic", **kwargs
):
    Wtmp = np.logspace(Wmin, Wmax)
    plt.plot(
        np.log10(Wtmp),
        np.log10(get_sfr(Wtmp * _kbol_cgs, Z=Z, Ytot=model) / _sfr_cgs),
        **kwargs,
    )
    if labels:
        plt.xlabel(r"$P_{\rm DE}\, [k_B{\rm\,cm^{-3}\,K}]$")
        plt.ylabel(r"$\Sigma_{\rm SFR}\, [M_\odot\,{\rm kpc^{-2}\,yr^{-1}}]$")


def add_PSFR_model_lines(Wmin=2, Wmax=8, model="classic"):
    for y, ls in zip(
        [1.0e3, f"tigress-{model}", f"tigress-{model}-decomp"], ["-", ":", "--"]
    ):
        add_PSFR_model_line(
            Wmin=Wmin,
            Wmax=Wmax,
            model=y,
            ls=ls,
            color="k",
            label=r"$\Upsilon:$" + "{}".format(y),
        )


def add_PSFR_ncr_model_lines(Wmin=2, Wmax=8):
    for Z in [0.1, 0.3, 1, 3]:
        add_PSFR_model_line(
            Wmin=Wmin,
            Wmax=Wmax,
            Z=Z,
            model="tigress-ncr-decomp",
            color=get_ncr_color(Z),
            label=r"$\Upsilon:$" + f"tigress-ncr Z={Z}",
        )


def add_yield_model_line(
    Wmin=2, Wmax=8, Z=None, labels=True, comp="th", model="tigress-classic", **kwargs
):
    Wtmp = np.logspace(Wmin, Wmax)
    plt.plot(
        np.log10(Wtmp),
        np.log10(
            get_feedback_yield_comp(Wtmp * _kbol_cgs, Z=Z, comp=comp, model=model)
        ),
        **kwargs,
    )


def add_yield_model_lines(Wmin=2, Wmax=8, comp="th"):
    for y, ls in zip(
        ["tigress-classic-decomp", "tigress-ncr-decomp"],
        [":", "--", "-"],
    ):
        add_yield_model_line(Wmin=Wmin, Wmax=Wmax, model=y, comp=comp, ls=ls, color="k")


def add_yield_ncr_model_lines(Wmin=2, Wmax=8, comp="th"):
    for Z in [0.1, 0.3, 1, 3]:
        add_yield_model_line(
            Wmin=Wmin,
            Wmax=Wmax,
            Z=Z,
            model="tigress-ncr-decomp",
            comp=comp,
            color=get_ncr_color(Z),
        )


def setup_axes(nrow=3, figsize=(12, 6.75), width_ratios=(2, 1)):
    fig = plt.figure(figsize=figsize, layout="constrained")
    spec = fig.add_gridspec(nrow, 2, width_ratios=width_ratios)

    main_ax = fig.add_subplot(spec[:, 0])
    side_axes = []
    for i in range(nrow):
        side_axes.append(fig.add_subplot(spec[i, 1]))
    return fig, main_ax, side_axes


# ========================================================================================
# loader
# ========================================================================================
def load_pretigress():
    # loading Kim & Ostriker 2014
    PRFM_KO15 = PRFM_data("KO15")
    data = PRFM_KO15
    data.SFR = np.array([1.59, 1.46, 0.91, 0.81, 0.74, 0.75, 0.56]) * 1.0e-3
    data.Pturb = np.array([5.24, 5.15, 3.4, 3.38, 2.85, 3.02, 2.18]) * 1.0e3
    data.Pth = np.array([1.78, 1.61, 1.24, 1.14, 1.04, 1.05, 0.83]) * 1.0e3
    data.dPimag = np.array([np.nan, np.nan, 0.94, 1.04, 1.03, 1.10, 1.16]) * 1.0e3
    data.oPimag = np.array([np.nan, np.nan, 0.46, 0.73, 0.87, 0.91, 1.83]) * 1.0e3

    data.SFR_std = np.array([0.38, 0.28, 0.16, 0.16, 0.15, 0.16, 0.14]) * 1.0e-3
    data.Pturb_std = np.array([3.05, 2.61, 2.42, 2.66, 1.46, 2.11, 1.12]) * 1.0e3
    data.Pth_std = np.array([0.40, 0.29, 0.18, 0.20, 0.18, 0.16, 0.16]) * 1.0e3
    data.dPimag_std = np.array([np.nan, np.nan, 0.21, 0.18, 0.17, 0.22, 0.18]) * 1.0e3
    data.oPimag_std = np.array([np.nan, np.nan, 0.15, 0.18, 0.23, 0.29, 0.45]) * 1.0e3
    data.convert_linear_log()
    data.get_Ptotal()
    data.get_yield()

    # KOK13
    PRFM_KOK13 = PRFM_data("KOK13")
    data = PRFM_KOK13
    data.log_SFR = np.array([-4.11, -3.43, -2.82, -2.18, -3.02])
    data.log_Pturb = np.array([2.52, 3.19, 3.69, 4.20, 3.60])
    data.log_Pth = np.array([2.23, 2.70, 3.23, 3.79, 3.64])

    data.log_SFR_std = np.array([0.32, 0.27, 0.12, 0.06, 0.30])
    data.log_Pturb_std = np.array([0.35, 0.35, 0.15, 0.09, 0.49])
    data.log_Pth_std = np.array([0.16, 0.15, 0.10, 0.04, 0.19])

    # cutout R-series
    for f in data.field_list:
        logf = f"log_{f}"
        logferr = f"log_{f}_std"
        if hasattr(data, logf):
            v = getattr(data, logf)
            setattr(data, logf, v[:-1])
        if hasattr(data, logferr):
            v = getattr(data, logferr)
            setattr(data, logferr, v[:-1])
    data.convert_log_linear()
    data.get_Ptotal()
    data.get_yield()

    # KKO11
    PRFM_KKO11 = PRFM_data("KKO11")
    data = PRFM_KKO11
    data.log_SFR = np.array(
        [
            -4.20,
            -3.52,
            -3.03,
            -2.74,
            -2.72,
            -2.38,
            -2.06,
            -3.85,
            -3.15,
            -2.79,
            -2.58,
            -2.24,
            -3.45,
            -2.93,
            -2.43,
            -2.31,
            -2.82,
            -2.66,
        ]
    )
    data.log_Pturb = np.array(
        [
            2.61,
            3.03,
            3.57,
            3.85,
            3.74,
            4.08,
            4.19,
            2.71,
            3.45,
            3.75,
            3.90,
            4.27,
            3.25,
            3.62,
            3.96,
            4.08,
            3.66,
            3.86,
        ]
    )
    data.log_Pth = np.array(
        [
            1.94,
            2.53,
            2.95,
            3.24,
            3.22,
            3.50,
            3.86,
            2.29,
            2.81,
            3.19,
            3.43,
            3.72,
            2.60,
            3.08,
            3.44,
            3.64,
            3.13,
            3.31,
        ]
    )

    data.log_SFR_std = np.array(
        [
            0.23,
            0.12,
            0.11,
            0.11,
            0.09,
            0.10,
            0.10,
            0.15,
            0.08,
            0.12,
            0.09,
            0.05,
            0.07,
            0.07,
            0.11,
            0.10,
            0.08,
            0.06,
        ]
    )
    data.log_Pturb_std = (
        np.array(
            [
                1.51,
                0.72,
                0.73,
                0.60,
                0.52,
                0.51,
                0.63,
                0.63,
                0.61,
                0.71,
                0.62,
                0.55,
                0.56,
                0.68,
                0.70,
                0.53,
                0.77,
                0.60,
            ]
        )
        * 0.5
    )
    data.log_Pth_std = (
        np.array(
            [
                0.35,
                0.22,
                0.16,
                0.15,
                0.08,
                0.11,
                0.07,
                0.23,
                0.21,
                0.13,
                0.13,
                0.09,
                0.23,
                0.12,
                0.12,
                0.07,
                0.20,
                0.14,
            ]
        )
        * 0.5
    )
    data.convert_log_linear()
    data.get_Ptotal()
    data.get_yield()

    return PRFM_KKO11, PRFM_KOK13, PRFM_KO15


def load_classic_data(from_table=False):
    # OK22: TIGRESS-classic
    data = PRFM_data("TIGRESS-classic")
    if from_table:
        # from table
        data.SFR = np.array([1.1, 5.37e-2, 2.67e-3, 6.21e-5, 1.17e-1, 5.41e-2, 2.16e-3])
        data.Pturb = np.array([1.26e6, 1.95e5, 5.71e3, 1.88e2, 6.41e5, 1.01e5, 4.78e3])
        data.Pth = np.array([1.13e5, 1.76e4, 5.02e3, 3.39e2, 6.60e4, 1.34e4, 2.36e3])
        data.Pimag = np.array([5.37e5, 2.22e4, 7.86e3, 1.67e2, 2.77e5, 1.92e4, 1.91e3])

        # ad hoc value
        data.SFR_std = data.SFR * 0.0
        data.Pturb_std = data.Pturb * 0.0
        data.Pth_std = data.Pth * 0.0
        data.Pimag_std = data.Pimag * 0.0

        data.convert_linear_log()
        data.get_Ptotal()
        data.get_yield()
    else:
        # my re caclulation
        data.SFR = np.array(
            [9.55e-01, 6.18e-01, 7.81e-02, 7.39e-02, 3.60e-03, 3.11e-03, 6.69e-05]
        )
        data.Pturb = np.array(
            [1.74e06, 6.16e05, 3.09e05, 1.68e05, 6.01e03, 5.15e03, 1.42e02]
        )
        data.Pth = np.array(
            [1.09e05, 7.00e04, 1.63e04, 1.35e04, 4.92e03, 2.32e03, 2.58e02]
        )
        data.oPimag = np.array(
            [1.44e05, 1.24e05, 8.77e03, 7.55e03, 4.37e03, 7.99e02, 1.88e02]
        )
        data.dPimag = np.array(
            [3.17e05, 1.65e05, 1.14e04, 1.06e04, 3.42e03, 1.13e03, 7.34e01]
        )

        data.SFR_std = np.array(
            [2.76e-01, 9.17e-02, 1.06e-02, 2.59e-02, 1.47e-03, 2.14e-03, 7.04e-05]
        )
        data.Pturb_std = np.array(
            [1.57e06, 2.47e05, 3.83e05, 1.43e05, 2.23e03, 2.77e03, 1.09e02]
        )
        data.Pth_std = np.array(
            [6.85e04, 1.54e04, 9.60e03, 3.80e03, 1.44e03, 1.30e03, 1.73e02]
        )
        data.oPimag_std = np.array(
            [1.60e05, 7.36e04, 1.21e04, 6.88e03, 1.25e03, 8.13e02, 2.08e02]
        )
        data.dPimag_std = np.array(
            [2.51e05, 7.71e04, 1.37e04, 1.09e04, 1.41e03, 7.40e02, 3.97e01]
        )
        data.convert_linear_log()
        data.get_Ptotal()
        data.get_yield()
    return data


def load_gc_data(from_table=False):
    # Galactic center from M21 and M23
    data = PRFM_data("TIGRESS-GC")
    if from_table:
        # data.field_list = ["SFR", "Pturb", "Pth", "Pimag"]
        # from paper table
        data.SFR = np.array([0.148, 0.648, 2.07, 7.94, 1.88, 6.67, 27.1, 94.8])
        data.Pturb = (
            np.array([0.137, 0.512, 1.49, 4.49, 1.44, 4.26, 12.7, 57.4]) * 1.0e6
        )
        data.Pth = np.array([0.128, 0.424, 0.966, 1.99, 1.01, 2.32, 3.82, 9.51]) * 1.0e6

        data.SFR_std = np.array([0.027, 0.099, 0.37, 0.86, 0.24, 0.98, 2.5, 7.4])
        data.Pturb_std = (
            np.array([0.0022, 0.079, 0.19, 0.78, 0.33, 0.82, 2.7, 19.7]) * 1.0e6
        )
        data.Pth_std = (
            np.array([0.0022, 0.051, 0.106, 0.32, 0.30, 0.42, 0.78, 2.06]) * 1.0e6
        )
        data.convert_linear_log()
        data.get_Ptotal()
        data.get_yield()
    else:
        # from recalculation for warm/cold
        data.field_list = ["SFR", "Pturb", "Pth", "Pimag"]
        # Read csv file
        df = pd.read_csv(os.path.join(dirpath, "prfm_ring.csv"), index_col="model")
        df.replace(0, np.nan, inplace=True)
        data.log_SFR = np.array(df["sigsfr_mean"])
        data.log_Pturb = np.array(df["ptrb_mean"])
        data.log_Pth = np.array(df["pthm_mean"])
        data.log_Pimag = np.array(df["pmag_mean"])

        data.log_SFR_std = np.array(df["sigsfr_std"])
        data.log_Pturb_std = np.array(df["ptrb_std"])
        data.log_Pth_std = np.array(df["pthm_std"])
        data.log_Pimag_std = np.array(df["pmag_std"])

        data.convert_log_linear()
        data.get_Ptotal()
        data.get_yield()
    return data


def load_ncr_data(fname="tigress_ncr_K24.nc"):
    dset = xr.open_dataset(os.path.join(dirpath, fname))
    data = PRFM_data("TIGRESS-NCR")
    data.Zdust = dset["Zdust"].sel(q="mean").data
    data.SFR = dset["sfr10"].sel(q="mean").data
    data.Pturb = dset["Pturb"].sel(q="mean").data
    data.Pth = dset["Pth"].sel(q="mean").data
    data.oPimag = dset["oPimag"].sel(q="mean").data
    data.dPimag = dset["dPimag"].sel(q="mean").data
    data.SFR_std = dset["sfr10"].sel(q="std").data
    data.Pturb_std = dset["Pturb"].sel(q="std").data
    data.Pth_std = dset["Pth"].sel(q="std").data
    data.oPimag_std = dset["oPimag"].sel(q="std").data
    data.dPimag_std = dset["dPimag"].sel(q="std").data
    data.convert_linear_log()
    data.get_Ptotal()
    data.get_yield()

    return data


def load_sim_data():
    """loading all simulation data as a dictionary

    Returns
    -------
    data : dict
        Dictionary containing all simulation data.
        PRFM class.
    """
    # pre-tigress (KKO11, KOK13, KO15)
    PRFM_KKO11, PRFM_KOK13, PRFM_KO15 = load_pretigress()
    # classic (OK22)
    PRFM_OK22 = load_classic_data()
    # GC (M21, M23)
    PRFM_GC = load_gc_data()
    # NCR (K24)
    PRFM_NCR = load_ncr_data()

    data = dict()
    colors = {
        "KKO11": "tab:gray",
        "KOK13": "tab:pink",
        "KO15": "tab:olive",
        "TIGRESS-classic": "tab:blue",
        "TIGRESS-GC": "tab:orange",
        "TIGRESS-NCR": "tab:red",
    }
    for d in [PRFM_KKO11, PRFM_KOK13, PRFM_KO15, PRFM_OK22, PRFM_GC, PRFM_NCR]:
        d.color = colors[d.paper]
        data[d.paper] = d
    return data
