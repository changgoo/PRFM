import astropy.constants as ac
import astropy.units as au
import numpy as np

Gconst_cgs = ac.G.cgs.value
pc_cgs = ac.pc.cgs.value
kbol_cgs = ac.k_B.cgs.value
surf_cgs = ac.M_sun.cgs.value/pc_cgs**2
sfr_cgs = (ac.M_sun/ac.kpc**2/au.yr).cgs.value
kms_cgs = 1.e5

sigma_eff_models = dict()
sigma_eff_models['tigress_mid'] = dict(sigma_0 = 9.8, expo=0.15, sigma_min=5)
sigma_eff_models['tigress_avg'] = dict(sigma_0 = 12, expo=0.22, sigma_min=5)

yield_models = dict()
yield_models['tigress-classic'] = dict(Y0 = 10**3.86, expo=-0.212)
yield_models['tigress-classic-decomp'] = dict(Yth0 = 10**4.45, expo_th=-0.506,
                                              Ytrb0 = 1.5*10**2.81, expo_trb=-0.06)
def get_weight_gas(Sigma_gas):
    """weight due to gas self-gravity
    """
    return 0.5*np.pi*Gconst_cgs*Sigma_gas**2

def get_weight_star(Sigma_gas,H_gas,Sigma_star,H_star):
    """weight due to stellar gravity
    """
    return np.pi*Gconst_cgs*Sigma_gas*Sigma_star*H_gas/(H_gas+H_star)

def get_weight_dm(Sigma_gas,H_gas,Omega_d,zeta_d=1/3.):
    """weight due to dark matter gravity
    """
    return zeta_d*Sigma_gas*H_gas*Omega_d**2

def get_pressure(Sigma_gas,H_gas,sigma_eff):
    """total pressure
    """
    return 0.5*Sigma_gas/H_gas*sigma_eff**2

@np.vectorize
def get_weights(Sigma_gas, Sigma_star, Omega_d, H_star,
                sigma_eff, zeta_d=1/3., method='analytic'):
    """calculate gas scale height and then each weight term

    ==========
    Parameters
    ==========
    method: str [analytic, numerical]
        how to calculate scale height

    =======
    returns
    =======
    H_gas:
        gas scale height
    W_gas:
        weight by gas
    W_star:
        weight by star
    W_dm:
        weight by dark matter
    """
    H_gas = get_scale_height(Sigma_gas,Sigma_star,Omega_d,H_star,
                             sigma_eff, zeta_d=zeta_d, method=method)

    W_gas = get_weight_gas(Sigma_gas)
    W_star = get_weight_star(Sigma_gas, H_gas, Sigma_star, H_star)
    W_dm = get_weight_dm(Sigma_gas, H_gas, Omega_d, zeta_d=zeta_d)

    return H_gas, W_gas, W_star, W_dm

@np.vectorize
def get_weight_contribution(Sigma_gas, Sigma_star, Omega_d, H_star,
                            sigma_eff, zeta_d=1/3., method='analytic'):
    """Wrapper function to calculate the ratio of each weight term to the total weight
    """
    H,wgas,wstar,wdm = get_weights(Sigma_gas,Sigma_star,Omega_d,H_star,
                                   sigma_eff, zeta_d=zeta_d, method=method)
    wtot = wgas+wstar+wdm

    return wgas/wtot,wstar/wtot,wdm/wtot

@np.vectorize
def get_scale_height_gas_only(*args,**kwargs):
    """analytic solution for gas only case
    """
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff = args
    return sigma_eff**2/(np.pi*Gconst_cgs*Sigma_gas)

@np.vectorize
def get_scale_height_star_only(*args,**kwargs):
    """analytic solution for star only case
    """
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff = args
    h = sigma_eff**2/(4*np.pi*Gconst_cgs*Sigma_star)
    h_star = H_star/h
    return h*(1+np.sqrt(1+2*h_star))

@np.vectorize
def get_scale_height_dm_only(*args,**kwargs):
    """analytic solution for dark matter only case
    """
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff = args
    zeta_d = kwargs['zeta_d']
    return sigma_eff/np.sqrt(2*zeta_d)/Omega_d

@np.vectorize
def get_scale_height_star_gas(*args,**kwargs):
    """analytic solution neglecting dark matter
    """
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff = args
    H_gas_only = get_scale_height_gas_only(*args,**kwargs)
    eta_star = H_star/H_gas_only
    s_star = Sigma_star/Sigma_gas
    h = ((1-eta_star)+np.sqrt((eta_star-1)**2+4*eta_star*(1+2.*s_star)))/(2*(1+2*s_star))
    return h*H_gas_only

@np.vectorize
def get_scale_height_star_gas_approx(*args,**kwargs):
    """analytic solution neglecting dark matter
    """
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff = args
    H_gas_only = get_scale_height_gas_only(*args,**kwargs)
    eta_star = H_star/H_gas_only
    s_star = Sigma_star/Sigma_gas
    h = 2/(1+np.sqrt(1+8*s_star/eta_star))
    return h*H_gas_only


@np.vectorize
def get_scale_height_dm_gas(*args,**kwargs):
    """analytic solution neglectic star
    """
    H_gas_only = get_scale_height_gas_only(*args,**kwargs)
    H_dm_only = get_scale_height_dm_only(*args,**kwargs)
    eta_D = H_dm_only/H_gas_only
    h = (-eta_D**2+np.sqrt(eta_D**4+4*eta_D**2))*0.5
    return h*H_gas_only

@np.vectorize
def get_sigma_eff(P, model='tigress_mid'):
    """pressure-dependent velocity dispersion model
    """
    sigma_eff_model = sigma_eff_models[model]
    sigma_0 = sigma_eff_model['sigma_0']
    expo = sigma_eff_model['expo']
    sigma_min = sigma_eff_model['sigma_min']

    veld = sigma_0*(1.e-4*P/kbol_cgs)**expo

    return np.clip(veld,sigma_min,None)*1.e5

@np.vectorize
def get_feedback_yield(P, model='tigress-classic'):
    """total feedback yield as a function of weight

    ===========
    Paramerters
    ===========
    P : float
        midplane pressure or weight
    model : str
        tigress-classic for OK22
        tigress-NCR for K23
    """
    yield_model = yield_models[model]
    if 'Y0' in yield_model:
        Y0 = yield_model['Y0']
        slope = yield_model['expo']
        Ytot = Y0*(P/kbol_cgs)**slope
    elif 'Yth0' in yield_model:
        Yth0 = yield_model['Yth0']
        expo_th = yield_model['expo_th']
        Ytrb0 = yield_model['Ytrb0']
        expo_trb = yield_model['expo_trb']
        Ytot = Yth0*(P/kbol_cgs)**expo_th + Ytrb0*(P/kbol_cgs)**expo_trb
    return Ytot

@np.vectorize
def get_feedback_yield_comp(P, comp='th', model='tigress-classic-decomp'):
    """total feedback yield as a function of weight

    ===========
    Paramerters
    ===========
    P : float
        midplane pressure or weight
    model : str
        tigress-classic for OK22
        tigress-NCR for K23
    """
    yield_model = yield_models[model]
    Yth0 = yield_model['Y{}0'.format(comp)]
    expo_th = yield_model['expo_{}'.format(comp)]
    return Yth0*(P/kbol_cgs)**expo_th

def get_sfr(P, Ytot='tigress-classic'):
    Ytot = get_feedback_yield(P,model=Ytot) if type(Ytot) == str else Ytot
    return P/(Ytot*1.e5)

def get_scale_height(Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff,
                     zeta_d=1/3., wgas=1, wstar=1, wdm=1, method='analytic'):
    """wrapper function to calculate gas scale height either numerically or analytically

    ===========
    Parameters
    ===========
    Sigma_gas: float
        gas surface density
    Sigma_star: float
        stellar surface density
    Omega_d: float
        galactic rotation speed
    H_star: float
        stellar scale height
    sigma_eff: str or float
        velocity dispersion

    zeta_d: float [1/3.]
        parameter for dm and gas distribution
    wgas: int [0 or 1]
        toggle weight term from gas
    wstar: int [0 or 1]
        toggle weight term from stars
    wdm: int [0 or 1]
        toggle weight term from dark matter
    method: str [analytic or numerical]
        calculate scale height either analytically or numerically
    """
    args = (Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff)
    kwargs = dict(zeta_d=zeta_d, wgas=wgas, wstar=wstar, wdm=wdm)
    w=(wgas<<2)+(wstar<<1)+(wdm<<0)

    if method == 'numerical':
        return get_scale_height_numerical(*args, **kwargs)
    else:
        if w == 7: # 111
            return get_scale_height_analytic(*args, zeta_d=zeta_d)
        elif w == 6: # 110
            return get_scale_height_star_gas(*args, zeta_d=zeta_d)
        elif w == 5: # 101
            return get_scale_height_dm_gas(*args, zeta_d=zeta_d)
        elif w == 4: # 100
            return get_scale_height_gas_only(*args, zeta_d=zeta_d)
        elif w == 3: # 011
            # no analytic solution provided
            return
        elif w == 2: # 010
            return get_scale_height_star_only(*args, zeta_d=zeta_d)
        elif w == 1: # 001
            return get_scale_height_dm_only(*args, zeta_d=zeta_d)
        else:
            # no such gas is expected
            raise("at least one term must be considered")

@np.vectorize
def get_scale_height_analytic(Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff,
                              zeta_d=1/3.):
    """Analytic solution of the cubic equation for the vertical dynamical equilibirum
    including all three weight terms.

    All inputs must be in c.g.s. units
    """
    args = (Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff)
    kwargs = dict(zeta_d=zeta_d)
    H_gas_only = get_scale_height_gas_only(*args,**kwargs)
    H_dm_only = get_scale_height_dm_only(*args,**kwargs)

    s_star = Sigma_star/Sigma_gas
    eta_star = H_star/H_gas_only
    eta_dm_sq = (H_dm_only/H_gas_only)**2

    a0 = -eta_star*eta_dm_sq
    a1 = (eta_star-1)*eta_dm_sq
    a2 = (1+2*s_star)*eta_dm_sq + eta_star

    Q = (3*a1-a2**2)/9.
    R = (9*a2*a1-27*a0-2*a2**3)/54.

    theta = np.arccos(R/np.sqrt(-Q**3))
    h = 2*np.sqrt(-Q)*np.cos(theta/3)-a2/3.

    return h*H_gas_only

from scipy.optimize import root, brentq
@np.vectorize
def get_scale_height_numerical(Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff,
                               zeta_d=1/3., wgas=1, wstar=1, wdm=1):
    """Numerical solution of the vertical dynamical equilibrium equation
    including all weight terms.

    All inputs must be in c.g.s. units
    """
    fun = lambda x: (get_pressure(Sigma_gas, x, sigma_eff)
                    - wgas*get_weight_gas(Sigma_gas)
                    - wstar*get_weight_star(Sigma_gas,x,Sigma_star,H_star)
                    - wdm*get_weight_dm(Sigma_gas, x, Omega_d, zeta_d=zeta_d))
    # soln=root(fun,H_gas_init)
    # return soln.x
    return brentq(fun,1.e-3*pc_cgs,1.e6*pc_cgs)

def get_self_consistent_solution(Sigma_gas, Sigma_star, Omega_d, H_star,
                                 sigma_eff='tigress_mid',
                                 zeta_d=1/3.,
                                 method='analytic',
                                 tol = 1.e-5,
                                 niter = 16
                                ):
    # make initial guess
    if type(sigma_eff) != str:
        sigma0 = sigma_eff
    else:
        if not (sigma_eff in sigma_eff_models):
            raise IndexError("sigma_eff models: ",
                             list(sigma_eff_models.keys()))
        sigma0 = 15.e5
    args = (Sigma_gas, Sigma_star, Omega_d, H_star, sigma0)
    kwargs = dict(zeta_d=zeta_d, method=method)
    Hprev, wgas, wstar, wdm = get_weights(*args, **kwargs)
    wtot_prev = wgas + wstar + wdm

    # return if velocity diseprsion is a constant
    if type(sigma_eff) != str: return wtot_prev, Hprev, sigma0*np.ones_like(Hprev)

    # iterative solve
    for i in range(niter):
        new_sigma_eff = get_sigma_eff(wtot_prev, model=sigma_eff)
        args = (Sigma_gas, Sigma_star, Omega_d, H_star, new_sigma_eff)
        Hnext, wgas, wstar, wdm = get_weights(*args, **kwargs)
        wtot_next = wgas + wstar + wdm
        # L1_norm = np.sum(np.abs(wtot_next/wtot_prev - 1))
        L1_norm = np.sum(np.abs(Hnext/Hprev - 1))
        if (L1_norm < tol): break
        wtot_prev = np.copy(wtot_next)
        Hprev = np.copy(Hnext)

    return wtot_next, Hnext, get_sigma_eff(wtot_next, model=sigma_eff)
