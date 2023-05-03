import astropy.constants as ac
import astropy.units as au
import numpy as np

Gconst_cgs = ac.G.cgs.value
pc_cgs = ac.pc.cgs.value
kbol_cgs = ac.k_B.cgs.value
surf_cgs = ac.M_sun.cgs.value/pc_cgs**2
kms_cgs = 1.e5

sigma_eff_models = dict()
sigma_eff_models['tigress_mid'] = dict(sigma_0 = 9.8, expo=0.15, sigma_min=10)
sigma_eff_models['tigress_avg'] = dict(sigma_0 = 12, expo=0.22, sigma_min=10)
def get_weight_gas(Sigma_gas):
    return 0.5*np.pi*Gconst_cgs*Sigma_gas**2

def get_weight_star(Sigma_gas,H_gas,Sigma_star,H_star):
    return np.pi*Gconst_cgs*Sigma_gas*Sigma_star*H_gas/(H_gas+H_star)

def get_weight_dm(Sigma_gas,H_gas,Omega_d,zeta_d=1/3.):
    return zeta_d*Sigma_gas*H_gas*Omega_d**2

def get_pressure(Sigma_gas,H_gas,sigma_eff):
    return 0.5*Sigma_gas/H_gas*sigma_eff**2

from scipy.optimize import root, brentq
@np.vectorize
def get_scale_height_numerical(Sigma_gas, Sigma_star, Omega_d, H_star,
                               sigma_eff=15.e5, zeta_d=1/3.,
                               wgas=1, wstar=1, wdm=1):
    fun = lambda x: (get_pressure(Sigma_gas, x, sigma_eff)
                    - wgas*get_weight_gas(Sigma_gas)
                    - wstar*get_weight_star(Sigma_gas,x,Sigma_star,H_star)
                    - wdm*get_weight_dm(Sigma_gas, x, Omega_d, zeta_d=zeta_d))
    # soln=root(fun,H_gas_init)
    # return soln.x
    return brentq(fun,pc_cgs,1.e6*pc_cgs)

@np.vectorize
def get_weights(Sigma_gas, Sigma_star, Omega_d, H_star,
                sigma_eff=15.e5, zeta_d=1/3.):
    H_gas = get_scale_height(Sigma_gas,Sigma_star,Omega_d,H_star,
                             sigma_eff=sigma_eff, zeta_d=zeta_d)

    wgas = get_weight_gas(Sigma_gas)
    wstar = get_weight_star(Sigma_gas, H_gas, Sigma_star, H_star)
    wdm = get_weight_dm(Sigma_gas, H_gas, Omega_d, zeta_d=zeta_d)

    return H_gas, wgas,wstar,wdm

@np.vectorize
def get_weight_contribution(Sigma_gas, Sigma_star, Omega_d, H_star,
                            sigma_eff=15.e5, zeta_d=1/3.):
    H, wgas,wstar,wdm = get_weights(Sigma_gas,Sigma_star,Omega_d,H_star,
                                    sigma_eff=sigma_eff, zeta_d=zeta_d)
    wtot = wgas+wstar+wdm

    return wgas/wtot,wstar/wtot,wdm/wtot

@np.vectorize
def get_scale_height_gas_only(*args,**kwargs):
    Sigma_gas, Sigma_star, Omega_d, H_star = args
    sigma_eff = kwargs['sigma_eff']
    return sigma_eff**2/(np.pi*Gconst_cgs*Sigma_gas)

@np.vectorize
def get_scale_height_star_only(*args,**kwargs):
    Sigma_gas, Sigma_star, Omega_d, H_star = args
    sigma_eff = kwargs['sigma_eff']

    h = sigma_eff**2/(4*np.pi*Gconst_cgs*Sigma_star)
    h_star = H_star/h
    return h*(1+np.sqrt(1+2*h_star))

@np.vectorize
def get_scale_height_dm_only(*args,**kwargs):
    Sigma_gas, Sigma_star, Omega_d, H_star = args
    sigma_eff = kwargs['sigma_eff']
    zeta_d = kwargs['zeta_d']
    return sigma_eff/np.sqrt(2*zeta_d)/Omega_d

@np.vectorize
def get_scale_height_star_gas(*args,**kwargs):
    Sigma_gas, Sigma_star, Omega_d, H_star = args
    H_gas_only = get_scale_height_gas_only(*args,**kwargs)
    eta_star = H_star/H_gas_only
    s_star = Sigma_star/Sigma_gas
    h = ((1-eta_star)+np.sqrt((eta_star-1)**2+4*eta_star*(1+2.*s_star)))/(2*(1+2*s_star))
    return h*H_gas_only

@np.vectorize
def get_scale_height_dm_gas(*args,**kwargs):
    H_gas_only = get_scale_height_gas_only(*args,**kwargs)
    H_dm_only = get_scale_height_dm_only(*args,**kwargs)
    eta_D = H_dm_only/H_gas_only
    h = (-eta_D**2+np.sqrt(eta_D**4+4*eta_D**2))*0.5
    return h*H_gas_only

@np.vectorize
def get_sigma_eff(P, model='tigress_mid'):
    sigma_eff_model = sigma_eff_models[model]
    sigma_0 = sigma_eff_model['sigma_0']
    expo = sigma_eff_model['expo']
    sigma_min = sigma_eff_model['sigma_min']

    veld = sigma_0*(1.e-4*P/kbol_cgs)**expo

    return np.clip(veld,sigma_min,None)*1.e5

@np.vectorize
def get_scale_height(*args,**kwargs):
    """ Function to calculate gas scale height in vertical dynamical equilibrium

    ===========
    Parameters
    ===========
    Sigma_gas: float

    Sigma_star: float

    Omega_d: float

    H_star: float

    sigma_eff: str or float

    """
    Sigma_gas, Sigma_star, Omega_d, H_star = args

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

def get_self_consistent_solution(*args, zeta_d=1/3.,
                                 sigma_eff='tigress_mid',
                                 L1_norm = np.inf,
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
    H, wgas, wstar, wdm = get_weights(*args, zeta_d=zeta_d,
                                      sigma_eff=sigma0)
    wtot_prev = wgas + wstar + wdm

    # return if velocity diseprsion is a constant
    if type(sigma_eff) != str: return wtot_prev, H, sigma0*np.ones_like(H)

    # iterative solve
    new_sigma_eff = get_sigma_eff(wtot_prev, model=sigma_eff)


    for i in range(niter):
        H, wgas, wstar, wdm = get_weights(*args, zeta_d=zeta_d,
                                          sigma_eff=new_sigma_eff)
        wtot_next = wgas + wstar + wdm
        L1_norm = np.sum(np.abs(wtot_next/wtot_prev - 1))
        new_sigma_eff = get_sigma_eff(wtot_prev, model=sigma_eff)
        if (L1_norm < tol): break

        wtot_prev = np.copy(wtot_next)

    return wtot_next, H, new_sigma_eff

