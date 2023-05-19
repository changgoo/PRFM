from argparse import ArgumentError
import astropy.constants as ac
import astropy.units as au
import numpy as np
from scipy.optimize import brentq

# setting up some convenient conversion factors
_Gconst_cgs = ac.G.cgs.value
_pc_cgs = ac.pc.cgs.value
_kbol_cgs = ac.k_B.cgs.value

_surf_cgs = (ac.M_sun/ac.pc**2).cgs.value
_rho_cgs = (ac.M_sun/ac.pc**3).cgs.value
_sfr_cgs = (ac.M_sun/ac.kpc**2/au.yr).cgs.value
_kms_cgs = (1*au.km/au.s).cgs.value

_sigma_eff_models = dict()
_sigma_eff_models['tigress_mid'] = dict(sigma_0 = 9.8, expo=0.15, sigma_min=5)
_sigma_eff_models['tigress_avg'] = dict(sigma_0 = 12, expo=0.22, sigma_min=5)

_yield_models = dict()
_yield_models['tigress-classic'] = dict(Y0 = 10**3.86, expo=-0.212)
_yield_models['tigress-classic-decomp'] = dict(Yth0 = 10**4.45, expo_th=-0.506,
                                              Ytrb0 = 1.5*10**2.81, expo_trb=-0.06)
_yield_models['tigress-ncr-decomp'] = dict(Yth0 = 10**4.7, expo_th=-0.5,
                                          Ytrb0 = 10**3.3, expo_trb=-0.1)
_yield_models['tigress-ncr-decomp-Z01'] = dict(Yth0 = 10**5.2, expo_th=-0.5,
                                          Ytrb0 = 10**3.3, expo_trb=-0.1)

def get_weight_gas(Sigma_gas):
    """weight due to gas self-gravity
    """
    return 0.5*np.pi*_Gconst_cgs*Sigma_gas**2

def get_weight_star(Sigma_gas,H_gas,Sigma_star,H_star):
    """weight due to stellar gravity
    """
    return np.pi*_Gconst_cgs*Sigma_gas*Sigma_star*H_gas/(H_gas+H_star)

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
    H,wgas,wstar,wdm = get_weights(Sigma_gas,Sigma_star,Omega_d,H_star,
                                   sigma_eff, zeta_d=zeta_d, method=method)
    wtot = wgas+wstar+wdm

    return wgas/wtot,wstar/wtot,wdm/wtot

@np.vectorize
def get_scale_height_gas_only(*args,**kwargs):
    """analytic solution for gas only case
    """
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff = args
    return sigma_eff**2/(np.pi*_Gconst_cgs*Sigma_gas)

@np.vectorize
def get_scale_height_star_only(*args,**kwargs):
    """analytic solution for star only case
    """
    Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff = args
    h = sigma_eff**2/(4*np.pi*_Gconst_cgs*Sigma_star)
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
    sigma_eff_model = _sigma_eff_models[model]
    sigma_0 = sigma_eff_model['sigma_0']
    expo = sigma_eff_model['expo']
    sigma_min = sigma_eff_model['sigma_min']

    veld = sigma_0*(1.e-4*P/_kbol_cgs)**expo

    return np.clip(veld,sigma_min,None)*1.e5

@np.vectorize
def get_feedback_yield(P, model='tigress-classic'):
    """total feedback yield as a function of weight

    Paramerters
    -----------
    P : float
        midplane pressure or weight
    model : str
        tigress-classic for OK22
        tigress-NCR for K23

    Returns
    -------
    Ytot : float
        feedback yield in km/s
    """
    yield_model = _yield_models[model]
    if 'Y0' in yield_model:
        Y0 = yield_model['Y0']
        slope = yield_model['expo']
        Ytot = Y0*(P/_kbol_cgs)**slope
    elif 'Yth0' in yield_model:
        Yth0 = yield_model['Yth0']
        expo_th = yield_model['expo_th']
        Ytrb0 = yield_model['Ytrb0']
        expo_trb = yield_model['expo_trb']
        Ytot = Yth0*(P/_kbol_cgs)**expo_th + Ytrb0*(P/_kbol_cgs)**expo_trb
    return Ytot

@np.vectorize
def get_feedback_yield_comp(P, comp='th', model='tigress-classic-decomp'):
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
    Y0 = yield_model['Y{}0'.format(comp)]
    expo = yield_model['expo_{}'.format(comp)]
    return Y0*(P/_kbol_cgs)**expo

def get_sfr(P, Ytot='tigress-classic'):
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
    Ytot = get_feedback_yield(P,model=Ytot) if type(Ytot) == str else Ytot
    return P/(Ytot*1.e5)

def get_scale_height(Sigma_gas, Sigma_star, Omega_d, H_star, sigma_eff,
                     zeta_d=1/3., method='analytic', wgas=1, wstar=1, wdm=1):
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
    return brentq(fun,1.e-3*_pc_cgs,1.e6*_pc_cgs)

def get_self_consistent_solution(Sigma_gas, Sigma_star, Omega_d, H_star,
                                 sigma_eff='tigress_mid',
                                 zeta_d=1/3.,
                                 method='analytic',
                                 tol = 1.e-5,
                                 niter = 16
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
    if type(sigma_eff) != str:
        sigma0 = sigma_eff
    else:
        if not (sigma_eff in _sigma_eff_models):
            raise IndexError("sigma_eff models: ",
                             list(_sigma_eff_models.keys()))
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

    def __init__(self,
                 Sigma_gas=None,
                 Sigma_star=None,
                 H_star=None,
                 rho_star=None,
                 Omega_d=None,
                 rho_dm=None,
                 sigma_eff='tigress-mid',
                 astro_units=True,
                 zeta_d=1/3., method='analytic'):
        # initialize method attribute
        self._method=method

        # set model name based on sigma_eff
        self.name = sigma_eff if type(sigma_eff) == str else 'constant'

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
        else:
            if H_star is None:
                if rho_star is None:
                    # this is limit for H_star >> H_gas
                    self._stellar_disk = "thick"
                    self._Sigma_star = Sigma_star
                elif Sigma_star is None:
                    # this is limit for H_star << H_gas
                    self._stellar_disk = "thin"
                    self._rho_star = rho_star
                else:
                    # both given
                    self._stellar_disk = "general"
                    self._H_star = Sigma_star/(2*rho_star)
            else:
                if rho_star is None:
                    # Sigma_star is given
                    rho_star = Sigma_star/(2*H_star)
                elif Sigma_star is None:
                    # rho_star is given
                    Sigma_star = rho_star*2*H_star
                else:
                    # both given
                    raise ArgumentError("cannot provide all three: Sigma_star, rho_star, H_star")

                self._stellar_disk = "general"
                self._Sigma_star = Sigma_star
                self._rho_star = rho_star
                self._H_star = H_star
        # unit conversion done before DM part
        if astro_units: self._astro_to_cgs()

        # setting up dark matter
        if (Omega_d is None) and (rho_dm is None):
            self._Omega_d = 0
            self._rho_dm = 0
        else:
            # this conversion assumes flat rotation
            _fourpiG = 4*np.pi*ac.G.cgs.value
            if rho_dm is None:
                # Omega_d is given
                if astro_units:
                    Omega_d = Omega_d*self.units['Omega_d'].cgs.value
                rho_dm = Omega_d**2/_fourpiG
            elif Omega_d is None:
                # rho_dm is given
                if astro_units:
                    rho_dm = rho_dm*self.units['rho_dm'].cgs.value
                Omega_d = np.sqrt(_fourpiG*rho_dm)
            else:
                # both given
                raise ArgumentError("cannot provide both Omega_d and rho_dm")

            self._Omega_d = Omega_d
            self._rho_dm = rho_dm
            self._zeta_d = zeta_d

        if type(sigma_eff) != str:
            if astro_units:
                u = self.units['sigma_eff'].cgs.value
                sigma_eff = sigma_eff*u
            self._sigma_eff = sigma_eff
        else:
            self._sigma_eff_model = sigma_eff

        # set arguments that will be passed to the functions
        self._args = (self._Sigma_gas,self._Sigma_star,self._Omega_d,
                      self._H_star,sigma_eff)

        # store parameters in astro-friendly units
        self._cgs_to_astro()

    def reset_arg_list(self,sigma_eff):
        # reset arguments with new velocity dispersion
        self._args = (self._Sigma_gas,self._Sigma_star,self._Omega_d,
                      self._H_star,sigma_eff)

    def _astro_to_cgs(self):
        for var in ['Sigma_gas','Sigma_star','H_star','rho_star',
                    'rho_dm','Omega_d','sigma_eff']:
            if hasattr(self,'_'+var):
                u_cgs = self.units[var].cgs.value
                v = getattr(self,'_'+var)*u_cgs
                setattr(self,'_'+var,v)

    def _cgs_to_astro(self):
        for var in ['Sigma_gas','Sigma_star','H_star','rho_star',
                    'rho_dm','Omega_d','sigma_eff']:
            if hasattr(self,'_'+var):
                u_cgs = self.units[var].cgs.value
                v = getattr(self,'_'+var)/u_cgs
                setattr(self,var,v)


    def _set_units(self):
        units = dict()
        units['Sigma_gas'] = units['Sigma_star'] = ac.M_sun/ac.pc**2
        units['rho_star'] = units['rho_dm'] = ac.M_sun/ac.pc**3
        units['Omega_d'] = au.km/au.s/ac.kpc
        units['H_star'] = units['H_gas'] = ac.pc
        units['Wtot'] = ac.k_B/au.cm**3*au.K
        units['Wgas'] = units['Wstar'] = units['Wdm'] = units['Wtot']
        units['Ytot'] = 1.*au.km/au.s
        units['sigma_eff'] = 1.*au.km/au.s
        units['Sigma_SFR'] = ac.M_sun/ac.kpc**2/au.yr

        return units

    def __repr__(self) -> str:
        args = "PRFM calculator is prepared for\n"
        for var in ['Sigma_gas','Sigma_star','rho_star',
                    'Omega_d','rho_dm','sigma_eff']:
            if not hasattr(self,'_'+var): continue
            v = np.atleast_1d(getattr(self,var))
            args += "  {}[0]: {}, N={}\n".format(var,v[0],len(v))
        if hasattr(self,'_sigma_eff_model'):
            args += "  sigma_eff_model: {}\n".format(self._sigma_eff_model)
        return args

    def get_scale_height(self,method=None,wgas=1,wstar=1,wdm=1):
        """Wrapper function to get the gas scale height

        Parameters
        ----------
        method : str
            override method if provided
        """

        if method is None: method=self._method

        H_gas = get_scale_height(*self._args,zeta_d=self._zeta_d,
                                 method=method,
                                 wgas=wgas,wstar=wstar,wdm=wdm)
        return H_gas/self.units['H_gas'].cgs.value

    def get_weight_contribution(self):
        """Wrapper function to get weight contribution
        """
        if not hasattr(self,'Wtot'): self.calc_weights()
        return self.Wgas/self.Wtot,self.Wstar/self.Wtot,self.Wdm/self.Wtot

    def calc_weights(self,method=None):
        """Wrapper function to calculate scale height and all weights

        Parameters
        ----------
        method : str
            override method if provided
        """

        if method is None: method=self._method

        results = get_weights(*self._args,zeta_d=self._zeta_d,method=method)

        for var,v in zip(['H_gas','Wgas','Wstar','Wdm'],results):
            u_cgs = self.units[var].cgs.value
            setattr(self,'_'+var,v) # results in cgs
            setattr(self,var,v/u_cgs) # results in astro units

        # total weight
        Wtot = results[1]+results[2]+results[3]
        u_cgs = self.units['Wtot'].cgs.value
        self._Wtot = Wtot
        self.Wtot = Wtot/u_cgs

        if hasattr(self,'_sigma_eff_model'):
            self._sigma_eff = get_sigma_eff(Wtot, model=self._sigma_eff_model)
            self.sigma_eff = self._sigma_eff/self.units['sigma_eff'].cgs.value

    def calc_self_consistent_solution(self,method=None,
                                      niter=16,tol=1.e-6):
        """Wrapper function to calculate self-consistent solutions

        Update H_gas, W_tot (and all components), sigma_eff
        """
        if method is None: method=self._method

        # for model sigma, need to give an initial guess
        if hasattr(self,'_sigma_eff_model'):
            self.reset_arg_list(15.e5)
        self.calc_weights(method=method)

        # that's it for constant sigma_eff
        # iterative solve for model sigma
        if hasattr(self,'_sigma_eff_model'):
            # iterative solve
            for i in range(niter):
                Hprev = np.copy(self._H_gas)
                self.reset_arg_list(self._sigma_eff)
                self.calc_weights(method=method)
                L1_norm = np.sum(np.abs(self._H_gas/Hprev - 1))
                if (L1_norm < tol): break

    def calc_sfr(self,Ytot=1.e3):
        """Wrapper function to calculate SFR surface density

        Update Sigma_SFR (astro units) and _Sigma_SFR (cgs)
        """
        if not hasattr(self,'_Wtot'):
            self.calc_self_consistent_solution()

        # results in c.g.s.
        Sigma_SFR = get_sfr(self._Wtot,Ytot=Ytot)
        u_cgs = self.units['Sigma_SFR'].cgs.value
        self._Sigma_SFR = Sigma_SFR
        self.Sigma_SFR = Sigma_SFR/u_cgs

    def check_solutions(self):
        self.calc_self_consistent_solution(method='analytic')
        sola = np.copy(self.Wtot)
        self.calc_self_consistent_solution(method='numerical')
        soln = np.copy(self.Wtot)

        L1 = np.mean(np.abs(sola-soln))
        print("difference between analytic and numerical solutions: {}".format(L1))
        return sola,soln

