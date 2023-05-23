# Feedback Yields

In the ISM, turbulent and thermal pressure support is quickly lost due to turbulence dissipation and cooling.
In order to maintain the vertical dynamical equilibrium, there must exist efficient recovery of pressure support.
Stellar feedback (esp. from massive young stars) is the major source of energy.
Turbulence is mainly driven by SN feedback, and heating is dominated by the photoelectric effect of FUV radiation on grains.
Assuming balance between gain and loss processes, one can define the feedback yield for each pressure component as a ratio between pressure and SFR surface density:

$$
  \Upsilon_{\rm th, turb, mag} = \frac{P_{\rm th, turb, mag}}{\Sigma_{\rm SFR}}
$$

Taking the sum of pressures,

$$
  P_{\rm tot} = P_{\rm th} + P_{\rm turb} + P_{\rm mag} =  (\Upsilon_{\rm th} + \Upsilon_{\rm turb} + \Upsilon_{\rm mag})\Sigma_{\rm SFR}.
$$

We can then obtain the SFR surface density from the total weight $\cal W_{\rm tot}$ and the calibrated values of $\Upsilon$, under the assumption
that the total pressure $P_{\rm tot}$ is in balance with the total weight $\cal W_{\rm tot}$.

An analytic form of the weight in vertical dynamical equilibrium is

$$
P_{\rm DE} = \frac{\pi G}{2}\Sigma_g^2 + [2 G (\rho_{sd})]^{1/2} \sigma_\mathrm{eff}
$$

for $\Sigma_g$ the gas surface density, $\rho_{sd}$ the midplane star+dark matter midplane density, and  $\sigma_\mathrm{eff}$ the effective vertical
velocity dispersion.

The numerically calibrated total feedback yield from the TIGRESS-classic suite ({cite}`2022ApJ...936..137O`; see Eqs (25)) is either

$$
  \Upsilon_{\rm tot} = 1030{\rm km/s} \left(\frac{P_{\rm DE}/k_B}{10^4 {\rm cm^{-3}\,K}}\right)^{-0.212}
$$

or the sum of thermal and kinetic+magnetic

$$
  \Upsilon_{\rm th} = 267{\rm km/s} \left(\frac{P_{\rm DE}/k_B}{10^4 {\rm cm^{-3}\,K}}\right)^{-0.506}
$$

$$
  \Upsilon_{\rm kin+mag} = 1.5*370{\rm km/s} \left(\frac{P_{\rm DE}/k_B}{10^4 {\rm cm^{-3}\,K}}\right)^{-0.060}
$$

where we assume that the magnetic contribution is 50\% of the kinetic.
