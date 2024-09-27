# Effective Equation of State

## TIGRESS-classic calibration

Using the TIGRESS-classic suite {cite:p}`2022ApJ...936..137O`, we  obtained a calibration for $\sigma_{\rm eff}$ as a function of the total pressure.

From the effective equation of state, in which a power-law fit is obtained between measured midplane pressure and density, one can solve to obtain

$$
\sigma_{\rm eff} = 9.8{\rm km/s} P_4^{0.15}
$$

where $P_4\equiv(P_{\rm tot}/k_B\,10^4{\rm\,cm^{-3}\,K})$. Note that under the assumption of the vertical dynamical equilibrium,
$P_{\rm tot}=W_{\rm tot}$.

Alternatively, from the direct power-law fit for the mass-weighted mean velocity dispersion in terms of pressure, we have

$$
\sigma_{\rm eff} = 12{\rm km/s} P_4^{0.22}
$$

<!-- In either case, we can first solve the cubic equation for a constant $\sigma_{\rm eff}$ to get $H_{\rm gas}$ and $P_{\rm tot}=W_{\rm tot}$. For the non-constant case, we then update $\sigma_{\rm eff}$ based on the new $P_{\rm tot}$ and then solve for $H_{\rm gas}$ and $P_{\rm tot}$ again. We iteratively solve until $H_{\rm gas}$ is converged. -->

## TIGRESS-NCR calibration

Using the TIGRESS-classic suite {cite:p}`2024ApJ...972...67K`, we  obtained a calibration for $\sigma_{\rm eff}$ as a function of the total pressure and gas metallicity, but the latter dependence is found to be very weak.

From the effective equation of state, in which a power-law fit is obtained between measured midplane pressure and density, one can solve to obtain

$$
\sigma_{\rm eff} = 8.9{\rm km/s} P_4^{0.08} {Z_g^\prime}^{-0.005}
$$

where $P_4\equiv(P_{\rm tot}/k_B\,10^4{\rm\,cm^{-3}\,K})$. Note that under the assumption of the vertical dynamical equilibrium,
$P_{\rm tot}=W_{\rm tot}$.

Alternatively, from the direct power-law fit for the mass-weighted mean velocity dispersion in terms of pressure, we have

$$
\sigma_{\rm eff} = 11.7{\rm km/s} P_4^{0.12}  {Z_g^\prime}^{-0.03}
$$

<!-- In either case, we can first solve the cubic equation for a constant $\sigma_{\rm eff}$ to get $H_{\rm gas}$ and $P_{\rm tot}=W_{\rm tot}$. For the non-constant case, we then update $\sigma_{\rm eff}$ based on the new $P_{\rm tot}$ and then solve for $H_{\rm gas}$ and $P_{\rm tot}$ again. We iteratively solve until $H_{\rm gas}$ is converged. -->
