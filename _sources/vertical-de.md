# Vertical Dynamical Equilibrium

## Basic fomulation

The vertical dynamical equilibrium in a disk consisting of gas, star, and dark matter can be written as

$$
  P_{\rm tot} = \frac{\Sigma_{\rm gas}}{2 H_{\rm gas}} \sigma_{\rm eff}^2
  = {\cal W}_\mathrm{g} + {\cal W}_* + {\cal W}_\mathrm{d}.
$$

where

$$
{\cal W}_\mathrm{g}= \frac{\pi G \Sigma_{\rm gas}^2}{2}
$$

$$
{\cal W}_* \approx \pi G \Sigma_{\rm gas} \Sigma_* \frac{H_{\rm gas}}{H_{\rm gas} + H_*},
$$

$$
{\cal W}_\mathrm{d}  = \zeta_\mathrm{d} \Sigma_{\rm gas} \Omega_d^2 H_{\rm gas}
$$

Given $\Sigma_{\rm gas}$, $\Sigma_*$, $H_*$, and $\Omega_d$, we solve for $H_{\rm gas}$ with an assumption of $\sigma_{\rm eff}$.
We use a constant $\sigma_{\rm eff}$ or a model calibarted by numerical simulations.

## Constant velocity dispersion
For a constant $\sigma_{\rm eff}$, the equation becomes a cubic to $H_{\rm gas}$. Using

$$
H_{\rm gas-only} = \frac{\sigma_{\rm eff}^2}{\pi G \Sigma_{\rm gas}}
$$

and

$$
H_{\rm dm-only} = \frac{\sigma_{\rm eff}}{(2\zeta_d)^{1/2}\Omega_d},
$$

the dimensionless form of the equation can be written as

$$
h^3 + [(1+2 s_*)\eta_d^2+\eta_*]h^2 + (\eta_*-1)\eta_d^2h -\eta_*\eta_d^2 = 0
$$

with
$s_* = \Sigma_*/\Sigma_{\rm gas}$, $\eta_* = H_*/H_{\rm gas-only}$, and $\eta_d = H_{\rm dm-only}/H_{\rm gas-only}$.
The only postive real solution of the cubic equation is

$$
H_{\rm gas} = 2\sqrt{-Q}\cos\left(\frac{\theta}{3}\right), \theta = \cos^{-1}\left(\frac{R}{\sqrt{-Q}}\right)
$$
where
$Q = (3a_1-a_2^2)/9$ and $R=(9a_2a_1-27a_0-2a_2^3)/54$ with $a_i$ for the coefficients of $h^i$ for $i=0$, 1, 2.

## Pressure-dependent velocity dispersion
We obtained $\sigma_{\rm eff}$ as a function of the total pressure calibrated in the TIGRESS-classic suite {cite}`2022ApJ...936..137O`.

From the effective equation of state determined by the midplane pressure and density directly, we obtain

$$
\sigma_{\rm eff} = 9.8{\rm km/s} P_4^{0.15}
$$

where $P_4\equiv(P_{\rm tot}/k_B\,10^4{\rm\,cm^{-3}\,K})$. Note that under the assumption of the vertical dynamical equilibrium,
$P_{\rm tot}=W_{\rm tot}$.

From the mass-weighted mean velocity dispersion, we have

$$
\sigma_{\rm eff} = 12{\rm km/s} P_4^{0.22}
$$

Either case, we can first solve the cubic equation for a contant $\sigma_{\rm eff}$ to get $H_{\rm gas}$ and $P_{\rm tot}=W_{\rm tot}$. We update $\sigma_{\rm eff}$ with new $P_{\rm tot}$ and then solve for $H_{\rm gas}$ and $P_{\rm tot}$ again. We iteratively solve until $H_{\rm gas}$ is converged.

## References

The most general form of the vertical dynamical equilibrium used here is derived in {cite}`Hassan_et_al`.
See also {cite}`2022ApJ...936..137O`, {cite}`2015ApJ...815...67K`, {cite}`2013ApJ...776....1K`,
{cite}`2011ApJ...743...25K`, {cite}`2011ApJ...731...41O`, {cite}`2010ApJ...721..975O`.