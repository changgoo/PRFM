# Vertical Dynamical Equilibrium

## Basic formulation

The vertical dynamical equilibrium in a disk consisting of gas, stars, and dark matter can be written as

$$
  P_{\rm tot} = \rho_{\rm gas} \sigma_{\rm eff}^2=\frac{\Sigma_{\rm gas}}{2 H_{\rm gas}} \sigma_{\rm eff}^2
  = {\cal W}_\mathrm{tot}={\cal W}_\mathrm{g} + {\cal W}_* + {\cal W}_\mathrm{d}.
$$

where

$$
{\cal W}_\mathrm{g}= \frac{\pi G \Sigma_{\rm gas}^2}{2},
$$

$$
{\cal W}_* \approx \pi G \Sigma_{\rm gas} \Sigma_* \frac{H_{\rm gas}}{H_{\rm gas} + H_*},
$$

and

$$
{\cal W}_\mathrm{d}  = \zeta_\mathrm{d} \Sigma_{\rm gas} \Omega_d^2 H_{\rm gas}
$$

are the weights (per unit area) of the gas in the disk within the gravitational potential of gas, stars, and dark matter, respectively.

Given $\Sigma_{\rm gas}$, $\Sigma_*$, $H_*$, and $\Omega_d$, together with an assumption for $\sigma_{\rm eff}$, we solve for $H_{\rm gas}$.
We use a constant $\sigma_{\rm eff}$ or a model for $\sigma_{\rm eff}$ calibrated from numerical simulations.

## Constant velocity dispersion
For given $\sigma_{\rm eff}$, the equation for $H_{\rm gas}$ becomes a cubic. Using

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
$h = H_{\rm gas}/H_{\rm gas-only}$,
$s_* = \Sigma_*/\Sigma_{\rm gas}$, $\eta_* = H_*/H_{\rm gas-only}$, and $\eta_d = H_{\rm dm-only}/H_{\rm gas-only}$.
The only positive real solution of the cubic equation is

$$
H_{\rm gas} = 2\sqrt{-Q}\cos\left(\frac{\theta}{3}\right), \theta = \cos^{-1}\left(\frac{R}{\sqrt{-Q}}\right)
$$
where
$Q = (3a_1-a_2^2)/9$ and $R=(9a_2a_1-27a_0-2a_2^3)/54$ with $a_i$ for the coefficients of $h^i$ for $i=0$, 1, 2.

### Analytic solutions for star-gas disks

By ignoring the dark mater contribution, one can obtain

$$ h = \frac{2\eta_*}{\eta_*-1+[(\eta_*+1)^2+8s_*\eta_*]^{1/2}} $$

or

$$ H_{\rm gas} = \frac{2\sigma_{\rm eff}^2}{\pi G \Sigma_{\rm gas} - \frac{\sigma_{\rm eff}^2}{H_*} +\left[\left(\pi G \Sigma_{\rm gas} + \frac{\sigma_{\rm eff}^2}{H_*}\right)^2+8\pi G \Sigma_*\frac{\sigma_{\rm eff}^2}{H_*}\right]^{1/2}} $$

The weight by stars is then

$$ {\cal W}_*\approx\pi G \Sigma_{\rm gas} \Sigma_* \frac{H_{\rm gas}}{H_{\rm gas} + H_*} =
\frac{2\pi G\Sigma_{\rm gas} \Sigma_* }{1+\pi G\Sigma_{\rm gas}\frac{H_*}{\sigma_{\rm eff}^2}+\left[\left(1+\pi G\Sigma_{\rm gas}\frac{H_*}{\sigma_{\rm eff}^2}\right)^2+8\pi G \Sigma_* \frac{H_*}{\sigma_{\rm eff}^2}\right]^{1/2}}. $$

Or, the total weight is

$$ {\cal W}_{\rm tot} = {\cal W}_{\rm g}+{\cal W}_* = {\cal W}_{\rm gas}\left(1+\frac{4s_*}{1+\eta_* + [(1+\eta_*)^2+8s_*\eta_*]^{1/2}}\right). $$

This is slightly different from Equation (12) in {cite}`Hassan_et_al` as the first two terms in the denominator of ${\cal W}_*$ has $1+\eta_*$ rather than 1. Equation (12) in {cite}`Hassan_et_al` is obtained by using $H_{\rm gas}$ in the star only limit (i.e., $\pi G \Sigma_{\rm gas} \rightarrow0$ or $\eta_*\rightarrow0$), which yields

$$ h\rightarrow \frac{2\eta_*}{-1+(1+8s_*\eta_*)^{1/2}} =
 \frac{1}{4s_*}(1+(1+8s_*\eta_*)^{1/2}) $$

or

$$ H_{\rm gas} = \frac{\sigma_{\rm eff}^2}{4\pi G \Sigma_*}\left(1+\left[1+8\pi G\Sigma_* \frac{H_*}{\sigma_{\rm eff}^2}\right]^{1/2}\right) $$

and

$${\cal W}_*\approx \frac{2\pi G\Sigma_{\rm gas} \Sigma_* }{1+\left[1+8\pi G \Sigma_* \frac{H_*}{\sigma_{\rm eff}^2}\right]^{1/2}}. $$

In the limit of thick stellar disk ($H_*\gg H_{\rm gas},H_{\rm gas-only}$ or $\eta_*\gg1$), a more formal approximation for the weight by stars is

$$ {\cal W}_* \rightarrow \pi G \Sigma_{\rm gas} \Sigma_*\frac{h}{\eta_*} = \frac{4\rho_*\sigma_{\rm eff}^2}{1+\left(1+\frac{16\rho_*\sigma_{\rm eff}^2}{\pi G \Sigma_{\rm gas}^2}\right)^{1/2}}\rightarrow\Sigma_{\rm gas}\sigma_{\rm eff}\sqrt{\pi G\rho_*}$$

The last limit assumes an general dominance of stellar component ($\Sigma_*>\Sigma_{\rm gas}$).

Here, the conversion between $\Sigma_*$, $H_*$, and $\rho_*$ is from the definition of $H_*$:

$$ H_* \equiv \frac{\Sigma_*}{2\rho_*}. $$

```{important}
In observations, a simplified weight formula has been widely used (e.g., {cite}`2020ApJ...892..148S`, {cite}`2021MNRAS.503.3643B`):

$$ {\cal W}_{\rm tot}\approx \frac{\pi G \Sigma_{\rm gas}^2}{2} + \Sigma_{\rm gas}\sigma_{\rm eff}\sqrt{2 G\rho_{\rm sd}}$$

where $\rho_{\rm sd} = \rho_* + \rho_{\rm dm}$.
This is equivalent to the limit of $H_*\gg H_{\rm gas}$ with $\Sigma_*>\Sigma_{\rm gas}$ discussed above.

In practical applications, often the relation between stellar disk scale length $l_*$ and scale height $z_*$; $z_*=l_*/7.3$ or $\rho_* = \Sigma_*/(4 z_*) = \Sigma_*/(0.54 l_*)$, which is based on the assumption of the isothermal (or ${\rm sech}^2$) stellar density profile. It is important to note that our definition of the stellar scale height yields $H_*=2z_*$. Therefore, when one wants to switch to the more general formula using $H_*$, it is important to use the correct empirical relation for

$$H_* = (2/7.3) l_*$$
```

## Pressure-dependent velocity dispersion
Using the TIGRESS-classic suite {cite}`2022ApJ...936..137O`, we  obtained a calibration for $\sigma_{\rm eff}$ as a function of the total pressure.

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

In either case, we can first solve the cubic equation for a constant $\sigma_{\rm eff}$ to get $H_{\rm gas}$ and $P_{\rm tot}=W_{\rm tot}$. For the non-constant case, we then update $\sigma_{\rm eff}$ based on the new $P_{\rm tot}$ and then solve for $H_{\rm gas}$ and $P_{\rm tot}$ again. We iteratively solve until $H_{\rm gas}$ is converged.

## References

The most general form of the vertical dynamical equilibrium used here is derived in {cite}`Hassan_et_al`.
See also {cite}`2022ApJ...936..137O`, {cite}`2015ApJ...815...67K`, {cite}`2013ApJ...776....1K`,
{cite}`2011ApJ...743...25K`, {cite}`2011ApJ...731...41O`, {cite}`2010ApJ...721..975O`.
