# Effective Equation of State

Using the TIGRESS-classic suite ({cite:t}`2022ApJ...936..137O`), we obtained a calibration for an effective equation of state -- a relation between the effective pressure and density $P_{\rm eff}=\rho\sigma_{\rm eff}^2$. We provided simple power-law fitting results using two characterizations. (1) The relation between the volume-weighted averages of the midplane total pressure and density of the warm and cold gas ($P_{\rm tot,2p}$ and $n_{\rm H, 2p}=\rho_{\rm 2p}/1.4m_H$ in Table 2 of OK22). (2) The relation between the mass-weighted vertical velocity dispersion of the warm and cold phase and the dynamical equilibrium pressure ($\overline{\sigma}_{\rm eff,2p}$ and $P_{\rm DE}$ in Table 2 of OK22). From (1), we calculated the midplane effective velocity dispersion $\sigma_{\rm eff,mid}\equiv(P_{\rm tot,2p}/\rho_{\rm 2p})^{1/2}$. From (2), we calculated the mass-weighted average density $\rho_{\rm avg}=P_{\rm DE}/\overline{\sigma}_{\rm eff,2p}^2$ and the press-density relation. Therefore, both provide same information; they just characterized differently.

In what follows, we will use $P_{\rm eff}$ and $n_H$ for effective pressure and number density. We do not distinguish $P_{\rm tot,2p}$, $P_{\rm DE}$ and $\mathcal{W}_{\rm tot}$ as they are proven to be consistent with each other. We will use the nomenclature $\overline{\sigma}_{\rm eff,2p}\to\sigma_{\rm eff,avg}$ to be consistent with the TIGRESS-NCR results.

```{important}
The power-law fit of the effective equation of state will not valid to the arbitrary low pressure and density. We show the fitting results and actual simulation data points. We caution that especially at low pressure and density regimes, the effective velocity dispersions cannot be much lower than the warm medium sound speed ($\sim 10{\rm km/s}$). It would be safe (and more physically meaningful) to apply a velocity floor of this level to the effective velocity dispersion and corresponding effective pressure when adopting this result. See [this example](eos_sims.ipynb).
```
## TIGRESS-classic calibration

See Table 2, Section 4.6, and Figure 14 in {cite:t}`2022ApJ...936..137O`.

### Midplane Measures

(1) Pressure-Density Relation

$$ P_{\rm eff} = 2\times10^4 k_B {\rm cm^{-3}\,K}\; n_H^{1.43} $$

or

$$ \log P_4 = 1.43 \log n_0 +0.3 $$

where $P_4\equiv(P_{\rm eff}/k_B\,10^4{\rm\,cm^{-3}\,K})$ and $n_0\equiv n_{\rm H}/cm^{-3}$.

(2) Effective Velocity Dispersion

$$
\sigma_{\rm eff,mid} = 9.8{\rm km/s}\; P_4^{0.15} = 10.9{\rm km/s}\; n_0^{0.2}
$$

### Mass-weighted Averages

(1) Pressure-Density Relation

$$ P_{\rm eff} = 5\times10^4 k_B {\rm cm^{-3}\,K}\; n_H^{1.8} $$

or

$$ \log P_4 = 1.8 \log n_0 + 0.7 $$

(2) Effective Velocity Dispersion

$$
\sigma_{\rm eff} = 12{\rm km/s}\; P_4^{0.22} = 17{\rm km/s}\; n_0^{0.4}
$$

## TIGRESS-NCR calibration

Using the TIGRESS-NCR suite ({cite:t}`2024ApJ...972...67K`), we  obtained a calibration for $\sigma_{\rm eff}$ as a function of the total pressure and gas metallicity, but the latter dependence is found to be very weak. So, here we omit the metallicity dependence.

See Tables 2 and 3, Figures 11 and 12, and Section 4.4.

### Midplane Measures

(1) Pressure-Density Relation

$$ P_{\rm eff} = 1.4\times10^4 k_B {\rm cm^{-3}\,K}\; n_H^{1.2} $$

or

$$ \log P_4 = 1.2 \log n_0 + 0.15 $$

(2) Effective Velocity Dispersion

$$
\sigma_{\rm eff,mid} = 8.9{\rm km/s}\; P_4^{0.08} = 9.2{\rm km/s}\; n_0^{0.1}
$$


### Mass-weighted Averages

(1) Pressure-Density Relation

$$ P_{\rm eff} = 3.1\times10^4 k_B {\rm cm^{-3}\,K}\; n_H^{1.3} $$

or

$$ \log P_4 = 1.3 \log n_0 + 0.5 $$

(2) Effective Velocity Dispersion

$$
\sigma_{\rm eff,avg} = 11.7{\rm km/s}\; P_4^{0.12} = 13.4{\rm km/s}\;n_0^{0.16}
$$

<!-- In either case, we can first solve the cubic equation for a constant $\sigma_{\rm eff}$ to get $H_{\rm gas}$ and $P_{\rm tot}=W_{\rm tot}$. For the non-constant case, we then update $\sigma_{\rm eff}$ based on the new $P_{\rm tot}$ and then solve for $H_{\rm gas}$ and $P_{\rm tot}$ again. We iteratively solve until $H_{\rm gas}$ is converged. -->
