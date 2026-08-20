# PRFM Theory Primer

The pressure-regulated, feedback-modulated (PRFM) framework connects galactic
environmental conditions to the pressure and star formation rate of the
multiphase interstellar medium. Its two central ingredients are vertical
dynamical equilibrium, which sets the required midplane pressure, and feedback
yield, which relates that pressure to the star formation rate. An effective
equation of state closes the model by specifying the turbulent, thermal, and
magnetic support calibrated from resolved ISM simulations.

## Generalized spherical gravity

The vertical-equilibrium model now treats dark matter halos and bulges as
general spherical components. Their density is connected to the vertical
harmonic frequency by

$$
\Omega_{\rm sph}^2 = 2\pi G a_d\rho_{\rm sph}.
$$

The dimensionless coefficient $a_d$ makes the assumed mass profile explicit:

| Profile convention | $a_d$ |
| --- | ---: |
| Flat rotation curve (default) | $2$ |
| NFW-like halo | $1$ |
| Hernquist bulge | $2/3$ |

The default preserves the package's previous conversion
$\Omega_{\rm sph}^2=4\pi G\rho_{\rm sph}$. When `PRFM` is initialized with
`rho_dm`, the `a_d` argument controls this conversion. When `Omega_d` is given
directly, it is already the vertical harmonic frequency and `a_d` does not
change the equilibrium solution. For separately modeled halo and bulge
components, their contributions add as $\Omega_{\rm sph}^2=\sum_i\Omega_i^2$.

The public conversion helpers are
`get_omega_spherical_from_density(rho, a_d)` and
`get_density_from_omega_spherical(Omega_sph, a_d)`. The existing `rho_dm` and
`Omega_d` argument names are retained for API compatibility.

## Contents

- [Vertical dynamical equilibrium](vertical-de.md) develops the gas, stellar,
  and spherical-component weight terms, including finite stellar-disk
  thickness.
- [Feedback yield](feedback-yield.md) relates equilibrium pressure to the star
  formation rate.
- [Effective equation of state](equation_of_state.md) summarizes the calibrated
  pressure-support model.
- [API reference](api.md) lists the corresponding package interfaces.

## Theoretical development

The formulation follows three stages represented by the papers under
`../references/`:

- Ostriker & Kim (2022) establishes the PRFM framework and the commonly used
  thick-stellar-disk vertical-equilibrium approximation, calibrated with
  TIGRESS-classic simulations.
- Hassan et al. (2024) generalizes the vertical-equilibrium solution to finite
  stellar scale height.
- Jeffreson et al. (2026) extends the external-gravity treatment to include
  bulges and other spherical components.
