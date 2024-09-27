# Welcome to the PRFM cookbook

[![Jupyter Book Badge](https://jupyterbook.org/badge.svg)](https://changgoo.github.io/PRFM/)

This is a cookbook to guide you in applying the PRFM theory for star formation.

## Motivation

Recent advances in theoretical and computational modeling of the star-forming ISM provide an opportunity to develop more sophisticated subgrid model treatments in cosmological simulations.  A key goal of the [Learning the Universe Simons Collaboration](https://www.learning-the-universe.org/) is to replace the current, _empirically-calibrated_ prescriptions for star formation with new subgrid models that are instead calibrated from _full-physics ISM simulations_. In this regard, the **pressure-regulated feedback-modulated (PRFM)** star formation theory has been developed ({cite:t}`2022ApJ...936..137O`; see also {cite:t}`2010ApJ...721..975O`, {cite:t}`2011ApJ...731...41O`, {cite:t}`2011ApJ...743...25K`).

The theoretical concepts of the PRFM theory assume two equilibria -- [vertical dynamical equilibrium](vertical-de) that connects galactic conditions (e.g., gas and stellar surface density, gas and stellar scale height, dark matter density) and the ISM pressure/stress and [thermal/dynamical equilibrium](feedback-yield) that connects the ISM pressure/stress with star formation rate surface density. The assumption of vertical dynamical equilibrium thus allows an access to the ISM energetical parameters (pressure/stress) from galactic properties that are observable and robustly measurable in low-resolution simulations. The ratio of the ISM pressure and SFR surface density -- also called **feedback yield** -- is the key quantity in the PRFM theory and encapsulates complicated interaction between ISM physics and stellar feedback. The feedback yield is thus best calibrated by the high-resolution, full physics ISM simulations, which we use the TIGRESS-classic and TIGRESS-NCR simulation suites here.
<!-- A vertical dynamical equilibrium state of the ISM in disk galaxies demands a certain level of total pressure to support disk weight. Each component of ISM pressure/stress is maintained by massive star feedback by heating the gas, driving turbulence, and activating galactic dynamo to offset losses through cooling and turbulence dissipation. From the former, one can connect galactic condtions with the required total pressure. From the latter, one can connect each pressure component (and hence total pressure) and star formation rate, which is characterized by [feedback yields](feedback-yield). -->



```{tableofcontents}
```
