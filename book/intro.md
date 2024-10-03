# Welcome to the PRFM cookbook

[![Jupyter Book Badge](https://jupyterbook.org/badge.svg)](https://changgoo.github.io/PRFM/)

This is a cookbook to guide you in applying the PRFM theory for star formation.

## Motivation

Recent advances in theoretical and computational modeling of the star-forming ISM provide an opportunity to develop more sophisticated subgrid model treatments in cosmological simulations.  A key goal of the [Learning the Universe Simons Collaboration](https://www.learning-the-universe.org/) is to replace the current, _empirically-calibrated_ prescriptions for star formation with new subgrid models that are instead calibrated from _full-physics ISM simulations_. In this regard, the **pressure-regulated feedback-modulated (PRFM)** star formation theory has been developed ({cite:t}`2022ApJ...936..137O`; see also {cite:t}`2010ApJ...721..975O`, {cite:t}`2011ApJ...731...41O`, {cite:t}`2011ApJ...743...25K`).

The theoretical concepts of the PRFM theory assume two equilibria: [vertical dynamical equilibrium](vertical-de) that connects galactic conditions (e.g., gas and stellar surface density, gas and stellar scale height, dark matter density) with the ISM pressure/stress and [thermal/dynamical equilibrium](feedback-yield) that connects the ISM pressure/stress with star formation rate (SFR) surface density. The assumption of vertical dynamical equilibrium thus allows the calculation of the ISM energetical properties (pressure/stress) from galactic conditions that are observable or robustly measurable in low-resolution galaxy simulations. The ratio of the ISM pressure and SFR surface density -- called [**feedback yield**](feedback-yield) -- is the key quantity in the PRFM theory and encapsulates complicated interaction between ISM physics and stellar feedback. The feedback yield is thus best calibrated by the high-resolution, full physics ISM simulations, which we use the TIGRESS-classic ({cite:t}`2022ApJ...936..137O`) and TIGRESS-NCR ({cite:t}`2024ApJ...972...67K`) simulation suites here.

The PRFM theory provides a new paradigm to undersand the observational correlations between galactic conditions and SFR surface density (see [this example](prfm_phangs.ipynb)) and to construct a physically-motivated subgrid recipe for SFR (see [this example](prfm_tng.ipynb)).

## Table of Contents

```{tableofcontents}
```
