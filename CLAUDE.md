# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**PRFM** is a Python toolkit implementing the Pressure-Regulated Feedback Model for calculating vertical equilibrium and star formation rates in galactic gas disks.

The main package is `prfm/` (NumPy-based). A JAX port lives on the `jax-experimental` branch (`prfm_jax/`) — see `prfm_jax/NOTES.md` there for status and next steps.

## Commands

**Install:**
```bash
pip install -r requirements.txt
pip install .
```

**Lint/Format:**
```bash
ruff check prfm/
ruff format prfm/
```

**Pre-commit (lint + format):**
```bash
pre-commit run --all-files
```

**Build documentation (Jupyter Book):**
```bash
jupyter-book build book/
```

**Run tests:**
```bash
pytest tests/ -v
```

## Architecture

### Core computation layer (`prfm/prfm.py`)

Pure functions organized around disk physics:

- **Scale height solvers**: `get_scale_height_gas_only`, `get_scale_height_star_only`, `get_scale_height_dm_only`, `get_scale_height_star_gas`, `get_scale_height_thin`, `get_scale_height_thick`, `get_scale_height_numerical`
- **Pressure/weight**: `get_weight_gas/star/dm`, `get_pressure`, `get_weights`
- **Feedback & SFR**: `get_sigma_eff`, `get_feedback_yield`, `get_feedback_yield_comp`, `get_sfr`, `get_self_consistent_solution`

Two embedded model dictionaries drive parameterized behavior:
- `_sigma_eff_models` — velocity dispersion models (e.g. `tigress-classic-mid`, `tigress-ncr-avg`)
- `_yield_models` — feedback yield models (e.g. `tigress-classic`, `tigress-ncr-decomp-all`)

### Object-oriented wrapper (`PRFM` class)

`PRFM` in `prfm/prfm.py` wraps the functional layer with:
- Automatic unit conversion (CGS ↔ astronomical)
- Multiple disk configurations: thin, thick, general
- Flexible input modes for stellar (surface or volume density) and dark matter (rotation speed or density)
- Built-in model management

### Data handling (`prfm/simulations.py`)

`PRFM_data` container for simulation outputs and observational data, with log/linear conversion and uncertainty propagation.

### Bundled data

- `prfm/prfm_ring.csv` — ring model parameters
- `prfm/tigress_ncr_K24.nc` — TIGRESS-NCR simulation data (Kim et al. 2024)
- `data/all_data_z0_selected.pkl` — processed galaxy/simulation data used in notebooks

