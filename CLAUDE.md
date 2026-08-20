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

**Run tests (PHANGS data-dependent tests skip automatically if data is absent):**
```bash
pytest tests/ -v
```

**Run a single test file:**
```bash
pytest tests/test_prfm.py -v
```

**Run all tests including integration (requires downloaded PHANGS data):**
```bash
pytest tests/ -v
```

Integration tests are gated by `HAS_DATA` in `tests/test_phangs.py` (a
`skipif` on the presence of `data/phangs_megatable/`), not by a registered
pytest marker. Therefore, `-m "not integration"` does not deselect them when
the PHANGS data are present.

**Download PHANGS megatable data (needed for integration tests + PHANGS notebooks):**
```bash
python scripts/download_phangs.py
```

**Run the TIGRESS-PHANGS simulation-suite workflow (this branch):**
```bash
project/scripts/run_workflow.sh -s 10 -n 32     # sample → plot → YAML → SLURM
project/scripts/run_workflow.sh -h              # full option list
```

## Theory and Physical Background

PRFM theory consists of vertical dynamical equilibrium and feedback yield,
with an effective-equation-of-state calibration. The primer is under `book/`:

- `book/vertical-de.md` -- vertical dynamical equilibrium
- `book/feedback-yield.md` -- feedback yields
- `book/equation_of_state.md` -- effective equation of state

The vertical-equilibrium treatment follows three papers stored under
`references/`: Ostriker & Kim (2022) for the thick-stellar-disk formulation,
Hassan et al. (2024) for finite stellar-disk thickness, and Jeffreson et al.
(2026) for the spherical-component generalization. Spherical density is
converted to vertical harmonic frequency with
`Omega_d**2 = 2*pi*G*a_d*rho_dm`: use `a_d=2` (default) for the flat-rotation
convention, `a_d=1` for an NFW-like halo, and `a_d=2/3` for a Hernquist bulge.

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
- Flexible stellar inputs (surface or volume density) and spherical-component
  inputs (vertical harmonic frequency or density, with configurable `a_d`)
- Built-in model management

### Data handling (`prfm/simulations.py`)

`PRFM_data` container for simulation outputs and observational data, with log/linear conversion and uncertainty propagation.

### Bundled data

- `prfm/prfm_ring.csv` — ring model parameters
- `prfm/tigress_ncr_K24.nc` — TIGRESS-NCR simulation data (Kim et al. 2024)
- `data/all_data_z0_selected.pkl` — processed galaxy/simulation data used in notebooks

## Development Rules

- **TDD**: Write tests before implementation; unit tests must not require network or disk data.
- **Plan first**: For any non-trivial feature, enter plan mode and get alignment before coding.
- **Ask clarifying questions** before starting work when requirements are ambiguous.
- **Type hints**: All Python functions must have type annotations.
- **Branch workflow**: New development goes on a feature branch; open a WIP PR with a task checklist immediately after planning.
- **Frequent commits**: Commit after each self-contained chunk of work — do not bundle unrelated changes.
- **Modular Design**: When writing python scripts, any helper functions should be in a python file such that reusable.

## matplotlib plotting tips

- **fontsize**: never use the absolute value. use medium, small, x-small, etc.
