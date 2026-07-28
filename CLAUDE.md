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

**Run tests (unit only, no data required):**
```bash
pytest tests/ -v -m "not integration"
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
pytest marker — `-m "not integration"` still selects them because the marker
is applied to each test.

**Download PHANGS megatable data (needed for integration tests + PHANGS notebooks):**
```bash
python scripts/download_phangs.py
```

**Run the TIGRESS-PHANGS simulation-suite workflow (this branch):**
```bash
project/scripts/run_workflow.sh -s 10 -n 32     # sample → plot → YAML → SLURM
project/scripts/run_workflow.sh -h              # full option list
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

### PHANGS observational layer (`prfm/phangs*.py`)

Added to apply PRFM to PHANGS megatable data (Sun et al. 2022/2023) and to
design TIGRESS-NCR simulation suites from it:

- `prfm/phangs.py` — load/download the PHANGS megatable (per-galaxy ECSV files
  under `data/phangs_megatable/`, fetched from CANFAR), join apertures
  (`annulus`/`gauss`/`hexagon`), and compute PRFM input columns.
- `prfm/phangs_sampling.py` — `SamplingConfig` + KDE-Sobol (and LHS/expanded)
  conditional sampling of the 5 design fields (`Sigma_gas`, `Sigma_star`,
  `H_star`, `Omega`, `qshear`) within a Σ_gas band. Fixed Sobol seed for
  reproducible, nested designs.
- `prfm/phangs_plot.py` — shared plotting/data-prep helpers for exploratory
  figures; not part of the public computation API.

`config/phangs_prfm.yml` drives megatable loading (data dir, aperture
join keys, canonical column choices) for `book/prfm_phangs.ipynb`.

## Simulation-suite workflow (`project/`, this branch only)

`project/scripts/run_workflow.sh` is an end-to-end pipeline: PHANGS-informed
KDE-Sobol sample → diagnostic plots → suite YAML → SLURM scripts for a
TIGRESS-NCR run. Stages are individually runnable via `-t sample,plot,yaml,slurm`.
See `project/WORKFLOW.md` for the full data-flow diagram, CSV→athinput key
mapping, machine-YAML layout, and reproducibility notes. The final SLURM step
shells out to `$ATHENA_TIGRESS_DIR/scripts/generate_slurm.py` in the separate
`Athena-TIGRESS` repo (set `ATHENA_TIGRESS_DIR` if not at `$HOME/Sources/Athena-TIGRESS`).

Ad-hoc PHANGS data-exploration scripts (download, column inspection, field and
correlation exploration) live in top-level `scripts/`.

### Bundled data

- `prfm/prfm_ring.csv` — ring model parameters
- `prfm/tigress_ncr_K24.nc` — TIGRESS-NCR simulation data (Kim et al. 2024)
- `prfm/ncr_eos.nc`, `prfm/classic_eos.nc` — effective EOS tables (σ_eff models)
- `data/all_data_z0_selected.pkl` — processed galaxy/simulation data used in notebooks
- `data/phangs_megatable/*.ecsv` — per-galaxy PHANGS megatables (downloaded, gitignored)

Note: `setup.py` `package_data` only ships `prfm_ring.csv` and
`tigress_ncr_K24.nc`; the `.nc` EOS tables are used from the source tree.

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