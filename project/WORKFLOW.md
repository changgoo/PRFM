# TIGRESS-PHANGS Simulation Suite Workflow

End-to-end pipeline from PHANGS-informed KDE-Sobol sampling to SLURM scripts
ready to submit for a TIGRESS-NCR simulation suite.

## One-command pipeline

```bash
project/scripts/run_workflow.sh           # defaults: Sigma_gas=10, n=32, R8_8pc
project/scripts/run_workflow.sh -s 20 -n 64
project/scripts/run_workflow.sh -h        # full option list
```

Options (all optional, sensible defaults):

| Flag | Meaning | Default |
|------|---------|---------|
| `-s SIGMA_GAS`  | Target Σ_gas [M_sun/pc²]            | `10` |
| `-d DELTA`      | Band half-width in dex               | `0.3` |
| `-n N`          | Sample size (powers of 2 preferred)  | `32` |
| `-b BASE`       | TIGRESS-NCR base model               | `R8_8pc` |
| `-m MACHINE_DIR`| Machine YAML directory               | `project/suites/machines/stellar` |
| `-q QUEUE`      | Queue name within machine dir        | `standard` |
| `-r RUN_BASE`   | Scratch root for RUNDIR              | value from machine YAML |

Environment: set `ATHENA_TIGRESS_DIR` if `Athena-TIGRESS` is not at
`$HOME/Sources/Athena-TIGRESS`.

## Pipeline steps

```
                    +--------------------+
    PHANGS data --> | 1. run_sampling.py |--+
    (config)        +--------------------+  |
                                            v
                            +--------------------------+
                            | design_Sgas<S>_n<NNNN>.csv|
                            +--------------------------+
                             |                        |
                             |                        v
                             |  +----------------------------+
                             |  | 2. plot_design_from_csv.py  |
                             |  +----------------------------+
                             |               |
                             |               v
                             |    diagnostic PNGs
                             v
                    +--------------------------+
                    | 3. csv_to_slurm_yaml.py  |
                    +--------------------------+
                                 |
                                 v
                    +--------------------------+
                    | design_Sgas<S>_n<NNNN>.yml|
                    +--------------------------+
                                 |
                                 v
                    +--------------------------------+
                    | 4. generate_slurm.py            |
                    |   (Athena-TIGRESS/scripts/)     |
                    +--------------------------------+
                                 |
                                 v
                    slurms/<STEM>/*.slurm  (one per row)
```

### Step 1 — sample

`project/scripts/run_sampling.py` fits a KDE to the PHANGS reference pixels in
the Σ_gas band and draws a scrambled Sobol sequence in the 5 design fields:
`Sigma_gas`, `Sigma_star`, `H_star`, `Omega`, `qshear`. Physical bounds
(`0 < q ≤ 1.5`) are enforced by rejection with sequential replacement.

Output: `project/output/design_Sgas<S>_n<NNNN>.csv`.

Also writes `project/output/fairness_summary.csv` with the per-field quantile
mismatch ε_a in dex.

### Step 2 — plot design against PHANGS

`project/scripts/plot_design_from_csv.py` visualises the exact CSV that will
be handed to TIGRESS-NCR. Three figures:

- `<stem>_selection.png` — full PHANGS distribution (gray) + selected Σ_gas
  band (blue) + design points (red). Answers "is my design where I think it is
  in the full observational parameter space?"
- `<stem>_corner.png` — reference vs design zoomed to the band. Answers "does
  my design cover the reference set well?"
- `<stem>_marginals.png` — 1-D marginals.

### Step 3 — CSV → suite YAML

`project/scripts/csv_to_slurm_yaml.py` maps CSV columns to TIGRESS-NCR
athinput keys:

| CSV column   | athinput key       |
|--------------|--------------------|
| `Sigma_gas`  | `problem/surf`     |
| `Sigma_star` | `problem/SurfS`    |
| `H_star`     | `problem/zstar`    |
| `Omega`      | `problem/Omega`    |
| `qshear`     | `problem/qshear`   |

Also emits fixed parameters that don't vary across the suite:

| Key                | Value                              |
|--------------------|------------------------------------|
| `problem/rhodm`    | `0.0064` M_sun/pc³ (R8 fiducial)   |
| `problem/Zgas`     | `1.0` (solar metallicity)          |
| `problem/Zdust`    | `1.0` (solar dust-to-gas)          |

Row `i` becomes one model with `suffix=rowNNNN`, so the resulting job IDs are
`<base>_NCR_row0000` … `<base>_NCR_rowNNNN` (short, sortable). Varying params
are written to `extra_overrides` (so they do NOT clutter the job ID); fixed
params go in the top-level `fixed_overrides` block.

Output: `project/suites/design_Sgas<S>_n<NNNN>.yml`.

### Step 4 — generate SLURM scripts

The driver calls
`$ATHENA_TIGRESS_DIR/scripts/generate_slurm.py` with:

- the suite YAML from step 3
- a machine directory (`--machine`) containing per-queue YAML files
- a queue name (`--queue`), default `standard`

Output: `project/suites/slurms/<STEM>/<jobid>.slurm`, one script per design
point.

## Machine YAML layout

`generate_slurm.py` expects `<machine>/<queue>.yml`:

```
project/suites/machines/
  stellar/
    standard.yml
    # gpu.yml         (add as needed)
```

Each queue file supplies scheduler/environment defaults (nodes, tasks_per_node,
walltime, modules, run_base, mail_user). Suite YAML values from step 3 take
precedence for anything both files set.

Add a new machine by creating `project/suites/machines/<name>/<queue>.yml`
following the Stellar template.

## Submitting jobs

Each generated SLURM script is a two-mode wrapper:

```bash
# On the cluster, from the slurm output directory:
sbatch R8_8pc_NCR_row0000.slurm -i     # fresh start
sbatch R8_8pc_NCR_row0000.slurm -r     # restart from latest .rst
```

If `slurm.auto_resubmit: true` is set (default in the Stellar template), each
job self-resubmits on wall-time exit until `tlim` is reached.

## Reproducibility

- **Sobol seed** is fixed to `42` by `run_sampling.py`; the same call always
  produces the same CSV.
- **Nesting**: the first 64 rows of the n=128 design are byte-identical to the
  n=64 design (verified by end-to-end test), so extending a suite from n=64 to
  n=128 requires running only the 64 new simulations.
- Row indices in the CSV correspond one-to-one with `rowNNNN` job suffixes,
  so provenance from a completed simulation back to the PHANGS-informed
  parameters is direct.

## Typical usage patterns

**First pilot at a new Σ_gas:**
```bash
project/scripts/run_workflow.sh -s 5 -n 32
```

**Scale-up (nested — old sims still valid):**
```bash
project/scripts/run_workflow.sh -s 10 -n 128
# rows 0-63 are identical to a prior n=64 design at Sigma_gas=10.
```

**Different resolution:**
```bash
project/scripts/run_workflow.sh -s 10 -n 32 -b R8s_4pc
```

**Different cluster:**
```bash
project/scripts/run_workflow.sh -m project/suites/machines/anvil -q standard
```

## Related documents

- `project/TIGRESS-PHANGS-ILI-project.md` — overall project description
- `project/TIGRESS-PHANGS-sampling-design.md` — sampling design rationale
- `project/TIGRESS-PHANGS-sobol-sampling.md` — Sobol formulation details
- `project/doc/paper1-suite/` — Paper I draft (simulation suite methodology)
- `project/TODO.md` — outstanding items
