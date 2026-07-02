#!/usr/bin/env bash
#
# End-to-end pipeline: PHANGS-informed KDE-Sobol sample -> design CSV ->
# diagnostic figures -> suite YAML -> TIGRESS-NCR SLURM scripts.
#
# Usage:
#   project/scripts/run_workflow.sh [options]
#
# Options (all have defaults):
#   -s SIGMA_GAS   Target Sigma_gas (M_sun/pc^2)              [default: 10]
#   -d DELTA        Band half-width (dex)                     [default: 0.3]
#   -n N            Sample size (rows in the design)          [default: 32]
#   -b BASE         TIGRESS-NCR base model                    [default: R8_8pc]
#   -m MACHINE_DIR  Machine YAML directory                    [default: project/suites/machines/stellar]
#   -q QUEUE        Queue name within machine dir             [default: standard]
#   -r RUN_BASE     Scratch root for RUNDIR (overrides YAML)  [default: unchanged]
#   -h              Show this help and exit
#
# Environment:
#   ATHENA_TIGRESS_DIR  Path to Athena-TIGRESS repo
#                       (default: $HOME/Sources/Athena-TIGRESS)
#
# Outputs:
#   project/output/design_Sgas<S>_n<NNNN>.csv                design CSV
#   project/output/design_Sgas<S>_n<NNNN>_{selection,corner,marginals}.png
#                                                            diagnostic figures
#   project/suites/design_Sgas<S>_n<NNNN>.yml                suite YAML
#   project/suites/slurms/design_Sgas<S>_n<NNNN>/*.slurm     SLURM scripts

set -euo pipefail

# ── defaults ────────────────────────────────────────────────────────────────
SIGMA_GAS=10
DELTA=0.3
N_SAMPLES=32
BASE=R8_8pc
MACHINE_DIR=project/suites/machines/stellar
QUEUE=standard
RUN_BASE=""

usage() { sed -n '2,25p' "$0"; exit 0; }

while getopts ":s:d:n:b:m:q:r:h" opt; do
    case "$opt" in
        s) SIGMA_GAS=$OPTARG ;;
        d) DELTA=$OPTARG ;;
        n) N_SAMPLES=$OPTARG ;;
        b) BASE=$OPTARG ;;
        m) MACHINE_DIR=$OPTARG ;;
        q) QUEUE=$OPTARG ;;
        r) RUN_BASE=$OPTARG ;;
        h) usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 2 ;;
    esac
done

# ── locate project root ─────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

ATHENA_DIR="${ATHENA_TIGRESS_DIR:-$HOME/Sources/Athena-TIGRESS}"
GENERATE_SLURM="$ATHENA_DIR/scripts/generate_slurm.py"

if [[ ! -f "$GENERATE_SLURM" ]]; then
    echo "ERROR: generate_slurm.py not found at $GENERATE_SLURM" >&2
    echo "       Set ATHENA_TIGRESS_DIR to override the location." >&2
    exit 1
fi

# ── derived paths ───────────────────────────────────────────────────────────
SG_TAG=$(printf "%.1f" "$SIGMA_GAS")
N_TAG=$(printf "%04d" "$N_SAMPLES")
STEM="design_Sgas${SG_TAG}_n${N_TAG}"
CSV_PATH="project/output/${STEM}.csv"
YAML_PATH="project/suites/${STEM}.yml"
SLURM_DIR="project/suites/slurms/${STEM}"

echo "──────────────────────────────────────────────────────────────"
echo "TIGRESS-PHANGS pilot workflow"
echo "  Sigma_gas target : ${SIGMA_GAS} M_sun/pc^2"
echo "  band half-width  : ${DELTA} dex"
echo "  sample size      : ${N_SAMPLES}"
echo "  base model       : ${BASE}"
echo "  machine          : ${MACHINE_DIR} (${QUEUE})"
echo "──────────────────────────────────────────────────────────────"

# ── step 1: sample ──────────────────────────────────────────────────────────
echo
echo "[1/4] KDE-Sobol sampling ..."
python project/scripts/run_sampling.py \
    --sigma-gas "$SIGMA_GAS" \
    --delta     "$DELTA" \
    --n-samples "$N_SAMPLES"

# ── step 2: plot sample against PHANGS distribution ─────────────────────────
echo
echo "[2/4] Plotting design against PHANGS distribution ..."
python project/scripts/plot_design_from_csv.py \
    "$CSV_PATH" \
    --sigma-gas "$SIGMA_GAS" \
    --delta     "$DELTA"

# ── step 3: convert CSV -> suite YAML ───────────────────────────────────────
echo
echo "[3/4] Converting CSV -> suite YAML ..."
python project/scripts/csv_to_slurm_yaml.py \
    "$CSV_PATH" \
    --base   "$BASE" \
    --output "$YAML_PATH"

if [[ -n "$RUN_BASE" ]]; then
    # Overwrite the run_base in the generated YAML in-place.
    python - <<PY
import yaml, sys
path = "$YAML_PATH"
with open(path) as f:
    cfg = yaml.safe_load(f)
cfg.setdefault("suite", {})["run_base"] = "$RUN_BASE"
with open(path, "w") as f:
    yaml.safe_dump(cfg, f, sort_keys=False, default_flow_style=False,
                   width=100, indent=2)
print(f"  overrode suite.run_base -> $RUN_BASE")
PY
fi

# ── step 4: generate SLURM scripts ──────────────────────────────────────────
echo
echo "[4/4] Generating SLURM scripts ..."
mkdir -p "$SLURM_DIR"
python "$GENERATE_SLURM" \
    "$YAML_PATH" \
    --machine "$MACHINE_DIR" \
    --queue   "$QUEUE" \
    --output-dir "$SLURM_DIR"

N_SLURM=$(ls -1 "$SLURM_DIR"/*.slurm 2>/dev/null | wc -l | tr -d ' ')
echo
echo "──────────────────────────────────────────────────────────────"
echo "Done. Summary:"
echo "  CSV       : $CSV_PATH"
echo "  Figures   : project/output/${STEM}_{selection,corner,marginals}.png"
echo "  Suite YAML: $YAML_PATH"
echo "  SLURM dir : $SLURM_DIR   ($N_SLURM scripts)"
echo "──────────────────────────────────────────────────────────────"
