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
#   -d DELTA       Band half-width (dex)                      [default: 0.3]
#   -n N           Sample size (rows in the design)           [default: 32]
#   -b BASE        TIGRESS-NCR base model                     [default: R8_8pc]
#   -m MACHINE_DIR Machine YAML directory                     [default: $ATHENA_TIGRESS_DIR/scripts/stellar]
#   -q QUEUE       Queue name within machine dir              [default: standard]
#   -t STEPS       Steps to run (comma-separated or 'all')    [default: all]
#                  Valid steps: sample, plot, yaml, slurm
#                  Example: -t plot,yaml   (skip sampling, skip slurm)
#                  Example: -t slurm       (only regenerate SLURM scripts)
#   -h             Show this help and exit
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
MACHINE_DIR=""            # if empty, set to $ATHENA_DIR/scripts/stellar below
QUEUE=standard
STEPS="all"

usage() { sed -n '2,35p' "$0"; exit 0; }

while getopts ":s:d:n:b:m:q:t:h" opt; do
    case "$opt" in
        s) SIGMA_GAS=$OPTARG ;;
        d) DELTA=$OPTARG ;;
        n) N_SAMPLES=$OPTARG ;;
        b) BASE=$OPTARG ;;
        m) MACHINE_DIR=$OPTARG ;;
        q) QUEUE=$OPTARG ;;
        t) STEPS=$OPTARG ;;
        h) usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 2 ;;
    esac
done

# ── parse and validate steps ────────────────────────────────────────────────
ALL_STEPS=(sample plot yaml slurm)

if [[ "$STEPS" == "all" ]]; then
    SELECTED_STEPS=("${ALL_STEPS[@]}")
else
    IFS=',' read -ra RAW_STEPS <<< "$STEPS"
    SELECTED_STEPS=()
    for step in "${RAW_STEPS[@]}"; do
        step_lc=$(echo "$step" | tr '[:upper:]' '[:lower:]' | xargs)
        valid=false
        for known in "${ALL_STEPS[@]}"; do
            if [[ "$step_lc" == "$known" ]]; then
                valid=true
                break
            fi
        done
        if [[ "$valid" == "false" ]]; then
            echo "ERROR: unknown step '$step'. Valid: ${ALL_STEPS[*]}, all" >&2
            exit 2
        fi
        SELECTED_STEPS+=("$step_lc")
    done
fi

# helper: is step selected?
has_step() {
    local target="$1"
    for step in "${SELECTED_STEPS[@]}"; do
        [[ "$step" == "$target" ]] && return 0
    done
    return 1
}

# ── locate project root ─────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

ATHENA_DIR="${ATHENA_TIGRESS_DIR:-$HOME/Sources/Athena-TIGRESS}"
GENERATE_SLURM="$ATHENA_DIR/scripts/generate_slurm.py"

# Default machine dir lives alongside generate_slurm.py in Athena-TIGRESS.
if [[ -z "$MACHINE_DIR" ]]; then
    MACHINE_DIR="$ATHENA_DIR/scripts/stellar"
fi

if has_step "slurm" && [[ ! -f "$GENERATE_SLURM" ]]; then
    echo "ERROR: generate_slurm.py not found at $GENERATE_SLURM" >&2
    echo "       Set ATHENA_TIGRESS_DIR to override the location." >&2
    exit 1
fi

if has_step "slurm" && [[ ! -d "$MACHINE_DIR" ]]; then
    echo "ERROR: machine directory not found: $MACHINE_DIR" >&2
    echo "       Pass -m <dir> to override, or set ATHENA_TIGRESS_DIR." >&2
    exit 1
fi

# ── derived paths ───────────────────────────────────────────────────────────
SG_TAG=$(printf "%.1f" "$SIGMA_GAS")
N_TAG=$(printf "%04d" "$N_SAMPLES")
STEM="design_Sgas${SG_TAG}_n${N_TAG}"
CSV_PATH="project/output/${STEM}.csv"
YAML_PATH="project/suites/${STEM}.yml"
SLURM_DIR="project/suites/slurms/${STEM}"

# ── input-existence checks for downstream steps ─────────────────────────────
require_csv() {
    if [[ ! -f "$CSV_PATH" ]]; then
        echo "ERROR: $CSV_PATH not found. Run step 'sample' first," >&2
        echo "       or use -t all to run the full workflow." >&2
        exit 1
    fi
}

require_yaml() {
    if [[ ! -f "$YAML_PATH" ]]; then
        echo "ERROR: $YAML_PATH not found. Run step 'yaml' first," >&2
        echo "       or use -t all to run the full workflow." >&2
        exit 1
    fi
}

if has_step "plot" && ! has_step "sample"; then require_csv; fi
if has_step "yaml" && ! has_step "sample"; then require_csv; fi
if has_step "slurm" && ! has_step "yaml" && ! has_step "sample"; then require_yaml; fi

# ── banner ──────────────────────────────────────────────────────────────────
echo "──────────────────────────────────────────────────────────────"
echo "TIGRESS-PHANGS pilot workflow"
echo "  steps            : ${SELECTED_STEPS[*]}"
echo "  Sigma_gas target : ${SIGMA_GAS} M_sun/pc^2"
echo "  band half-width  : ${DELTA} dex"
echo "  sample size      : ${N_SAMPLES}"
echo "  base model       : ${BASE}"
echo "  machine          : ${MACHINE_DIR} (${QUEUE})"
echo "──────────────────────────────────────────────────────────────"

# ── step: sample ────────────────────────────────────────────────────────────
if has_step "sample"; then
    echo
    echo "[sample] KDE-Sobol sampling ..."
    python project/scripts/run_sampling.py \
        --sigma-gas "$SIGMA_GAS" \
        --delta     "$DELTA" \
        --n-samples "$N_SAMPLES"
fi

# ── step: plot ──────────────────────────────────────────────────────────────
if has_step "plot"; then
    echo
    echo "[plot] Plotting design against PHANGS distribution ..."
    python project/scripts/plot_design_from_csv.py \
        "$CSV_PATH" \
        --sigma-gas "$SIGMA_GAS" \
        --delta     "$DELTA"
fi

# ── step: yaml ──────────────────────────────────────────────────────────────
if has_step "yaml"; then
    echo
    echo "[yaml] Converting CSV -> suite YAML ..."
    python project/scripts/csv_to_slurm_yaml.py \
        "$CSV_PATH" \
        --base   "$BASE" \
        --output "$YAML_PATH"
fi

# ── step: slurm ─────────────────────────────────────────────────────────────
if has_step "slurm"; then
    echo
    echo "[slurm] Generating SLURM scripts ..."
    mkdir -p "$SLURM_DIR"
    python "$GENERATE_SLURM" \
        "$YAML_PATH" \
        --machine "$MACHINE_DIR" \
        --queue   "$QUEUE" \
        --output-dir "$SLURM_DIR"
fi

# ── summary ─────────────────────────────────────────────────────────────────
echo
echo "──────────────────────────────────────────────────────────────"
echo "Done. Summary:"
has_step "sample" && echo "  CSV       : $CSV_PATH"
has_step "plot"   && echo "  Figures   : project/output/${STEM}_{selection,corner,marginals}.png"
has_step "yaml"   && echo "  Suite YAML: $YAML_PATH"
if has_step "slurm"; then
    N_SLURM=$(ls -1 "$SLURM_DIR"/*.slurm 2>/dev/null | wc -l | tr -d ' ')
    echo "  SLURM dir : $SLURM_DIR   ($N_SLURM scripts)"
fi
echo "──────────────────────────────────────────────────────────────"
