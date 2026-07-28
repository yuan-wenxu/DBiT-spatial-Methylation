#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: dbitm.sh <assay> <step> --input <path> [--config <path>]

Main control script for DBiT-spatial-Methylation pipeline.

Arguments:
  assay          Assay type: taps | taps-v2 | emseq
  step           Pipeline step: fastp | barcode | align | spike-align | pool | call
  --input PATH   Raw FASTQ directory path
  --config PATH  Optional config file (default: config/dbitm.config.sh)

Execution mode is controlled by RUN_MODE in the config file:
  RUN_MODE=local   Run step directly (default)
  RUN_MODE=hpc     Submit step via sbatch

Examples:
  dbitm.sh taps fastp --input /data/raw
  dbitm.sh emseq align --input /data/raw --config my_config.sh
  dbitm.sh taps spike-align --input /data/raw --config my_config.sh
EOF
    exit 1
}

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
STEPS_DIR="$SCRIPT_DIR/steps"

declare -A STEP_SCRIPTS=(
    [fastp]=01.fastp.sh
    [barcode]=02.barcode.sh
    [align]=03.align.sh
    [spike-align]=03.spike_align.sh
    [pool]=04.pool.sh
    [call]=05.call.sh
)

# ── parse positional arguments ──
if [[ $# -lt 2 ]]; then
    usage
fi
assay=$1
shift
step=$1
shift

case "$assay" in
    taps|taps-v2|emseq) ;;
    *) echo "[dbitm] error: unsupported assay: $assay" >&2; usage ;;
esac

if [[ -z "${STEP_SCRIPTS[$step]:-}" ]]; then
    echo "[dbitm] error: unsupported step: $step" >&2
    echo "[dbitm] available steps: ${!STEP_SCRIPTS[*]}" >&2
    usage
fi

# ── parse optional arguments ──
input=""
config=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            input=$2
            shift 2
            ;;
        --config)
            config=$2
            shift 2
            ;;
        *)
            echo "[dbitm] error: unknown option: $1" >&2
            usage
            ;;
    esac
done

if [[ -z "$input" ]]; then
    echo "[dbitm] error: --input is required" >&2
    usage
fi
if [[ ! -d "$input" ]]; then
    echo "[dbitm] error: input directory not found: $input" >&2
    exit 1
fi
input=$(realpath "$input")

# resolve config
DEFAULT_CONFIG="$REPO_DIR/config/dbitm.config.sh"
if [[ -z "$config" ]]; then
    config=$DEFAULT_CONFIG
fi
if [[ "$config" != /* ]]; then
    config=$(realpath "$config")
fi
if [[ ! -f "$config" ]]; then
    echo "[dbitm] error: config file not found: $config" >&2
    exit 1
fi

if [[ "$config" != "$DEFAULT_CONFIG" ]]; then
    echo "[dbitm] using custom config: $config"
fi

export DBITM_PROJECT_ROOT=$REPO_DIR
export DBITM_CONFIG=$config

STEP_SCRIPT="$STEPS_DIR/${STEP_SCRIPTS[$step]}"
echo "[dbitm] assay: $assay | step: $step | input: $input"
echo "[dbitm] config: $config"

source "$config"
RUN_MODE=${RUN_MODE:-local}
case "$RUN_MODE" in
    local|hpc) ;;
    *) echo "[dbitm] error: RUN_MODE must be local or hpc, got: $RUN_MODE" >&2; exit 1 ;;
esac
echo "[dbitm] run mode: $RUN_MODE"

if [[ "$RUN_MODE" == hpc ]]; then

    step_upper=$(echo "$step" | tr '[:lower:]-' '[:upper:]_')
    job_name_var="${step_upper}_NAME"
    threads_var="${step_upper}_THREADS"
    partition_var="${step_upper}_PARTITION"
    mem_var="${step_upper}_MEM"
    time_var="${step_upper}_TIME"

    required_resource_vars=(
        "$job_name_var"
        "$threads_var"
        "$mem_var"
        "$time_var"
    )
    for resource_var in "${required_resource_vars[@]}"; do
        if [[ -z ${!resource_var:-} ]]; then
            echo "[dbitm] error: $resource_var is required for '$step' in hpc mode" >&2
            exit 1
        fi
    done

    job_name=${!job_name_var}
    threads=${!threads_var}
    partition=${!partition_var:-}
    mem=${!mem_var}
    time=${!time_var}

    sbatch_args=(
        --job-name="$job_name"
        --cpus-per-task="$threads"
        --mem="$mem"
        --time="$time"
        --output="${SBATCH_OUTPUT:-%x_%j.out}"
        --error="${SBATCH_ERROR:-%x_%j.err}"
    )
    [[ -n "$partition" ]] && sbatch_args+=(--partition="$partition")
    [[ "${SBATCH_REQUEUE:-}" == true ]] && sbatch_args+=(--requeue)

    echo "[dbitm] submitting $step via sbatch (job-name=$job_name cpus=$threads mem=$mem time=$time)..."
    sbatch "${sbatch_args[@]}" \
        --wrap="export DBITM_PROJECT_ROOT=$REPO_DIR DBITM_CONFIG=$config; '$STEP_SCRIPT' '$assay' '$input'"
else
    echo "[dbitm] running $step directly..."
    "$STEP_SCRIPT" "$assay" "$input"
fi
