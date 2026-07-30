#!/usr/bin/env bash
set -euo pipefail

usage() {
    local status=${1:-1}
    cat <<'EOF'
Usage: dbitm.sh <assay> <step> --input <path> [--config <path>]

Main control script for DBiT-spatial-Methylation pipeline.

Arguments:
  assay          Assay type: taps | taps-v2 | emseq
  step           Pipeline step: fastp | barcode | align | spike-align | pool | call | all
  --input PATH   Raw FASTQ directory path
  --config PATH  Optional config file (default: config/dbitm.config.sh)
  -h, --help     Show this help message and exit

Execution mode is controlled by RUN_MODE in the config file:
  RUN_MODE=local   Run step directly (default)
  RUN_MODE=hpc     Submit step via sbatch

The all step runs/submits:
  fastp -> barcode -> (align + spike-align) -> pool -> call

Examples:
  dbitm.sh taps fastp --input /data/raw
  dbitm.sh emseq align --input /data/raw --config my_config.sh
  dbitm.sh taps spike-align --input /data/raw --config my_config.sh
  dbitm.sh emseq all --input /data/raw --config my_config.sh
EOF
    exit "$status"
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
ALL_STEPS=(fastp barcode align spike-align pool call)

# ── parse positional arguments ──
case "${1:-}" in
    -h|--help) usage 0 ;;
esac
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

if [[ "$step" != all && -z "${STEP_SCRIPTS[$step]:-}" ]]; then
    echo "[dbitm] error: unsupported step: $step" >&2
    echo "[dbitm] available steps: ${ALL_STEPS[*]} all" >&2
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
        -h|--help)
            usage 0
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

echo "[dbitm] assay: $assay | step: $step | input: $input"
echo "[dbitm] config: $config"

source "$config"

case "$RUN_MODE" in
    local|hpc) ;;
    *) echo "[dbitm] error: RUN_MODE must be local or hpc, got: $RUN_MODE" >&2; exit 1 ;;
esac
echo "[dbitm] run mode: $RUN_MODE"

get_step_script() {
    local step_name=$1
    printf '%s/%s\n' "$STEPS_DIR" "${STEP_SCRIPTS[$step_name]}"
}

run_step_local() {
    local step_name=$1
    local step_script
    step_script=$(get_step_script "$step_name")
    echo "[dbitm] running $step_name directly..."
    "$step_script" "$assay" "$input"
}

SUBMITTED_JOB_ID=

validate_hpc_resources() {
    local step_name=$1
    local step_upper
    local resource_suffix resource_var

    step_upper=$(echo "$step_name" | tr '[:lower:]-' '[:upper:]_')
    for resource_suffix in NAME THREADS MEM TIME; do
        resource_var="${step_upper}_${resource_suffix}"
        if [[ -z ${!resource_var:-} ]]; then
            echo "[dbitm] error: $resource_var is required for '$step_name' in hpc mode" >&2
            return 1
        fi
    done
}

submit_step() {
    local step_name=$1
    local dependency=${2:-}
    local step_upper
    local job_name_var threads_var partition_var mem_var time_var
    local job_name threads partition mem time
    local step_script wrapped_command submission job_id
    local -a sbatch_args

    validate_hpc_resources "$step_name"
    step_upper=$(echo "$step_name" | tr '[:lower:]-' '[:upper:]_')
    job_name_var="${step_upper}_NAME"
    threads_var="${step_upper}_THREADS"
    partition_var="${step_upper}_PARTITION"
    mem_var="${step_upper}_MEM"
    time_var="${step_upper}_TIME"

    job_name=${!job_name_var}
    threads=${!threads_var}
    partition=${!partition_var:-}
    mem=${!mem_var}
    time=${!time_var}
    step_script=$(get_step_script "$step_name")

    sbatch_args=(
        --parsable
        --job-name="$job_name"
        --cpus-per-task="$threads"
        --mem="$mem"
        --time="$time"
        --output="${SBATCH_OUTPUT:-%x_%j.out}"
        --error="${SBATCH_ERROR:-%x_%j.err}"
    )
    [[ -n "$partition" ]] && sbatch_args+=(--partition="$partition")
    [[ "${SBATCH_REQUEUE:-}" == true ]] && sbatch_args+=(--requeue)
    [[ -n "$dependency" ]] && sbatch_args+=(--dependency="afterok:$dependency")

    printf -v wrapped_command \
        'export DBITM_PROJECT_ROOT=%q DBITM_CONFIG=%q; %q %q %q' \
        "$REPO_DIR" "$config" "$step_script" "$assay" "$input"

    echo "[dbitm] submitting $step_name via sbatch (job-name=$job_name cpus=$threads mem=$mem time=$time dependency=${dependency:-none})..."
    submission=$(sbatch "${sbatch_args[@]}" --wrap="$wrapped_command")
    job_id=${submission%%;*}
    if [[ ! "$job_id" =~ ^[0-9]+$ ]]; then
        echo "[dbitm] error: unable to parse sbatch job ID for '$step_name': $submission" >&2
        return 1
    fi
    SUBMITTED_JOB_ID=$job_id
    echo "[dbitm] submitted $step_name job-id=$job_id"
}

run_all_local() {
    local step_name
    echo "[dbitm] running complete pipeline locally..."
    for step_name in "${ALL_STEPS[@]}"; do
        run_step_local "$step_name"
    done
    echo "[dbitm] complete pipeline finished successfully"
}

submit_all_hpc() {
    local fastp_job_id barcode_job_id align_job_id spike_align_job_id
    local pool_job_id call_job_id step_name

    echo "[dbitm] submitting complete pipeline with Slurm dependencies..."
    for step_name in "${ALL_STEPS[@]}"; do
        validate_hpc_resources "$step_name"
    done

    submit_step fastp
    fastp_job_id=$SUBMITTED_JOB_ID

    submit_step barcode "$fastp_job_id"
    barcode_job_id=$SUBMITTED_JOB_ID

    submit_step align "$barcode_job_id"
    align_job_id=$SUBMITTED_JOB_ID

    submit_step spike-align "$barcode_job_id"
    spike_align_job_id=$SUBMITTED_JOB_ID

    submit_step pool "$align_job_id:$spike_align_job_id"
    pool_job_id=$SUBMITTED_JOB_ID

    submit_step call "$pool_job_id"
    call_job_id=$SUBMITTED_JOB_ID

    echo "[dbitm] complete pipeline submitted successfully"
    echo "[dbitm] job IDs: fastp=$fastp_job_id barcode=$barcode_job_id align=$align_job_id spike-align=$spike_align_job_id pool=$pool_job_id call=$call_job_id"
}

if [[ "$RUN_MODE" == hpc ]]; then
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "[dbitm] error: sbatch executable not found" >&2
        exit 1
    fi
    if [[ "$step" == all ]]; then
        submit_all_hpc
    else
        submit_step "$step"
    fi
else
    if [[ "$step" == all ]]; then
        run_all_local
    else
        run_step_local "$step"
    fi
fi
