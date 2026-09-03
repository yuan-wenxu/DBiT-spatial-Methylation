#!/usr/bin/env bash
set -euo pipefail

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ ${use_scratch:-false} == true && -n ${scratch_run:-} && -d $scratch_run ]]; then
        if (( status != 0 )) && [[ -n ${run_output:-} && -d $run_output && -n ${output_dir:-} ]]; then
            echo "[dbitm] saturation: recovering scratch results after exit status $status: $output_dir" >&2
            if mkdir -p "$output_dir" && cp -a "$run_output/." "$output_dir/"; then
                echo "[dbitm] saturation: scratch results recovered: $output_dir" >&2
            else
                echo "[dbitm] saturation: warning: failed to recover scratch results: $run_output" >&2
            fi
        fi
        rm -rf -- "$scratch_run" || echo "[dbitm] saturation: warning: failed to clean scratch directory: $scratch_run" >&2
    fi
    exit "$status"
}

enable_cleanup() {
    trap cleanup_scratch EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM
    trap 'exit 129' HUP
}

if (( $# < 2 || $# > 3 )); then
    echo "Usage: 07.saturation.sh <assay> <raw_fastq_folder> [--dry-run]" >&2
    exit 1
fi
assay=$1
raw_path=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] saturation: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
case "$assay" in
    taps|taps-v2|emseq|cabernet|smc) ;;
    *) echo "[dbitm] saturation: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] saturation: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] saturation: config file not found: $config_file" >&2
    exit 1
fi
source "$config_file"

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
output_dir=$final_dir/saturation
saturation_script=$REPO_DIR/script/steps/python/07.saturation.py
if [[ ! -f "$saturation_script" ]]; then
    echo "[dbitm] saturation: Python script not found: $saturation_script" >&2
    exit 1
fi

declare -a saturation_args=(
    --work-dir "$final_dir"
    --reads-threshold "$SATURATION_READS_THRESHOLD"
    --pred-fraction "$SATURATION_PRED_FRACTION"
    --linear-r2-threshold "$SATURATION_LINEAR_R2_THRESHOLD"
)

echo "====== dbitm saturation ======"
echo "[dbitm] work directory: $final_dir"
echo "[dbitm] output directory: $output_dir"
echo "[dbitm] reads threshold: $SATURATION_READS_THRESHOLD"
echo "[dbitm] prediction fraction: $SATURATION_PRED_FRACTION"
echo "[dbitm] linear R2 threshold: $SATURATION_LINEAR_R2_THRESHOLD"
echo "[dbitm] config: $config_file"

if [[ "$dry_run" == true ]]; then
    if [[ -n ${SCRATCH_ROOT:-} && "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] saturation: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        python "$saturation_script" "${saturation_args[@]}" \
        --output-dir "$output_dir" \
        --dry-run
    [[ -z ${SCRATCH_ROOT:-} ]] || echo "[dbitm] planned scratch root: $SCRATCH_ROOT"
    echo "[dbitm] dry-run: no files will be written"
    echo "====== dbitm saturation dry-run finished ======"
    exit 0
fi

use_scratch=false
run_output=$output_dir
if [[ -n ${SCRATCH_ROOT:-} ]]; then
    if [[ "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] saturation: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    scratch_root=$(realpath -m "$SCRATCH_ROOT")
    run_id=${SLURM_JOB_ID:-saturation_$$}
    scratch_run=$scratch_root/dbitm/$run_id
    run_output=$scratch_run/saturation
    use_scratch=true
    enable_cleanup
fi
mkdir -p "$run_output"

pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
    python "$saturation_script" "${saturation_args[@]}" \
    --output-dir "$run_output" \
    2>&1 | tee "$run_output/saturation.log"

if [[ "$use_scratch" == true ]]; then
    mkdir -p "$output_dir"
    cp -a "$run_output/." "$output_dir/"
fi
echo "[dbitm] saturation result: $output_dir"
echo "====== dbitm saturation finished ======"
