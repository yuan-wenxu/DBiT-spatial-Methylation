#!/usr/bin/env bash
set -euo pipefail

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ ${use_scratch:-false} == true && -n ${scratch_run:-} && -d $scratch_run ]]; then
        if (( status != 0 )) && [[ -n ${run_output:-} && -d $run_output && -n ${output_dir:-} ]]; then
            echo "[dbitm] summary: recovering scratch results after exit status $status: $output_dir" >&2
            if mkdir -p "$output_dir" && cp -a "$run_output/." "$output_dir/"; then
                echo "[dbitm] summary: scratch results recovered: $output_dir" >&2
            else
                echo "[dbitm] summary: warning: failed to recover scratch results: $run_output" >&2
            fi
        fi
        rm -rf -- "$scratch_run" || echo "[dbitm] summary: warning: failed to clean scratch directory: $scratch_run" >&2
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
    echo "Usage: 08.summary.sh <assay> <raw_fastq_folder> [--dry-run]" >&2
    exit 1
fi
assay=$1
raw_path=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] summary: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
case "$assay" in
    taps|taps-v2|emseq|cabernet) ;;
    *) echo "[dbitm] summary: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] summary: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] summary: config file not found: $config_file" >&2
    exit 1
fi
source "$config_file"

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
output_dir=$final_dir/summary
summary_script=$REPO_DIR/script/steps/python/08.summary.py
if [[ ! -f "$summary_script" ]]; then
    echo "[dbitm] summary: Python script not found: $summary_script" >&2
    exit 1
fi

minimum_mapq=${CALL_MIN_MAPPING_QUALITY}
if [[ ! "$minimum_mapq" =~ ^[0-9]+$ ]]; then
    echo "[dbitm] summary: CALL_MIN_MAPPING_QUALITY must be >= 0" >&2
    exit 1
fi

context_mode=${CALL_CONTEXT_MODE}
case "$context_mode" in
    cg|ch|both) ;;
    *) echo "[dbitm] summary: CALL_CONTEXT_MODE must be cg, ch, or both" >&2; exit 1 ;;
esac

declare -a spike_names=()
spike_declaration=$(declare -p CALL_SPIKE_IN_REFERENCES 2>/dev/null || true)
if [[ -n "$spike_declaration" ]]; then
    if [[ "$spike_declaration" != "declare -A "* ]]; then
        echo "[dbitm] summary: CALL_SPIKE_IN_REFERENCES must be an associative array" >&2
        exit 1
    fi
    declare -n spike_references=CALL_SPIKE_IN_REFERENCES
    if (( ${#spike_references[@]} > 0 )); then
        while IFS= read -r spike_name; do
            [[ -n ${spike_references[$spike_name]} ]] || continue
            if [[ ! "$spike_name" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
                echo "[dbitm] summary: invalid spike-in name: $spike_name" >&2
                exit 1
            fi
            spike_names+=("$spike_name")
        done < <(printf '%s\n' "${!spike_references[@]}" | LC_ALL=C sort)
    fi
fi

declare -a summary_args=(
    --work-dir "$final_dir"
    --min-mapping-quality "$minimum_mapq"
    --context-mode "$context_mode"
)
for spike_name in "${spike_names[@]}"; do
    summary_args+=(--spike-in-name "$spike_name")
done

echo "====== dbitm summary ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] work directory: $final_dir"
echo "[dbitm] output directory: $output_dir"
echo "[dbitm] context mode: $context_mode"
echo "[dbitm] spike-ins: ${spike_names[*]:-none}"
echo "[dbitm] config: $config_file"

if [[ "$dry_run" == true ]]; then
    if [[ -n ${SCRATCH_ROOT:-} && "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] summary: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        python "$summary_script" "${summary_args[@]}" \
        --output-dir "$output_dir" \
        --dry-run
    [[ -z ${SCRATCH_ROOT:-} ]] || echo "[dbitm] planned scratch root: $SCRATCH_ROOT"
    echo "[dbitm] dry-run: no files will be written"
    echo "====== dbitm summary dry-run finished ======"
    exit 0
fi

use_scratch=false
run_output=$output_dir
if [[ -n ${SCRATCH_ROOT:-} ]]; then
    if [[ "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] summary: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    scratch_root=$(realpath -m "$SCRATCH_ROOT")
    run_id=${SLURM_JOB_ID:-summary_$$}
    scratch_run=$scratch_root/dbitm/$run_id
    run_output=$scratch_run/summary
    use_scratch=true
    enable_cleanup
fi
mkdir -p "$run_output"

pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
    python "$summary_script" "${summary_args[@]}" \
    --output-dir "$run_output" \
    2>&1 | tee "$run_output/summary.log"

if [[ "$use_scratch" == true ]]; then
    mkdir -p "$output_dir"
    cp -a "$run_output/." "$output_dir/"
fi
echo "[dbitm] summary result: $output_dir"
echo "====== dbitm summary finished ======"
