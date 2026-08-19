#!/usr/bin/env bash
set -euo pipefail

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ ${use_scratch:-false} == true && -n ${scratch_run:-} && -d $scratch_run ]]; then
        if (( status != 0 )) && [[ -n ${run_output:-} && -d $run_output && -n ${output_dir:-} ]]; then
            echo "[dbitm] methscan: recovering scratch results after exit status $status: $output_dir" >&2
            if mkdir -p "$output_dir" && cp -a "$run_output/." "$output_dir/"; then
                echo "[dbitm] methscan: scratch results recovered: $output_dir" >&2
            else
                echo "[dbitm] methscan: warning: failed to recover scratch results: $run_output" >&2
            fi
        fi
        rm -rf -- "$scratch_run" || echo "[dbitm] methscan: warning: failed to clean scratch directory: $scratch_run" >&2
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
    echo "Usage: 09.methscan.sh <assay> <raw_fastq_folder> [--dry-run]" >&2
    exit 1
fi
assay=$1
raw_path=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] methscan: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
case "$assay" in
    taps|taps-v2|emseq|cabernet) ;;
    *) echo "[dbitm] methscan: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] methscan: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] methscan: config file not found: $config_file" >&2
    exit 1
fi
source "$config_file"

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
coverage_dir=$final_dir/coverage/host
output_dir=$final_dir/methscan

context_mode=${CALL_CONTEXT_MODE}
case "$context_mode" in
    cg) methscan_contexts=(CG) ;;
    ch) methscan_contexts=(CA CC CT) ;;
    both) methscan_contexts=(CG CA CC CT) ;;
    *) echo "[dbitm] methscan: CALL_CONTEXT_MODE must be cg, ch, or both" >&2; exit 1 ;;
esac

methscan_chunksize=${METHSCAN_CHUNKSIZE}
declare -A methscan_min_sites_by_context=(
    [CG]="${METHSCAN_CG_MIN_SITES}"
    [CA]="${METHSCAN_CA_MIN_SITES}"
    [CC]="${METHSCAN_CC_MIN_SITES}"
    [CT]="${METHSCAN_CT_MIN_SITES}"
)
methscan_threads=${METHSCAN_THREADS}
if [[ ! "$methscan_chunksize" =~ ^[1-9][0-9]*$ ]]; then
    echo "[dbitm] methscan: METHSCAN_CHUNKSIZE must be greater than zero" >&2
    exit 1
fi
for methscan_context in "${methscan_contexts[@]}"; do
    methscan_min_sites=${methscan_min_sites_by_context[$methscan_context]}
    if [[ ! "$methscan_min_sites" =~ ^[1-9][0-9]*$ ]]; then
        echo "[dbitm] methscan: METHSCAN_${methscan_context}_MIN_SITES must be greater than zero" >&2
        exit 1
    fi
done
if [[ ! "$methscan_threads" =~ ^[1-9][0-9]*$ ]]; then
    echo "[dbitm] methscan: METHSCAN_THREADS must be greater than zero" >&2
    exit 1
fi

echo "====== dbitm MethSCAn matrices ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] contexts: ${methscan_contexts[*]}"
echo "[dbitm] coverage directory: $coverage_dir"
echo "[dbitm] output directory: $output_dir"
for methscan_context in "${methscan_contexts[@]}"; do
    echo "[dbitm] $methscan_context min sites: ${methscan_min_sites_by_context[$methscan_context]}"
done
echo "[dbitm] threads: $methscan_threads"
echo "[dbitm] config: $config_file"

if [[ "$dry_run" == true ]]; then
    if [[ -n ${SCRATCH_ROOT:-} && "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] methscan: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    for methscan_context in "${methscan_contexts[@]}"; do
        context_output=$output_dir/$methscan_context
        echo "[dbitm] planned command: methscan prepare $coverage_dir/**/*.${methscan_context}.cov -> $context_output/compact"
        echo "[dbitm] planned command: methscan filter $context_output/compact -> $context_output/filter"
        echo "[dbitm] planned command: methscan smooth $context_output/filter"
        echo "[dbitm] planned command: methscan scan $context_output/filter -> $context_output/VMRs.bed"
        echo "[dbitm] planned command: methscan matrix $context_output/VMRs.bed -> $context_output/matrix"
    done
    [[ -z ${SCRATCH_ROOT:-} ]] || echo "[dbitm] planned scratch root: $SCRATCH_ROOT"
    echo "[dbitm] dry-run: no files will be written"
    echo "====== dbitm MethSCAn matrices dry-run finished ======"
    exit 0
fi

if [[ ! -d "$coverage_dir" ]]; then
    echo "[dbitm] methscan: host coverage directory not found: $coverage_dir" >&2
    echo "[dbitm] methscan: run the call step first." >&2
    exit 1
fi
for methscan_context in "${methscan_contexts[@]}"; do
    if ! find "$coverage_dir" -type f -name "*.${methscan_context}.cov" -print -quit | grep -q .; then
        echo "[dbitm] methscan: no host *.${methscan_context}.cov files found under $coverage_dir" >&2
        exit 1
    fi
done

use_scratch=false
run_output=$output_dir
if [[ -n ${SCRATCH_ROOT:-} ]]; then
    if [[ "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] methscan: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    scratch_root=$(realpath -m "$SCRATCH_ROOT")
    run_id=${SLURM_JOB_ID:-methscan_$$}
    scratch_run=$scratch_root/dbitm/$run_id
    run_output=$scratch_run/methscan
    use_scratch=true
    enable_cleanup
fi

mkdir -p "$run_output"
for methscan_context in "${methscan_contexts[@]}"; do
    context_run_output=$run_output/$methscan_context
    compact_dir=$context_run_output/compact
    filtered_dir=$context_run_output/filter
    vmrs_bed=$context_run_output/VMRs.bed
    matrix_dir=$context_run_output/matrix
    methscan_log=$context_run_output/methscan.log
    methscan_min_sites=${methscan_min_sites_by_context[$methscan_context]}
    mapfile -d '' -t context_cov_files < <(
        find "$coverage_dir" -type f -name "*.${methscan_context}.cov" -print0 | sort -z
    )
    mkdir -p "$context_run_output"
    {
        echo "[dbitm] MethSCAn context: $methscan_context"
        echo "[dbitm] MethSCAn inputs: ${#context_cov_files[@]}"
        echo "[dbitm] MethSCAn min sites: $methscan_min_sites"
        echo "[dbitm] MethSCAn run directory: $context_run_output"
    } | tee "$methscan_log"

    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        methscan prepare "${context_cov_files[@]}" "$compact_dir" \
        --input-format bismark \
        --chunksize "$methscan_chunksize" \
        2>&1 | tee -a "$methscan_log"
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        methscan filter "$compact_dir" "$filtered_dir" \
        --min-sites "$methscan_min_sites" \
        2>&1 | tee -a "$methscan_log"
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        methscan smooth "$filtered_dir" \
        2>&1 | tee -a "$methscan_log"
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        methscan scan "$filtered_dir" "$vmrs_bed" \
        --threads "$methscan_threads" \
        2>&1 | tee -a "$methscan_log"
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        methscan matrix "$vmrs_bed" "$filtered_dir" "$matrix_dir" \
        --threads "$methscan_threads" \
        2>&1 | tee -a "$methscan_log"
done

if [[ "$use_scratch" == true ]]; then
    mkdir -p "$output_dir"
    cp -a "$run_output/." "$output_dir/"
fi
for methscan_context in "${methscan_contexts[@]}"; do
    echo "[dbitm] MethSCAn $methscan_context matrix result: $output_dir/$methscan_context/matrix"
done
echo "====== dbitm MethSCAn matrices finished ======"
