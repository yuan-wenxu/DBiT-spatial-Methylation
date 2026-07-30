#!/usr/bin/env bash
set -euo pipefail

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ ${use_scratch:-false} == true && -n ${scratch_run:-} && -d $scratch_run ]]; then
        rm -rf -- "$scratch_run" || echo "[dbitm] call: warning: failed to clean scratch directory: $scratch_run" >&2
    fi
    exit "$status"
}

enable_cleanup() {
    trap cleanup_scratch EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM
    trap 'exit 129' HUP
}

if [[ $# -ne 2 ]]; then
    echo "Usage: 05.call.sh <assay> <raw_fastq_folder>" >&2
    exit 1
fi
assay=$1
raw_path=$2
case "$assay" in
    taps|taps-v2|emseq) ;;
    *) echo "[dbitm] call: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] call: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] call: config file not found: $config_file" >&2
    exit 1
fi
source "$config_file"

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
pooled_bam=$final_dir/pooled/pooled.byCB.bam
if [[ "$assay" == emseq ]]; then
    caller_script=$REPO_DIR/script/steps/05.methy_caller_emseq.sh
else
    caller_script=$REPO_DIR/script/steps/python/05.methy_caller_taps.py
fi

if [[ ! -f "$pooled_bam" ]]; then
    echo "[dbitm] call: pooled BAM not found: $pooled_bam" >&2
    echo "[dbitm] call: run the pool step first." >&2
    exit 1
fi
if [[ ! -f "$caller_script" ]]; then
    echo "[dbitm] call: caller script not found: $caller_script" >&2
    exit 1
fi

# resolve reference
reference=${CALL_REFERENCE:-}
if [[ -z "$reference" ]]; then
    echo "[dbitm] call: CALL_REFERENCE is required" >&2
    exit 1
fi
if [[ "$reference" != /* ]]; then
    reference=$REPO_DIR/$reference
fi
reference=$(realpath -m "$reference")

use_scratch=false
run_bam=$pooled_bam
run_output=$final_dir/coverage

echo "====== dbitm call ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] reference: $reference"
echo "[dbitm] pooled BAM: $pooled_bam"
echo "[dbitm] chromosomes: $CALL_CHROMOSOMES"
caller_context_mode=$CALL_CONTEXT_MODE
echo "[dbitm] context mode: $caller_context_mode"
echo "[dbitm] config: $config_file"
echo "[dbitm] output directory: $final_dir/coverage"

if [[ -n ${SCRATCH_ROOT:-} ]]; then
    if [[ "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] call: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    echo "[dbitm] scratch root: $SCRATCH_ROOT"
    mkdir -p "$SCRATCH_ROOT"
    scratch_root=$(realpath "$SCRATCH_ROOT")
    run_id=${SLURM_JOB_ID:-call_$$}
    scratch_run=$scratch_root/dbitm/$run_id
    scratch_input=$scratch_run/pooled_input
    run_bam=$scratch_input/pooled.byCB.bam
    run_output=$scratch_run/coverage
    use_scratch=true
    enable_cleanup
    mkdir -p "$scratch_input"
    echo "[dbitm] copying pooled BAM + index to scratch: $scratch_input"
    cp -a "$pooled_bam" "$scratch_input/"
    cp -a "$pooled_bam.bai" "$scratch_input/" 2>/dev/null || true
    cp -a "$pooled_bam.csi" "$scratch_input/" 2>/dev/null || true
fi

mkdir -p "$run_output"

echo "[dbitm] starting methylation calling..."
if [[ "$assay" == emseq ]]; then
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        bash "$caller_script" \
        --bam "$run_bam" \
        --reference "$reference" \
        --out-dir "$run_output" \
        --chromosomes "$CALL_CHROMOSOMES" \
        --context-mode "$caller_context_mode" \
        --min-base-quality "$CALL_MIN_BASE_QUALITY" \
        --min-mapping-quality "$CALL_MIN_MAPPING_QUALITY" \
        --r1-left-trim "$CALL_R1_LEFT_TRIMMING" \
        --r1-right-trim "$CALL_R1_RIGHT_TRIMMING" \
        --r2-left-trim "$CALL_R2_LEFT_TRIMMING" \
        --r2-right-trim "$CALL_R2_RIGHT_TRIMMING" \
        --jobs "$CALL_JOBS" \
        --max-spots "${EMSEQ_CALL_MAX_SPOTS:-10000}" \
        > "$run_output/call.log" 2>&1
else
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        python "$caller_script" \
        --bam "$run_bam" \
        --reference "$reference" \
        --out-dir "$run_output" \
        --chromosomes "$CALL_CHROMOSOMES" \
        --context-mode "$caller_context_mode" \
        --min-base-quality "$CALL_MIN_BASE_QUALITY" \
        --min-mapping-quality "$CALL_MIN_MAPPING_QUALITY" \
        --max-depth "$CALL_MAX_DEPTH" \
        --batch-size "$CALL_BATCH_SIZE" \
        --r1-left-trim "$CALL_R1_LEFT_TRIMMING" \
        --r1-right-trim "$CALL_R1_RIGHT_TRIMMING" \
        --r2-left-trim "$CALL_R2_LEFT_TRIMMING" \
        --r2-right-trim "$CALL_R2_RIGHT_TRIMMING" \
        --jobs "$CALL_JOBS" \
        > "$run_output/call.log" 2>&1
fi

echo "[dbitm] call finished successfully"
if [[ "$use_scratch" == true ]]; then
    rm -rf -- "$final_dir/coverage"
    mkdir -p "$final_dir"
    echo "[dbitm] copying coverage result from scratch to $final_dir/coverage"
    cp -a "$run_output" "$final_dir/coverage"
fi
echo "[dbitm] call log: $final_dir/coverage/call.log"
echo "[dbitm] call result: $final_dir/coverage"
echo "====== dbitm call finished ======"
