#!/usr/bin/env bash
set -euo pipefail

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ ${use_scratch:-false} == true && -n ${scratch_run:-} && -d $scratch_run ]]; then
        if (( status != 0 )) && [[ -n ${run_output:-} && -d $run_output && -n ${final_dir:-} ]]; then
            echo "[dbitm] call: recovering scratch results after exit status $status: $final_dir/coverage" >&2
            if mkdir -p "$final_dir/coverage" && cp -a "$run_output/." "$final_dir/coverage/"; then
                echo "[dbitm] call: scratch results recovered: $final_dir/coverage" >&2
            else
                echo "[dbitm] call: warning: failed to recover scratch results: $run_output" >&2
            fi
        fi
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

load_trims() {
    local cutoff_path=$1
    awk -F '\t' '
        FNR == 1 { next }
        $1 == "R1" { r1_left=$7; r1_right=$8; seen_r1=1 }
        $1 == "R2" { r2_left=$7; r2_right=$8; seen_r2=1 }
        END {
            if (!seen_r1 || !seen_r2 ||
                r1_left !~ /^[0-9]+$/ || r1_right !~ /^[0-9]+$/ ||
                r2_left !~ /^[0-9]+$/ || r2_right !~ /^[0-9]+$/) exit 2
            printf "%d\t%d\t%d\t%d\n", r1_left, r1_right, r2_left, r2_right
        }
    ' "$cutoff_path"
}

if (( $# < 2 || $# > 3 )); then
    echo "Usage: 06.call.sh <assay> <raw_fastq_folder> [--dry-run]" >&2
    exit 1
fi
assay=$1
raw_path=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] call: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
case "$assay" in
    taps|taps-v2|emseq|cabernet) ;;
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

caller_context_mode=${CALL_CONTEXT_MODE:-both}
case "$caller_context_mode" in
    cg|ch|both) ;;
    *) echo "[dbitm] call: CALL_CONTEXT_MODE must be cg, ch, or both" >&2; exit 1 ;;
esac

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
pooled_bam=$final_dir/pooled/pooled.cb.bam
caller_script=$REPO_DIR/script/steps/python/06.methy_caller.py

if [[ ! -f "$pooled_bam" ]]; then
    if [[ "$dry_run" == true ]]; then
        echo "[dbitm] call: dry-run expects future pooled BAM: $pooled_bam"
    else
        echo "[dbitm] call: pooled BAM not found: $pooled_bam" >&2
        echo "[dbitm] call: run the pool step first." >&2
        exit 1
    fi
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

barcode_whitelist=${BARCODE_WHITELIST:-$REPO_DIR/docs/barcodes/barcodes50.tsv}
if [[ "$barcode_whitelist" != /* ]]; then
    barcode_whitelist=$REPO_DIR/$barcode_whitelist
fi
if [[ ! -f "$barcode_whitelist" ]]; then
    echo "[dbitm] call: barcode whitelist not found: $barcode_whitelist" >&2
    exit 1
fi
barcode_whitelist=$(realpath "$barcode_whitelist")

cutoff_dir=$final_dir/mbias
cutoff_path=$cutoff_dir/host.mbias.cutoffs.tsv
cutoff_available=true
if [[ ! -f "$cutoff_path" ]]; then
    if [[ "$dry_run" == true ]]; then
        echo "[dbitm] call: dry-run expects future host M-bias cutoff: $cutoff_path"
        cutoff_available=false
    else
        echo "[dbitm] call: host M-bias cutoff file not found: $cutoff_path" >&2
        echo "[dbitm] call: run the mbias step first." >&2
        exit 1
    fi
fi

echo "====== dbitm call ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] reference: $reference"
echo "[dbitm] barcode whitelist: $barcode_whitelist"
echo "[dbitm] pooled BAM: $pooled_bam"
echo "[dbitm] chromosomes: $CALL_CHROMOSOMES"
echo "[dbitm] context mode: $caller_context_mode"
echo "[dbitm] host M-bias cutoff: $cutoff_path"
echo "[dbitm] config: $config_file"
echo "[dbitm] output directory: $final_dir/coverage"

if [[ "$cutoff_available" == false ]]; then
    call_r1_left_trim=auto
    call_r1_right_trim=auto
    call_r2_left_trim=auto
    call_r2_right_trim=auto
elif [[ ! -s "$cutoff_path" ]]; then
    echo "[dbitm] call: no suitable host M-bias cutoff; continuing without trimming"
    call_r1_left_trim=0
    call_r1_right_trim=0
    call_r2_left_trim=0
    call_r2_right_trim=0
else
    if ! trim_values=$(load_trims "$cutoff_path"); then
        echo "[dbitm] call: invalid host M-bias cutoff: $cutoff_path" >&2
        exit 1
    fi
    IFS=$'\t' read -r \
        call_r1_left_trim call_r1_right_trim \
        call_r2_left_trim call_r2_right_trim \
        <<< "$trim_values"
    for trim_name in \
        call_r1_left_trim call_r1_right_trim call_r2_left_trim call_r2_right_trim
    do
        trim_value=${!trim_name:-}
        if [[ ! "$trim_value" =~ ^[0-9]+$ ]]; then
            echo "[dbitm] call: invalid or missing $trim_name under: $cutoff_dir" >&2
            exit 1
        fi
    done
fi

echo "[dbitm] trimming: R1=${call_r1_left_trim},${call_r1_right_trim} R2=${call_r2_left_trim},${call_r2_right_trim}"

use_scratch=false
run_bam=$pooled_bam
run_output=$final_dir/coverage

if [[ "$dry_run" == true ]]; then
    if [[ -n ${SCRATCH_ROOT:-} && "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] call: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    echo "[dbitm] dry-run: no files will be written"
    echo "[dbitm] planned command: 06.methy_caller.py -> $final_dir/coverage"
    [[ -z ${SCRATCH_ROOT:-} ]] || echo "[dbitm] planned scratch root: $SCRATCH_ROOT"
    echo "====== dbitm call dry-run finished ======"
    exit 0
fi

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
    run_bam=$scratch_input/pooled.cb.bam
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
pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
    python "$caller_script" \
    --assay "$assay" \
    --bam "$run_bam" \
    --reference "$reference" \
    --out-dir "$run_output" \
    --barcode-whitelist "$barcode_whitelist" \
    --chromosomes "$CALL_CHROMOSOMES" \
    --context-mode "$caller_context_mode" \
    --min-base-quality "$CALL_MIN_BASE_QUALITY" \
    --min-mapping-quality "$CALL_MIN_MAPPING_QUALITY" \
    --max-depth "$CALL_MAX_DEPTH" \
    --batch-size "$CALL_BATCH_SIZE" \
    --r1-left-trim "$call_r1_left_trim" \
    --r1-right-trim "$call_r1_right_trim" \
    --r2-left-trim "$call_r2_left_trim" \
    --r2-right-trim "$call_r2_right_trim" \
    --jobs "$CALL_JOBS" \
    > "$run_output/call.log" 2>&1

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
