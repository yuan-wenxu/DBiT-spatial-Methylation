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

load_target_trims() {
    local label=$1
    local cutoff_path=$2
    local trim_values trim_name trim_value
    if [[ ! -s "$cutoff_path" ]]; then
        echo "[dbitm] call: no suitable host M-bias cutoff for $label; continuing without trimming"
        target_r1_left=0
        target_r1_right=0
        target_r2_left=0
        target_r2_right=0
        return 0
    fi
    if ! trim_values=$(load_trims "$cutoff_path"); then
        echo "[dbitm] call: invalid host M-bias cutoff ($label): $cutoff_path" >&2
        return 1
    fi
    IFS=$'\t' read -r \
        target_r1_left target_r1_right target_r2_left target_r2_right \
        <<< "$trim_values"
    for trim_name in target_r1_left target_r1_right target_r2_left target_r2_right; do
        trim_value=${!trim_name:-}
        if [[ ! "$trim_value" =~ ^[0-9]+$ ]]; then
            echo "[dbitm] call: invalid or missing trim for $label under: $cutoff_path" >&2
            return 1
        fi
    done
}

copy_bam_with_index() {
    local bam_path=$1
    local destination=$2
    local candidate
    cp -a "$bam_path" "$destination/"
    for candidate in "$bam_path.bai" "${bam_path%.bam}.bai" "$bam_path.csi"; do
        if [[ -f "$candidate" ]]; then
            cp -a "$candidate" "$destination/"
        fi
    done
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
    taps|taps-v2|emseq|cabernet|smc) ;;
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
    cg) call_contexts=(CG) ;;
    ch) call_contexts=(CA CC CT) ;;
    both) call_contexts=(CG CA CC CT) ;;
    *) echo "[dbitm] call: CALL_CONTEXT_MODE must be cg, ch, or both" >&2; exit 1 ;;
esac

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
pooled_dir=$final_dir/pooled
smc_filter_dir=$final_dir/smc_xxx_filter
cutoff_dir=$final_dir/mbias
caller_script=$REPO_DIR/script/steps/python/06.methy_caller.py
merge_script=$REPO_DIR/script/steps/python/06.merge_cov.py
if [[ ! -f "$caller_script" ]]; then
    echo "[dbitm] call: caller script not found: $caller_script" >&2
    exit 1
fi
if [[ "$assay" == smc && ! -f "$merge_script" ]]; then
    echo "[dbitm] call: coverage merge script not found: $merge_script" >&2
    exit 1
fi

declare -a call_labels=()
declare -a call_bams=()
declare -a cutoff_paths=()
if [[ "$assay" == smc ]]; then
    call_labels=(watson crick)
    call_bams=(
        "$smc_filter_dir/pooled.watson.cb.bam"
        "$smc_filter_dir/pooled.crick.cb.bam"
    )
    cutoff_paths=(
        "$cutoff_dir/host.watson.mbias.cutoffs.tsv"
        "$cutoff_dir/host.crick.mbias.cutoffs.tsv"
    )
else
    call_labels=(host)
    call_bams=("$pooled_dir/pooled.cb.bam")
    cutoff_paths=("$cutoff_dir/host.mbias.cutoffs.tsv")
fi

for call_index in "${!call_labels[@]}"; do
    label=${call_labels[$call_index]}
    bam_path=${call_bams[$call_index]}
    cutoff_path=${cutoff_paths[$call_index]}
    if [[ ! -f "$bam_path" ]]; then
        if [[ "$dry_run" == true ]]; then
            echo "[dbitm] call: dry-run expects future input BAM ($label): $bam_path"
        else
            echo "[dbitm] call: input BAM not found ($label): $bam_path" >&2
            if [[ "$assay" == smc ]]; then
                echo "[dbitm] call: run the smc-filter step first." >&2
            else
                echo "[dbitm] call: run the pool step first." >&2
            fi
            exit 1
        fi
    fi
    if [[ ! -f "$cutoff_path" ]]; then
        if [[ "$dry_run" == true ]]; then
            echo "[dbitm] call: dry-run expects future host M-bias cutoff ($label): $cutoff_path"
        else
            echo "[dbitm] call: host M-bias cutoff file not found ($label): $cutoff_path" >&2
            echo "[dbitm] call: run the mbias step first." >&2
            exit 1
        fi
    fi
done

reference=${CALL_REFERENCE:-}
if [[ -z "$reference" ]]; then
    echo "[dbitm] call: CALL_REFERENCE is required" >&2
    exit 1
fi
if [[ "$reference" != /* ]]; then
    reference=$REPO_DIR/$reference
fi
reference=$(realpath -m "$reference")

if [[ "$assay" == smc ]]; then
    barcode_whitelist=${SMC_BARCODE_WHITELIST:-$REPO_DIR/docs/barcodes/barcodes-smc.tsv}
else
    barcode_whitelist=${BARCODE_WHITELIST:-$REPO_DIR/docs/barcodes/barcodes50.tsv}
fi
if [[ "$barcode_whitelist" != /* ]]; then
    barcode_whitelist=$REPO_DIR/$barcode_whitelist
fi
if [[ ! -f "$barcode_whitelist" ]]; then
    echo "[dbitm] call: barcode whitelist not found: $barcode_whitelist" >&2
    exit 1
fi
barcode_whitelist=$(realpath "$barcode_whitelist")

echo "====== dbitm call ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] reference: $reference"
echo "[dbitm] barcode whitelist: $barcode_whitelist"
for call_index in "${!call_labels[@]}"; do
    echo "[dbitm] call target ${call_labels[$call_index]} BAM: ${call_bams[$call_index]}"
    echo "[dbitm] call target ${call_labels[$call_index]} M-bias cutoff: ${cutoff_paths[$call_index]}"
done
echo "[dbitm] chromosomes: $CALL_CHROMOSOMES"
echo "[dbitm] context mode: $caller_context_mode"
echo "[dbitm] config: $config_file"
echo "[dbitm] output directory: $final_dir/coverage"
if [[ "$assay" == smc ]]; then
    echo "[dbitm] SmC filtered BAM directory: $smc_filter_dir"
fi

if [[ "$dry_run" == true ]]; then
    if [[ -n ${SCRATCH_ROOT:-} && "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] call: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    for label in "${call_labels[@]}"; do
        echo "[dbitm] planned call target: $label"
    done
    if [[ "$assay" == smc ]]; then
        echo "[dbitm] planned merge: Watson + Crick calls -> $final_dir/coverage/host"
    fi
    echo "[dbitm] dry-run: no files will be written"
    [[ -z ${SCRATCH_ROOT:-} ]] || echo "[dbitm] planned scratch root: $SCRATCH_ROOT"
    echo "====== dbitm call dry-run finished ======"
    exit 0
fi

use_scratch=false
run_output=$final_dir/coverage
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
    run_output=$scratch_run/coverage
    use_scratch=true
    enable_cleanup
    mkdir -p "$scratch_input"
    echo "[dbitm] copying pooled BAMs + indexes to scratch: $scratch_input"
    for call_index in "${!call_bams[@]}"; do
        copy_bam_with_index "${call_bams[$call_index]}" "$scratch_input"
        call_bams[$call_index]=$scratch_input/$(basename "${call_bams[$call_index]}")
    done
fi
mkdir -p "$run_output"

echo "[dbitm] starting methylation calling..."
for call_index in "${!call_labels[@]}"; do
    label=${call_labels[$call_index]}
    bam_path=${call_bams[$call_index]}
    cutoff_path=${cutoff_paths[$call_index]}
    load_target_trims "$label" "$cutoff_path"
    if [[ "$assay" == smc ]]; then
        target_output=$run_output/$label
        target_log=$run_output/call.$label.log
    else
        target_output=$run_output
        target_log=$run_output/call.log
    fi
    echo "[dbitm] call target: $label"
    echo "[dbitm] trimming $label: R1=${target_r1_left},${target_r1_right} R2=${target_r2_left},${target_r2_right}"
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        python "$caller_script" \
        --assay "$assay" \
        --bam "$bam_path" \
        --reference "$reference" \
        --out-dir "$target_output" \
        --barcode-whitelist "$barcode_whitelist" \
        --chromosomes "$CALL_CHROMOSOMES" \
        --context-mode "$caller_context_mode" \
        --min-base-quality "$CALL_MIN_BASE_QUALITY" \
        --min-mapping-quality "$CALL_MIN_MAPPING_QUALITY" \
        --max-depth "$CALL_MAX_DEPTH" \
        --batch-size "$CALL_BATCH_SIZE" \
        --r1-left-trim "$target_r1_left" \
        --r1-right-trim "$target_r1_right" \
        --r2-left-trim "$target_r2_left" \
        --r2-right-trim "$target_r2_right" \
        --jobs "$CALL_JOBS" \
        > "$target_log" 2>&1
done

if [[ "$assay" == smc ]]; then
    manifest_tmp=$run_output/.spot_manifest.tsv.tmp.$$
    if ! awk -F '\t' '
        BEGIN { OFS="\t" }
        FNR == 1 {
            if (!header_written) {
                print
                header_written=1
            }
            next
        }
        $1 != "" && !seen[$1]++ { print }
    ' \
        "$run_output/watson/spot_manifest.tsv" \
        "$run_output/crick/spot_manifest.tsv" \
        > "$manifest_tmp"; then
        rm -f -- "$manifest_tmp"
        echo "[dbitm] call: failed to merge Watson/Crick spot manifests" >&2
        exit 1
    fi
    mv -f -- "$manifest_tmp" "$run_output/spot_manifest.tsv"
    mkdir -p "$run_output/host"
    mkdir -p "$run_output/watson/host" "$run_output/crick/host"
    for context in "${call_contexts[@]}"; do
        watson_cov=$run_output/watson/host/host.$context.cov
        crick_cov=$run_output/crick/host/host.$context.cov
        [[ -f "$watson_cov" ]] || : > "$watson_cov"
        [[ -f "$crick_cov" ]] || : > "$crick_cov"
        pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
            python "$merge_script" \
            --input "$watson_cov" \
            --input "$crick_cov" \
            --output "$run_output/host/host.$context.cov" \
            --chromosomes "$CALL_CHROMOSOMES" \
            --spatial \
            --reference "$reference" \
            --batch-size "$CALL_BATCH_SIZE" \
            --jobs "$CALL_JOBS"
    done
    {
        echo "===== Watson call ====="
        cat "$run_output/call.watson.log"
        echo "===== Crick call ====="
        cat "$run_output/call.crick.log"
    } > "$run_output/call.log"
fi

echo "[dbitm] call finished successfully"
if [[ "$use_scratch" == true ]]; then
    mkdir -p "$final_dir/coverage"
    if [[ "$assay" == smc ]]; then
        rm -rf -- "$final_dir/coverage/host" "$final_dir/coverage/watson" "$final_dir/coverage/crick"
    else
        rm -rf -- "$final_dir/coverage/host"
    fi
    echo "[dbitm] copying coverage result from scratch to $final_dir/coverage"
    cp -a "$run_output/." "$final_dir/coverage/"
fi
echo "[dbitm] call log: $final_dir/coverage/call.log"
if [[ "$assay" == smc ]]; then
    echo "[dbitm] Watson call result: $final_dir/coverage/watson"
    echo "[dbitm] Crick call result: $final_dir/coverage/crick"
    echo "[dbitm] merged call result: $final_dir/coverage/host"
    echo "[dbitm] merged spot manifest: $final_dir/coverage/spot_manifest.tsv"
else
    echo "[dbitm] call result: $final_dir/coverage"
fi
echo "====== dbitm call finished ======"
