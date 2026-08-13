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
mbias_mode=${MBIAS_MODE:-all}
case "$mbias_mode" in
    all|host|spike) ;;
    *) echo "[dbitm] call: MBIAS_MODE must be all, host, or spike" >&2; exit 1 ;;
esac

declare -a cutoff_targets=()
if [[ "$mbias_mode" == all || "$mbias_mode" == host ]]; then
    cutoff_targets+=(host)
fi
if [[ "$mbias_mode" == all || "$mbias_mode" == spike ]]; then
    spike_declaration=$(declare -p CALL_SPIKE_IN_REFERENCES 2>/dev/null || true)
    if [[ "$spike_declaration" != "declare -A "* ]]; then
        echo "[dbitm] call: CALL_SPIKE_IN_REFERENCES must be defined with declare -A" >&2
        exit 1
    fi
    declare -n call_spike_references=CALL_SPIKE_IN_REFERENCES
    declare -a call_spike_names=()
    if (( ${#call_spike_references[@]} > 0 )); then
        mapfile -t call_spike_names < <(
            printf '%s\n' "${!call_spike_references[@]}" | LC_ALL=C sort
        )
    fi
    for spike_name in "${call_spike_names[@]}"; do
        spike_reference=${call_spike_references[$spike_name]}
        spike_reference=${spike_reference#"${spike_reference%%[![:space:]]*}"}
        spike_reference=${spike_reference%"${spike_reference##*[![:space:]]}"}
        if [[ -n "$spike_reference" ]]; then
            cutoff_targets+=("$spike_name")
        fi
    done
fi
if (( ${#cutoff_targets[@]} == 0 )); then
    echo "[dbitm] call: no host or spike-in targets were selected for M-bias cutoffs" >&2
    exit 1
fi

declare -a cutoff_paths=()
declare -a missing_cutoff_targets=()
for cutoff_target in "${cutoff_targets[@]}"; do
    cutoff_path=$cutoff_dir/$cutoff_target.mbias.cutoffs.tsv
    if [[ ! -f "$cutoff_path" ]]; then
        if [[ "$dry_run" == true ]]; then
            echo "[dbitm] call: dry-run expects future M-bias cutoff ($cutoff_target): $cutoff_path"
            missing_cutoff_targets+=("$cutoff_target")
            continue
        else
            echo "[dbitm] call: M-bias cutoff file not found ($cutoff_target): $cutoff_path" >&2
            echo "[dbitm] call: run the mbias step first." >&2
            exit 1
        fi
    fi
    cutoff_paths+=("$cutoff_path")
done
if (( ${#missing_cutoff_targets[@]} > 0 )); then
    call_r1_left_trim=auto
    call_r1_right_trim=auto
    call_r2_left_trim=auto
    call_r2_right_trim=auto
else
    if ! trim_values=$(awk -F '\t' '
    FNR == 1 { next }
    $1 == "R1" {
        if ($7 !~ /^[0-9]+$/ || $8 !~ /^[0-9]+$/) invalid = 1
        if (!seen_r1 || $7 + 0 > r1_left) {
            r1_left = $7 + 0
            r1_left_source = FILENAME
        }
        if (!seen_r1 || $8 + 0 > r1_right) {
            r1_right = $8 + 0
            r1_right_source = FILENAME
        }
        seen_r1 = 1
    }
    $1 == "R2" {
        if ($7 !~ /^[0-9]+$/ || $8 !~ /^[0-9]+$/) invalid = 1
        if (!seen_r2 || $7 + 0 > r2_left) {
            r2_left = $7 + 0
            r2_left_source = FILENAME
        }
        if (!seen_r2 || $8 + 0 > r2_right) {
            r2_right = $8 + 0
            r2_right_source = FILENAME
        }
        seen_r2 = 1
    }
    END {
        if (invalid || !seen_r1 || !seen_r2) exit 2
        printf "%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
            r1_left, r1_right, r2_left, r2_right,
            r1_left_source, r1_right_source,
            r2_left_source, r2_right_source
    }
' "${cutoff_paths[@]}"); then
        echo "[dbitm] call: invalid M-bias cutoff file under: $cutoff_dir" >&2
        exit 1
    fi
    IFS=$'\t' read -r \
        call_r1_left_trim call_r1_right_trim \
        call_r2_left_trim call_r2_right_trim \
        r1_left_source r1_right_source r2_left_source r2_right_source \
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

use_scratch=false
run_bam=$pooled_bam
run_output=$final_dir/coverage

echo "====== dbitm call ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] reference: $reference"
echo "[dbitm] barcode whitelist: $barcode_whitelist"
echo "[dbitm] pooled BAM: $pooled_bam"
echo "[dbitm] chromosomes: $CALL_CHROMOSOMES"
echo "[dbitm] context mode: $caller_context_mode"
echo "[dbitm] M-bias cutoff files: ${#cutoff_paths[@]}/${#cutoff_targets[@]}"
echo "[dbitm] combined trimming: R1=${call_r1_left_trim},${call_r1_right_trim} R2=${call_r2_left_trim},${call_r2_right_trim}"
if (( ${#missing_cutoff_targets[@]} == 0 )); then
    echo "[dbitm] cutoff sources: R1-left=$(basename "$r1_left_source") R1-right=$(basename "$r1_right_source") R2-left=$(basename "$r2_left_source") R2-right=$(basename "$r2_right_source")"
fi
echo "[dbitm] config: $config_file"
echo "[dbitm] output directory: $final_dir/coverage"

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
