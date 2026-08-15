#!/usr/bin/env bash
set -euo pipefail

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ ${use_scratch:-false} == true && -n ${scratch_run:-} && -d $scratch_run ]]; then
        if (( status != 0 )) && [[ -n ${run_output:-} && -d $run_output && -n ${output_dir:-} ]]; then
            echo "[dbitm] spike-call: recovering scratch results after exit status $status: $output_dir" >&2
            if mkdir -p "$output_dir" && cp -a "$run_output/." "$output_dir/"; then
                echo "[dbitm] spike-call: scratch results recovered: $output_dir" >&2
            else
                echo "[dbitm] spike-call: warning: failed to recover scratch results: $run_output" >&2
            fi
        fi
        rm -rf -- "$scratch_run" || echo "[dbitm] spike-call: warning: failed to clean scratch directory: $scratch_run" >&2
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
    echo "Usage: 06.spike_call.sh <assay> <raw_fastq_folder> [--dry-run]" >&2
    exit 1
fi
assay=$1
raw_path=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] spike-call: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
case "$assay" in
    taps|taps-v2|emseq|cabernet) ;;
    *) echo "[dbitm] spike-call: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] spike-call: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] spike-call: config file not found: $config_file" >&2
    exit 1
fi
source "$config_file"

spike_call_mode=${SPIKE_CALL_MODE:-all}
case "$spike_call_mode" in
    all|mito|spike) ;;
    *) echo "[dbitm] spike-call: SPIKE_CALL_MODE must be all, mito, or spike" >&2; exit 1 ;;
esac
caller_context_mode=${CALL_CONTEXT_MODE:-both}
case "$caller_context_mode" in
    cg|ch|both) ;;
    *) echo "[dbitm] spike-call: CALL_CONTEXT_MODE must be cg, ch, or both" >&2; exit 1 ;;
esac

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
pooled_dir=$final_dir/pooled
cutoff_dir=$final_dir/mbias
output_dir=$final_dir/coverage
caller_script=$REPO_DIR/script/steps/python/06.spike_caller.py
if [[ ! -f "$caller_script" ]]; then
    echo "[dbitm] spike-call: caller script not found: $caller_script" >&2
    exit 1
fi

resolve_reference() {
    local reference_path=$1
    if [[ "$reference_path" != /* ]]; then
        reference_path=$REPO_DIR/$reference_path
    fi
    realpath -m "$reference_path"
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

find_bam_index() {
    local bam_path=$1
    local candidate
    for candidate in "$bam_path.bai" "${bam_path%.bam}.bai" "$bam_path.csi"; do
        if [[ -f "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

use_scratch=false
run_output=$output_dir
if [[ "$dry_run" == false && -n ${SCRATCH_ROOT:-} ]]; then
    if [[ "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] spike-call: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    scratch_root=$(realpath -m "$SCRATCH_ROOT")
    run_id=${SLURM_JOB_ID:-spike_call_$$}
    scratch_run=$scratch_root/dbitm/$run_id
    scratch_input=$scratch_run/input
    run_output=$scratch_run/coverage
    use_scratch=true
    enable_cleanup
    mkdir -p "$scratch_input" "$run_output"
fi

target_count=0
log_path=$run_output/spike_call.log
if [[ "$dry_run" == false ]]; then
    mkdir -p "$run_output"
    : > "$log_path"
fi

run_target() {
    local output_name=$1
    local bam_path=$2
    local reference_path=$3
    local chromosomes=$4
    local cutoff_name=$5
    local cutoff_path=$cutoff_dir/$cutoff_name.mbias.cutoffs.tsv
    local trims r1_left r1_right r2_left r2_right bam_index
    local run_bam=$bam_path
    local run_reference=$reference_path

    if [[ ! "$output_name" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        echo "[dbitm] spike-call: invalid output name: $output_name" >&2
        return 1
    fi
    if [[ ! -f "$reference_path" || ! -f "$reference_path.fai" ]]; then
        echo "[dbitm] spike-call: reference FASTA/.fai not found ($output_name): $reference_path" >&2
        return 1
    fi
    if [[ "$dry_run" == true ]]; then
        [[ -f "$bam_path" ]] || echo "[dbitm] spike-call: dry-run expects future BAM ($output_name): $bam_path"
        [[ -f "$cutoff_path" ]] || echo "[dbitm] spike-call: dry-run expects future M-bias cutoff ($output_name): $cutoff_path"
        echo "[dbitm] spike-call dry-run target: $output_name"
        echo "[dbitm] spike-call dry-run chromosomes: $chromosomes"
        echo "[dbitm] spike-call dry-run output: $output_dir/$output_name.{CG,CH}.cov"
        target_count=$((target_count + 1))
        return 0
    fi
    if [[ ! -f "$bam_path" ]]; then
        echo "[dbitm] spike-call: pooled BAM not found ($output_name): $bam_path" >&2
        return 1
    fi
    if ! bam_index=$(find_bam_index "$bam_path"); then
        echo "[dbitm] spike-call: BAM index not found ($output_name): $bam_path" >&2
        return 1
    fi
    if [[ ! -f "$cutoff_path" ]]; then
        echo "[dbitm] spike-call: M-bias cutoff not found ($output_name): $cutoff_path" >&2
        return 1
    fi
    if [[ ! -s "$cutoff_path" ]]; then
        echo "[dbitm] spike-call: no M-bias data for $output_name; skipping target"
        return 0
    fi
    if ! trims=$(load_trims "$cutoff_path"); then
        echo "[dbitm] spike-call: invalid M-bias cutoff ($output_name): $cutoff_path" >&2
        return 1
    fi
    IFS=$'\t' read -r r1_left r1_right r2_left r2_right <<< "$trims"

    if [[ "$use_scratch" == true ]]; then
        local target_input=$scratch_input/$output_name
        mkdir -p "$target_input"
        run_bam=$target_input/input.bam
        run_reference=$target_input/reference.fa
        cp -L -- "$bam_path" "$run_bam"
        case "$bam_index" in
            *.csi) cp -L -- "$bam_index" "$run_bam.csi" ;;
            *) cp -L -- "$bam_index" "$run_bam.bai" ;;
        esac
        cp -L -- "$reference_path" "$run_reference"
        cp -L -- "$reference_path.fai" "$run_reference.fai"
    fi

    echo "[dbitm] spike-call target: $output_name"
    echo "[dbitm] spike-call trimming: R1=$r1_left,$r1_right R2=$r2_left,$r2_right"
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        python "$caller_script" \
        --assay "$assay" \
        --bam "$run_bam" \
        --reference "$run_reference" \
        --output-dir "$run_output" \
        --output-name "$output_name" \
        --chromosomes "$chromosomes" \
        --context-mode "$caller_context_mode" \
        --min-base-quality "$CALL_MIN_BASE_QUALITY" \
        --min-mapping-quality "$CALL_MIN_MAPPING_QUALITY" \
        --max-depth "$CALL_MAX_DEPTH" \
        --r1-left-trim "$r1_left" \
        --r1-right-trim "$r1_right" \
        --r2-left-trim "$r2_left" \
        --r2-right-trim "$r2_right" \
        >> "$log_path" 2>&1
    target_count=$((target_count + 1))
}

echo "====== dbitm spike-call ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] mode: $spike_call_mode"
echo "[dbitm] context mode: $caller_context_mode"
echo "[dbitm] pooled directory: $pooled_dir"
echo "[dbitm] output directory: $output_dir"
echo "[dbitm] config: $config_file"

if [[ "$spike_call_mode" == all || "$spike_call_mode" == mito ]]; then
    mito_chromosomes=${CALL_MITO_CHROMOSOMES:-}
    if [[ -n "$mito_chromosomes" ]]; then
        host_reference=${CALL_REFERENCE:-}
        if [[ -z "$host_reference" ]]; then
            echo "[dbitm] spike-call: CALL_REFERENCE is required for mito calling" >&2
            exit 1
        fi
        host_reference=$(resolve_reference "$host_reference")
        run_target host_mito "$pooled_dir/pooled.cb.bam" "$host_reference" "$mito_chromosomes" host
    fi
fi

if [[ "$spike_call_mode" == all || "$spike_call_mode" == spike ]]; then
    spike_declaration=$(declare -p CALL_SPIKE_IN_REFERENCES 2>/dev/null || true)
    if [[ "$spike_declaration" != "declare -A "* ]]; then
        echo "[dbitm] spike-call: CALL_SPIKE_IN_REFERENCES must be defined with declare -A" >&2
        exit 1
    fi
    declare -n spike_references=CALL_SPIKE_IN_REFERENCES
    declare -a spike_names=()
    if (( ${#spike_references[@]} > 0 )); then
        mapfile -t spike_names < <(
            printf '%s\n' "${!spike_references[@]}" | LC_ALL=C sort
        )
    fi
    for spike_name in "${spike_names[@]}"; do
        spike_reference=${spike_references[$spike_name]}
        spike_reference=${spike_reference#"${spike_reference%%[![:space:]]*}"}
        spike_reference=${spike_reference%"${spike_reference##*[![:space:]]}"}
        [[ -n "$spike_reference" ]] || continue
        spike_reference=$(resolve_reference "$spike_reference")
        if [[ ! -f "$spike_reference.fai" ]]; then
            echo "[dbitm] spike-call: reference .fai not found ($spike_name): $spike_reference.fai" >&2
            exit 1
        fi
        spike_chromosomes=$(cut -f 1 "$spike_reference.fai" | paste -sd, -)
        if [[ -z "$spike_chromosomes" ]]; then
            echo "[dbitm] spike-call: reference contains no contigs ($spike_name): $spike_reference" >&2
            exit 1
        fi
        run_target "$spike_name" "$pooled_dir/pooled.$spike_name.bam" "$spike_reference" "$spike_chromosomes" "$spike_name"
    done
fi

if (( target_count == 0 )); then
    echo "[dbitm] spike-call: no mito or spike-in targets were called; nothing to do"
    echo "====== dbitm spike-call finished ======"
    exit 0
fi

if [[ "$dry_run" == true ]]; then
    [[ -z ${SCRATCH_ROOT:-} ]] || echo "[dbitm] planned scratch root: $SCRATCH_ROOT"
    echo "[dbitm] dry-run: no files will be written"
    echo "====== dbitm spike-call dry-run finished ======"
    exit 0
fi

if [[ "$use_scratch" == true ]]; then
    mkdir -p "$output_dir"
    cp -a "$run_output/." "$output_dir/"
fi
echo "[dbitm] spike-call targets: $target_count"
echo "[dbitm] spike-call log: $output_dir/spike_call.log"
echo "====== dbitm spike-call finished ======"
