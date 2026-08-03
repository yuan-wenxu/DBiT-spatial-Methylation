#!/usr/bin/env bash
set -euo pipefail

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ ${use_scratch:-false} == true && -n ${scratch_run:-} && -d $scratch_run ]]; then
        rm -rf -- "$scratch_run" || echo "[dbitm] mbias: warning: failed to clean scratch directory: $scratch_run" >&2
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
    echo "Usage: 05.mbias.sh <assay> <raw_fastq_folder>" >&2
    exit 1
fi
assay=$1
raw_path=$2
case "$assay" in
    taps|taps-v2|emseq) ;;
    *) echo "[dbitm] mbias: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] mbias: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] mbias: config file not found: $config_file" >&2
    exit 1
fi
source "$config_file"
MBIAS_R1_ORIGINAL_LENGTH=${MBIAS_R1_ORIGINAL_LENGTH:-150}

case "$MBIAS_MODE" in
    all|host|spike) ;;
    *) echo "[dbitm] mbias: MBIAS_MODE must be all, host, or spike" >&2; exit 1 ;;
esac
for integer_name in MBIAS_HOST_MAX_RECORDS MBIAS_MAX_CYCLE MBIAS_MIN_BASE_QUALITY \
    MBIAS_MIN_MAPPING_QUALITY MBIAS_SAMPLING_SEED MBIAS_R1_ORIGINAL_LENGTH; do
    integer_value=${!integer_name}
    if [[ ! "$integer_value" =~ ^[0-9]+$ ]]; then
        echo "[dbitm] mbias: $integer_name must be a non-negative integer" >&2
        exit 1
    fi
done
if (( MBIAS_MAX_CYCLE == 0 )); then
    echo "[dbitm] mbias: MBIAS_MAX_CYCLE must be greater than zero" >&2
    exit 1
fi

resolve_reference() {
    local reference_path=$1
    if [[ "$reference_path" != /* ]]; then
        reference_path=$REPO_DIR/$reference_path
    fi
    realpath -m "$reference_path"
}

validate_reference() {
    local label=$1
    local reference_path=$2
    if [[ ! -f "$reference_path" ]]; then
        echo "[dbitm] mbias: reference FASTA not found ($label): $reference_path" >&2
        exit 1
    fi
}

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
pooled_dir=$final_dir/pooled
output_dir=$final_dir/mbias
run_output=$output_dir
python_script=$REPO_DIR/script/steps/python/05.mbias.py
cutoff_script=$REPO_DIR/script/steps/python/05.mbias_cutoff.py
if [[ ! -f "$python_script" ]]; then
    echo "[dbitm] mbias: implementation not found: $python_script" >&2
    exit 1
fi
if [[ ! -f "$cutoff_script" ]]; then
    echo "[dbitm] mbias: cutoff implementation not found: $cutoff_script" >&2
    exit 1
fi
if [[ ! -d "$pooled_dir" ]]; then
    echo "[dbitm] mbias: pooled output directory not found: $pooled_dir" >&2
    echo "[dbitm] mbias: run the pool step first." >&2
    exit 1
fi

use_scratch=false
if [[ -n ${SCRATCH_ROOT:-} ]]; then
    if [[ "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] mbias: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    echo "[dbitm] scratch root: $SCRATCH_ROOT"
    mkdir -p "$SCRATCH_ROOT"
    scratch_root=$(realpath "$SCRATCH_ROOT")
    run_id=${SLURM_JOB_ID:-mbias_$$}
    scratch_run=$scratch_root/dbitm/$run_id
    scratch_input=$scratch_run/pooled
    scratch_references=$scratch_run/references
    run_output=$scratch_run/mbias
    use_scratch=true
    enable_cleanup
    mkdir -p "$scratch_input" "$scratch_references" "$run_output"
    if [[ -d "$output_dir" ]]; then
        echo "[dbitm] copying existing M-bias outputs to scratch: $run_output"
        cp -a "$output_dir/." "$run_output/"
    fi
fi

mkdir -p "$run_output"
log_path=$run_output/mbias.log
: > "$log_path"

echo "====== dbitm mbias ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] mode: $MBIAS_MODE"
echo "[dbitm] pooled directory: $pooled_dir"
echo "[dbitm] output directory: $output_dir"
echo "[dbitm] config: $config_file"

infer_target_cutoffs() {
    local label=$1
    local cutoff_path=$run_output/$label.mbias.cutoffs.tsv
    if [[ -f "$cutoff_path" ]]; then
        echo "[dbitm] M-bias cutoffs already exist, skipping: $cutoff_path"
        return 0
    fi
    echo "[dbitm] inferring M-bias trimming cutoffs: $label"
    if ! pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        python "$cutoff_script" \
        --mbias-tsv "$run_output/$label.mbias.tsv" \
        --output "$cutoff_path" \
        --r1-original-length "$MBIAS_R1_ORIGINAL_LENGTH" \
        >> "$log_path" 2>&1; then
        if [[ "$use_scratch" == true ]]; then
            mkdir -p "$output_dir"
            cp -a "$log_path" "$output_dir/mbias.log" || true
        fi
        echo "[dbitm] mbias: cutoff inference failed; see log: $output_dir/mbias.log" >&2
        return 1
    fi
}

run_mbias_target() {
    local label=$1
    local bam_path=$2
    local reference_path=$3
    local subsample_fraction=$4
    local max_records=$5
    local scratch_bam scratch_reference
    local target_tsv=$run_output/$label.mbias.tsv
    local target_png=$run_output/$label.mbias.png

    if [[ ! -f "$bam_path" ]]; then
        echo "[dbitm] mbias: pooled BAM not found ($label): $bam_path" >&2
        exit 1
    fi
    validate_reference "$label" "$reference_path"
    if [[ -f "$target_tsv" && -f "$target_png" ]]; then
        echo "[dbitm] mbias target already exists, skipping: $label"
        echo "[mbias] skip_existing_outputs=$label" >> "$log_path"
        infer_target_cutoffs "$label"
        return 0
    fi
    if [[ -f "$target_tsv" || -f "$target_png" ]]; then
        echo "[dbitm] mbias: partial output exists for $label under: $run_output" >&2
        return 1
    fi
    if [[ "$use_scratch" == true ]]; then
        scratch_bam=$scratch_input/$(basename "$bam_path")
        scratch_reference=$scratch_references/$label.fa
        echo "[dbitm] copying pooled BAM to scratch ($label): $scratch_bam"
        cp -L -- "$bam_path" "$scratch_bam"
        echo "[dbitm] copying reference FASTA + index to scratch ($label): $scratch_reference"
        cp -L -- "$reference_path" "$scratch_reference"
        cp -L -- "$reference_path.fai" "$scratch_reference.fai"
        bam_path=$scratch_bam
        reference_path=$scratch_reference
    fi
    echo "[dbitm] mbias target: $label"
    if ! pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        python "$python_script" \
        --assay "$assay" \
        --bam "$bam_path" \
        --reference "$reference_path" \
        --label "$label" \
        --output-dir "$run_output" \
        --subsample-fraction "$subsample_fraction" \
        --max-records "$max_records" \
        --seed "$MBIAS_SAMPLING_SEED" \
        --max-cycle "$MBIAS_MAX_CYCLE" \
        --r1-original-length "$MBIAS_R1_ORIGINAL_LENGTH" \
        --min-base-quality "$MBIAS_MIN_BASE_QUALITY" \
        --min-mapping-quality "$MBIAS_MIN_MAPPING_QUALITY" \
        >> "$log_path" 2>&1; then
        if [[ "$use_scratch" == true ]]; then
            mkdir -p "$output_dir"
            cp -a "$log_path" "$output_dir/mbias.log" || true
        fi
        echo "[dbitm] mbias: target failed ($label); see log: $output_dir/mbias.log" >&2
        return 1
    fi
    infer_target_cutoffs "$label"
}

target_count=0
if [[ "$MBIAS_MODE" == all || "$MBIAS_MODE" == host ]]; then
    host_reference=${CALL_REFERENCE:-}
    if [[ -z "$host_reference" ]]; then
        echo "[dbitm] mbias: CALL_REFERENCE is required for host M-bias" >&2
        exit 1
    fi
    host_reference=$(resolve_reference "$host_reference")
    run_mbias_target \
        host "$pooled_dir/pooled.cb.bam" "$host_reference" \
        "$MBIAS_HOST_SUBSAMPLE_FRACTION" "$MBIAS_HOST_MAX_RECORDS"
    target_count=$((target_count + 1))
fi

if [[ "$MBIAS_MODE" == all || "$MBIAS_MODE" == spike ]]; then
    spike_declaration=$(declare -p CALL_SPIKE_IN_REFERENCES 2>/dev/null || true)
    if [[ "$spike_declaration" != "declare -A "* ]]; then
        echo "[dbitm] mbias: CALL_SPIKE_IN_REFERENCES must be defined with declare -A" >&2
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
        if [[ -z "$spike_reference" ]]; then
            continue
        fi
        spike_reference=$(resolve_reference "$spike_reference")
        run_mbias_target \
            "$spike_name" "$pooled_dir/pooled.$spike_name.bam" "$spike_reference" 1 0
        target_count=$((target_count + 1))
    done
fi

if (( target_count == 0 )); then
    echo "[dbitm] mbias: no host or spike-in targets were selected" >&2
    exit 1
fi

if [[ "$use_scratch" == true ]]; then
    mkdir -p "$output_dir"
    echo "[dbitm] copying M-bias result from scratch to $output_dir"
    cp -a "$run_output/." "$output_dir/"
fi

echo "[dbitm] mbias targets completed: $target_count"
echo "[dbitm] mbias log: $output_dir/mbias.log"
echo "[dbitm] mbias result: $output_dir"
echo "====== dbitm mbias finished ======"
