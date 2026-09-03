#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 || $# > 3 )); then
    echo "Usage: 05.mbias.sh <assay> <raw_fastq_folder> [--dry-run]" >&2
    exit 1
fi
assay=$1
raw_path=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] mbias: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
case "$assay" in
    taps|taps-v2|emseq|cabernet|smc) ;;
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
MBIAS_CUTOFF_RATE_TOLERANCE=${MBIAS_CUTOFF_RATE_TOLERANCE:-0.05}

case "$MBIAS_MODE" in
    all|host|spike) ;;
    *) echo "[dbitm] mbias: MBIAS_MODE must be all, host, or spike" >&2; exit 1 ;;
esac
for integer_name in MBIAS_HOST_MAX_RECORDS MBIAS_MAX_CYCLE MBIAS_MIN_CYCLE_COVERAGE \
    MBIAS_MIN_BASE_QUALITY MBIAS_MIN_MAPPING_QUALITY MBIAS_SAMPLING_SEED MBIAS_R1_ORIGINAL_LENGTH; do
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
    if [[ "$dry_run" == true ]]; then
        echo "[dbitm] mbias: dry-run expects future pooled output directory: $pooled_dir"
    else
        echo "[dbitm] mbias: pooled output directory not found: $pooled_dir" >&2
        echo "[dbitm] mbias: run the pool step first." >&2
        exit 1
    fi
fi

echo "====== dbitm mbias ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] mode: $MBIAS_MODE"
echo "[dbitm] pooled directory: $pooled_dir"
echo "[dbitm] output directory: $output_dir"
echo "[dbitm] host M-bias chromosomes: ${CALL_CHROMOSOMES:-all}"
echo "[dbitm] config: $config_file"

log_path=/dev/null
if [[ "$dry_run" == false ]]; then
    mkdir -p "$run_output"
    log_path=$run_output/mbias.log
    : > "$log_path"
fi

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
        --rate-tolerance "$MBIAS_CUTOFF_RATE_TOLERANCE" \
        >> "$log_path" 2>&1; then
        echo "[dbitm] mbias: cutoff inference failed; see log: $output_dir/mbias.log" >&2
        return 1
    fi
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

run_mbias_target() {
    local label=$1
    local bam_path=$2
    local reference_path=$3
    local subsample_fraction=$4
    local max_records=$5
    local chromosomes=${6:-}
    local chrom_args=()
    if [[ -n "$chromosomes" ]]; then
        chrom_args=(--chromosomes "$chromosomes")
    fi
    local target_tsv=$run_output/$label.mbias.tsv
    local target_png=$run_output/$label.mbias.png

    if [[ "$dry_run" == true ]]; then
        validate_reference "$label" "$reference_path"
        echo "[dbitm] dry-run target: $label"
        echo "[dbitm] dry-run BAM: $bam_path"
        echo "[dbitm] dry-run reference: $reference_path"
        echo "[dbitm] dry-run sampling: fraction=$subsample_fraction max-records=$max_records"
        echo "[dbitm] dry-run chromosomes: ${chromosomes:-all}"
        return 0
    fi

    if [[ ! -f "$bam_path" ]]; then
        echo "[dbitm] mbias: pooled BAM not found ($label): $bam_path" >&2
        exit 1
    fi
    validate_reference "$label" "$reference_path"
    if [[ -n "$chromosomes" ]] && ! find_bam_index "$bam_path" > /dev/null; then
        echo "[dbitm] mbias: chromosome-filtered M-bias requires an indexed BAM ($label): $bam_path" >&2
        exit 1
    fi
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
        --min-cycle-coverage "$MBIAS_MIN_CYCLE_COVERAGE" \
        --r1-original-length "$MBIAS_R1_ORIGINAL_LENGTH" \
        --min-base-quality "$MBIAS_MIN_BASE_QUALITY" \
        --min-mapping-quality "$MBIAS_MIN_MAPPING_QUALITY" \
        "${chrom_args[@]}" \
        >> "$log_path" 2>&1; then
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
    if [[ -z "${CALL_CHROMOSOMES:-}" ]]; then
        echo "[dbitm] mbias: CALL_CHROMOSOMES is required for host M-bias; set it to the calling chromosome list" >&2
        exit 1
    fi
    host_reference=$(resolve_reference "$host_reference")
    run_mbias_target \
        host "$pooled_dir/pooled.cb.bam" "$host_reference" \
        "$MBIAS_HOST_SUBSAMPLE_FRACTION" "$MBIAS_HOST_MAX_RECORDS" \
        "$CALL_CHROMOSOMES"
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

echo "[dbitm] mbias targets completed: $target_count"
if [[ "$dry_run" == true ]]; then
    echo "[dbitm] dry-run: no files will be written"
    echo "====== dbitm mbias dry-run finished ======"
else
    echo "[dbitm] mbias log: $output_dir/mbias.log"
    echo "[dbitm] mbias result: $output_dir"
    echo "====== dbitm mbias finished ======"
fi
