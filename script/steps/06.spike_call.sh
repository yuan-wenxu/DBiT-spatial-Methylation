#!/usr/bin/env bash
set -euo pipefail

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
    taps|taps-v2|emseq|cabernet|smc) ;;
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
smc_filter_dir=$final_dir/smc_xxx_filter
cutoff_dir=$final_dir/mbias
output_dir=$final_dir/coverage
caller_script=$REPO_DIR/script/steps/python/06.spike_caller.py
merge_script=$REPO_DIR/script/steps/python/06.merge_cov.py
if [[ ! -f "$caller_script" ]]; then
    echo "[dbitm] spike-call: caller script not found: $caller_script" >&2
    exit 1
fi
if [[ "$assay" == smc && ! -f "$merge_script" ]]; then
    echo "[dbitm] spike-call: coverage merge script not found: $merge_script" >&2
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

run_output=$output_dir

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
    local trims r1_left r1_right r2_left r2_right
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
        echo "[dbitm] spike-call: input BAM not found ($output_name): $bam_path" >&2
        if [[ "$assay" == smc ]]; then
            echo "[dbitm] spike-call: run the smc-filter step first." >&2
        fi
        return 1
    fi
    if ! find_bam_index "$bam_path" > /dev/null; then
        echo "[dbitm] spike-call: BAM index not found ($output_name): $bam_path" >&2
        return 1
    fi
    if [[ ! -f "$cutoff_path" ]]; then
        echo "[dbitm] spike-call: M-bias cutoff not found ($output_name): $cutoff_path" >&2
        return 1
    fi
    if [[ ! -s "$cutoff_path" ]]; then
        echo "[dbitm] spike-call: no suitable M-bias cutoff for $output_name; continuing without trimming"
        r1_left=0
        r1_right=0
        r2_left=0
        r2_right=0
    else
        if ! trims=$(load_trims "$cutoff_path"); then
            echo "[dbitm] spike-call: invalid M-bias cutoff ($output_name): $cutoff_path" >&2
            return 1
        fi
        IFS=$'\t' read -r r1_left r1_right r2_left r2_right <<< "$trims"
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

run_conversion_split_target() {
    local output_name=$1
    local pooled_suffix=$2
    local reference_path=$3
    local chromosomes=$4
    local cutoff_name=$5
    local conversion_class context branch_cov conversion_bam
    local -a contexts=()

    if [[ "$assay" != smc ]]; then
        run_target \
            "$output_name" "$pooled_dir/pooled.$pooled_suffix.bam" \
            "$reference_path" "$chromosomes" "$cutoff_name"
        return
    fi

    for conversion_class in watson crick; do
        conversion_bam=$pooled_dir/pooled.$conversion_class.$pooled_suffix.bam
        if [[ "$assay" == smc ]]; then
            conversion_bam=$smc_filter_dir/$(basename "$conversion_bam")
        fi
        run_target \
            "$output_name.$conversion_class" \
            "$conversion_bam" \
            "$reference_path" "$chromosomes" \
            "$cutoff_name.$conversion_class"
    done
    if [[ "$dry_run" == true ]]; then
        echo "[dbitm] spike-call planned merge: $output_name Watson + Crick"
        return
    fi

    case "$caller_context_mode" in
        cg) contexts=(CG) ;;
        ch) contexts=(CA CC CT) ;;
        both) contexts=(CG CA CC CT) ;;
    esac
    for context in "${contexts[@]}"; do
        for conversion_class in watson crick; do
            branch_cov=$run_output/$output_name.$conversion_class.$context.cov
            [[ -f "$branch_cov" ]] || : > "$branch_cov"
        done
        pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
            python "$merge_script" \
            --input "$run_output/$output_name.watson.$context.cov" \
            --input "$run_output/$output_name.crick.$context.cov" \
            --output "$run_output/$output_name.$context.cov" \
            --chromosomes "$chromosomes" \
            >> "$log_path" 2>&1
    done
}

echo "====== dbitm spike-call ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] mode: $spike_call_mode"
echo "[dbitm] context mode: $caller_context_mode"
echo "[dbitm] pooled directory: $pooled_dir"
if [[ "$assay" == smc ]]; then
    echo "[dbitm] SmC filtered BAM directory: $smc_filter_dir"
fi
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
        run_conversion_split_target \
            host_mito cb "$host_reference" "$mito_chromosomes" host
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
        run_conversion_split_target \
            "$spike_name" "$spike_name" "$spike_reference" \
            "$spike_chromosomes" "$spike_name"
    done
fi

if (( target_count == 0 )); then
    echo "[dbitm] spike-call: no mito or spike-in targets were called; nothing to do"
    echo "====== dbitm spike-call finished ======"
    exit 0
fi

if [[ "$dry_run" == true ]]; then
    echo "[dbitm] dry-run: no files will be written"
    echo "====== dbitm spike-call dry-run finished ======"
    exit 0
fi

echo "[dbitm] spike-call targets: $target_count"
echo "[dbitm] spike-call log: $output_dir/spike_call.log"
echo "====== dbitm spike-call finished ======"
