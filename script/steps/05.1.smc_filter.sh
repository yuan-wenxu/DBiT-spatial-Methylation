#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 || $# > 3 )); then
    echo "Usage: 05.1.smc_filter.sh <assay> <raw_fastq_folder> [--dry-run]" >&2
    exit 1
fi
assay=$1
raw_path=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] smc-filter: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
if [[ "$assay" != smc ]]; then
    echo "[dbitm] smc-filter: this stage is only available for assay=smc" >&2
    exit 1
fi
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] smc-filter: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] smc-filter: config file not found: $config_file" >&2
    exit 1
fi
source "$config_file"

filter_script=$REPO_DIR/script/steps/python/05.1.smc_xxx_filter.py
if [[ ! -f "$filter_script" ]]; then
    echo "[dbitm] smc-filter: filter script not found: $filter_script" >&2
    exit 1
fi
smc_filter_threads=${SMC_FILTER_THREADS:-}
if [[ ! "$smc_filter_threads" =~ ^[1-9][0-9]*$ ]]; then
    echo "[dbitm] smc-filter: SMC_FILTER_THREADS must be greater than zero" >&2
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
            return 0
        fi
    done
    return 1
}

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
pooled_dir=$final_dir/pooled
cutoff_dir=$final_dir/mbias
output_dir=$final_dir/smc_xxx_filter
target_count=0
declare -a target_logs=()

run_target() {
    local label=$1
    local bam_path=$2
    local reference_path=$3
    local cutoff_path=$4
    local output_bam=$output_dir/$(basename "$bam_path")
    local target_log=$output_dir/$label.log
    local trims r1_left r1_right r2_left r2_right

    if [[ ! "$label" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        echo "[dbitm] smc-filter: invalid target label: $label" >&2
        return 1
    fi
    if [[ "$dry_run" == true ]]; then
        [[ -f "$bam_path" ]] || echo "[dbitm] smc-filter: dry-run expects future BAM ($label): $bam_path"
        [[ -f "$cutoff_path" ]] || echo "[dbitm] smc-filter: dry-run expects future M-bias cutoff ($label): $cutoff_path"
        echo "[dbitm] smc-filter dry-run target: $label"
        echo "[dbitm] smc-filter dry-run input: $bam_path"
        echo "[dbitm] smc-filter dry-run output: $output_bam"
        target_count=$((target_count + 1))
        return 0
    fi
    if [[ ! -f "$bam_path" ]]; then
        echo "[dbitm] smc-filter: pooled BAM not found ($label): $bam_path" >&2
        return 1
    fi
    if ! find_bam_index "$bam_path"; then
        echo "[dbitm] smc-filter: BAM index not found ($label): $bam_path" >&2
        return 1
    fi
    if [[ ! -f "$reference_path" || ! -f "$reference_path.fai" ]]; then
        echo "[dbitm] smc-filter: reference FASTA/.fai not found ($label): $reference_path" >&2
        return 1
    fi
    if [[ ! -f "$cutoff_path" ]]; then
        echo "[dbitm] smc-filter: M-bias cutoff not found ($label): $cutoff_path" >&2
        return 1
    fi
    if [[ ! -s "$cutoff_path" ]]; then
        echo "[dbitm] smc-filter: no stable M-bias cutoff for $label; using zero trims"
        r1_left=0
        r1_right=0
        r2_left=0
        r2_right=0
    else
        if ! trims=$(load_trims "$cutoff_path"); then
            echo "[dbitm] smc-filter: invalid M-bias cutoff ($label): $cutoff_path" >&2
            return 1
        fi
        IFS=$'\t' read -r r1_left r1_right r2_left r2_right <<< "$trims"
    fi

    echo "[dbitm] smc-filter target: $label"
    echo "[dbitm] smc-filter trimming: R1=$r1_left,$r1_right R2=$r2_left,$r2_right"
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        python "$filter_script" \
        --bam "$bam_path" \
        --reference "$reference_path" \
        --output-bam "$output_bam" \
        --r1-left-trim "$r1_left" \
        --r1-right-trim "$r1_right" \
        --r2-left-trim "$r2_left" \
        --r2-right-trim "$r2_right" \
        --jobs "$smc_filter_threads" \
        > "$target_log" 2>&1
    target_logs+=("$target_log")
    target_count=$((target_count + 1))
    echo "[dbitm] smc-filter output: $output_bam"
}

host_reference=${CALL_REFERENCE:-}
if [[ -z "$host_reference" ]]; then
    echo "[dbitm] smc-filter: CALL_REFERENCE is required" >&2
    exit 1
fi
host_reference=$(resolve_reference "$host_reference")

echo "====== dbitm smc-filter ======"
echo "[dbitm] pooled directory: $pooled_dir"
echo "[dbitm] M-bias directory: $cutoff_dir"
echo "[dbitm] output directory: $output_dir"
echo "[dbitm] filter processes per BAM: $smc_filter_threads"
echo "[dbitm] config: $config_file"

if [[ "$dry_run" == false ]]; then
    mkdir -p "$output_dir"
fi
for conversion_class in watson crick; do
    run_target \
        "host.$conversion_class" \
        "$pooled_dir/pooled.$conversion_class.cb.bam" \
        "$host_reference" \
        "$cutoff_dir/host.$conversion_class.mbias.cutoffs.tsv"
done

spike_call_mode=${SPIKE_CALL_MODE:-all}
case "$spike_call_mode" in
    all|mito|spike) ;;
    *) echo "[dbitm] smc-filter: SPIKE_CALL_MODE must be all, mito, or spike" >&2; exit 1 ;;
esac
if [[ "$spike_call_mode" == all || "$spike_call_mode" == spike ]]; then
    spike_declaration=$(declare -p CALL_SPIKE_IN_REFERENCES 2>/dev/null || true)
    if [[ "$spike_declaration" != "declare -A "* ]]; then
        echo "[dbitm] smc-filter: CALL_SPIKE_IN_REFERENCES must be defined with declare -A" >&2
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
        for conversion_class in watson crick; do
            run_target \
                "$spike_name.$conversion_class" \
                "$pooled_dir/pooled.$conversion_class.$spike_name.bam" \
                "$spike_reference" \
                "$cutoff_dir/$spike_name.$conversion_class.mbias.cutoffs.tsv"
        done
    done
fi

if [[ "$dry_run" == true ]]; then
    echo "[dbitm] dry-run: no files will be written"
    echo "====== dbitm smc-filter dry-run finished ======"
    exit 0
fi

combined_log=$output_dir/filter.log
combined_temp=$output_dir/.filter.log.tmp.$$
{
    for target_log in "${target_logs[@]}"; do
        echo "===== $(basename "$target_log" .log) ====="
        cat "$target_log"
    done
} > "$combined_temp"
mv -f -- "$combined_temp" "$combined_log"

echo "[dbitm] smc-filter targets: $target_count"
echo "[dbitm] smc-filter log: $combined_log"
echo "====== dbitm smc-filter finished ======"
