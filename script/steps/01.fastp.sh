#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 || $# > 3 )); then
    echo "Usage: 01.fastp.sh <assay> <raw_fastq_folder> [--dry-run]" >&2
    exit 1
fi
assay=$1
raw_path=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] fastp: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
case "$assay" in
    taps|taps-v2|emseq|cabernet|smc) ;;
    *) echo "[dbitm] fastp: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] fastp: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] fastp: config file not found: $config_file" >&2
    exit 1
fi
source "$config_file"

raw_abs=$(realpath "$raw_path")

if [[ ! "$FASTP_THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "[dbitm] fastp: FASTP_THREADS must be >= 1" >&2
    exit 1
fi

final_dir=$(dirname "$raw_abs")/dbitm
run_input=$raw_abs
run_output=$final_dir/fastp
merge_dir=$final_dir/merge
merged_r1=$merge_dir/R1.merged.fastq.gz
merged_r2=$merge_dir/R2.merged.fastq.gz

declare -a r1_files=()
declare -a r2_files=()

collect_fastq_pairs() {
    local input_dir=$1
    local path filename pair_key mate i
    declare -a fastq_candidates=()
    declare -a r1_keys=()
    declare -a r2_keys=()

    shopt -s nullglob nocaseglob
    fastq_candidates=(
        "$input_dir"/*.fastq.gz
        "$input_dir"/*.fq.gz
        "$input_dir"/*.fastq
        "$input_dir"/*.fq
    )
    shopt -u nullglob nocaseglob

    if (( ${#fastq_candidates[@]} == 0 )); then
        echo "[dbitm] fastp: no FASTQ files found: $input_dir" >&2
        exit 1
    fi

    for path in "${fastq_candidates[@]}"; do
        [[ -f "$path" ]] || continue
        filename=$(basename "$path")
        if [[ ! "$filename" =~ ^((.*[^[:alnum:]])?)[Rr]([12])([^[:alnum:]].*|$)$ ]]; then
            echo "[dbitm] fastp: FASTQ filename must contain an independent R1 or R2 token: $filename" >&2
            exit 1
        fi
        pair_key=${BASH_REMATCH[1]}R#${BASH_REMATCH[4]}
        mate=${BASH_REMATCH[3]}
        if [[ "$mate" == 1 ]]; then
            r1_files+=("$path")
            r1_keys+=("$pair_key")
        else
            r2_files+=("$path")
            r2_keys+=("$pair_key")
        fi
    done

    if (( ${#r1_files[@]} == 0 || ${#r2_files[@]} == 0 )); then
        echo "[dbitm] fastp: input directory must contain at least one R1/R2 FASTQ pair" >&2
        exit 1
    fi
    if (( ${#r1_files[@]} != ${#r2_files[@]} )); then
        echo "[dbitm] fastp: found ${#r1_files[@]} R1 files but ${#r2_files[@]} R2 files" >&2
        exit 1
    fi
    for (( i = 0; i < ${#r1_files[@]}; i++ )); do
        if [[ "${r1_keys[$i]}" != "${r2_keys[$i]}" ]]; then
            echo "[dbitm] fastp: FASTQ pair mismatch: $(basename "${r1_files[$i]}") and $(basename "${r2_files[$i]}")" >&2
            exit 1
        fi
    done
}

collect_fastq_pairs "$run_input"

echo "====== dbitm fastp ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] input directory: $raw_abs"
echo "[dbitm] config: ${config_file}"
echo "[dbitm] threads: $FASTP_THREADS"
echo "[dbitm] output directory: $final_dir/fastp"

if [[ "$dry_run" == true ]]; then
    echo "[dbitm] dry-run: no files will be written"
    echo "[dbitm] input pairs: ${#r1_files[@]}"
    for (( i = 0; i < ${#r1_files[@]}; i++ )); do
        echo "[dbitm] input pair $((i + 1)) R1: ${r1_files[$i]}"
        echo "[dbitm] input pair $((i + 1)) R2: ${r2_files[$i]}"
    done
    if (( ${#r1_files[@]} > 1 )); then
        if [[ -f "$merged_r1" && -f "$merged_r2" ]]; then
            echo "[dbitm] planned merge: reuse existing R1/R2 merged FASTQs"
        else
            echo "[dbitm] planned merge: ${#r1_files[@]} pairs -> one R1/R2 pair"
        fi
    fi
    echo "[dbitm] planned command: fastp -> $final_dir/fastp"
    echo "====== dbitm fastp dry-run finished ======"
    exit 0
fi

mkdir -p "$final_dir"
rm -rf -- "$run_output"
mkdir -p "$run_output"
if (( ${#r1_files[@]} > 1 )); then
    mkdir -p "$merge_dir"
fi

merge_fastqs() {
    local output=$1
    shift
    local input
    local all_gzip=true

    for input in "$@"; do
        if [[ "${input,,}" != *.gz ]]; then
            all_gzip=false
            break
        fi
    done
    if [[ "$all_gzip" == true ]]; then
        if ! cat -- "$@" > "$output"; then
            rm -f -- "$output"
            return 1
        fi
    else
        if ! {
                for input in "$@"; do
                    if [[ "${input,,}" == *.gz ]]; then
                        pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default pigz -cd -- "$input"
                    else
                        cat -- "$input"
                    fi
                done
            } | pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default pigz -p "$FASTP_THREADS" > "$output"
        then
            rm -f -- "$output"
            return 1
        fi
    fi
}

r1=${r1_files[0]}
r2=${r2_files[0]}
echo "[dbitm] input pairs: ${#r1_files[@]}"
if (( ${#r1_files[@]} > 1 )); then
    if [[ -f "$merged_r1" && -f "$merged_r2" ]]; then
        echo "[dbitm] reusing merged R1: $merged_r1"
        echo "[dbitm] reusing merged R2: $merged_r2"
    else
        if [[ -f "$merged_r1" || -f "$merged_r2" ]]; then
            echo "[dbitm] incomplete merged FASTQ pair found; regenerating both files"
        fi
        echo "[dbitm] merging ${#r1_files[@]} R1 files: $merged_r1"
        merge_fastqs "$merged_r1" "${r1_files[@]}"
        echo "[dbitm] merging ${#r2_files[@]} R2 files: $merged_r2"
        merge_fastqs "$merged_r2" "${r2_files[@]}"
    fi
    r1=$merged_r1
    r2=$merged_r2
fi

echo "[dbitm] input R1: $r1"
echo "[dbitm] input R2: $r2"
echo "[dbitm] starting fastp..."

pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default fastp \
    -i "$r1" \
    -I "$r2" \
    -o "$run_output/R1.filtered.fastq.gz" \
    -O "$run_output/R2.filtered.fastq.gz" \
    -w "$FASTP_THREADS" \
    --disable_adapter_trimming \
    -h "$run_output/fastp.html" \
    -j "$run_output/fastp.json" > "$run_output/fastp.log" 2>&1

echo "[dbitm] fastp finished successfully"
echo "[dbitm] fastp log: $final_dir/fastp/fastp.log"

echo "[dbitm] fastp result: $final_dir/fastp"
echo "====== dbitm fastp finished ======"
