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

echo "====== dbitm fastp ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] input directory: $raw_abs"
echo "[dbitm] config: ${config_file}"
echo "[dbitm] threads: $FASTP_THREADS"
echo "[dbitm] output directory: $final_dir/fastp"

if [[ "$dry_run" == true ]]; then
    declare -a dry_fastq_candidates=()
    declare -a dry_fastq_files=()
    shopt -s nullglob nocaseglob
    dry_fastq_candidates=(
        "$raw_abs"/*.fastq.gz
        "$raw_abs"/*.fq.gz
        "$raw_abs"/*.fastq
        "$raw_abs"/*.fq
    )
    shopt -u nullglob nocaseglob
    for path in "${dry_fastq_candidates[@]}"; do
        [[ -f "$path" ]] && dry_fastq_files+=("$path")
    done
    if (( ${#dry_fastq_files[@]} != 2 )); then
        echo "[dbitm] fastp: input directory must contain exactly two FASTQ files, found ${#dry_fastq_files[@]}: $raw_abs" >&2
        exit 1
    fi
    dry_r1=
    dry_r2=
    for path in "${dry_fastq_files[@]}"; do
        filename=$(basename "$path")
        if [[ "$filename" =~ (^|[^[:alnum:]])[Rr]1([^[:alnum:]]|$) ]]; then
            [[ -z "$dry_r1" ]] || { echo "[dbitm] fastp: input directory contains two R1 FASTQ files" >&2; exit 1; }
            dry_r1=$path
        elif [[ "$filename" =~ (^|[^[:alnum:]])[Rr]2([^[:alnum:]]|$) ]]; then
            [[ -z "$dry_r2" ]] || { echo "[dbitm] fastp: input directory contains two R2 FASTQ files" >&2; exit 1; }
            dry_r2=$path
        else
            echo "[dbitm] fastp: FASTQ filename must contain an independent R1 or R2 token: $filename" >&2
            exit 1
        fi
    done
    [[ -n "$dry_r1" && -n "$dry_r2" ]] || { echo "[dbitm] fastp: input directory must contain one R1 and one R2 FASTQ file" >&2; exit 1; }
    echo "[dbitm] dry-run: no files will be written"
    echo "[dbitm] input R1: $dry_r1"
    echo "[dbitm] input R2: $dry_r2"
    echo "[dbitm] planned command: fastp -> $final_dir/fastp"
    echo "====== dbitm fastp dry-run finished ======"
    exit 0
fi

mkdir -p "$final_dir"
rm -rf -- "$run_output"
mkdir -p "$run_output"

declare -a fastq_candidates=()
declare -a fastq_files=()
shopt -s nullglob nocaseglob
fastq_candidates=(
    "$run_input"/*.fastq.gz
    "$run_input"/*.fq.gz
    "$run_input"/*.fastq
    "$run_input"/*.fq
)
shopt -u nullglob nocaseglob
for path in "${fastq_candidates[@]}"; do
    [[ -f "$path" ]] && fastq_files+=("$path")
done

if (( ${#fastq_files[@]} != 2 )); then
    echo "[dbitm] fastp: input directory must contain exactly two FASTQ files, found ${#fastq_files[@]}: $run_input" >&2
    exit 1
fi

r1=""
r2=""
for path in "${fastq_files[@]}"; do
    filename=$(basename "$path")
    if [[ "$filename" =~ (^|[^[:alnum:]])[Rr]1([^[:alnum:]]|$) ]]; then
        [[ -z "$r1" ]] || { echo "[dbitm] fastp: input directory contains two R1 FASTQ files" >&2; exit 1; }
        r1=$path
    elif [[ "$filename" =~ (^|[^[:alnum:]])[Rr]2([^[:alnum:]]|$) ]]; then
        [[ -z "$r2" ]] || { echo "[dbitm] fastp: input directory contains two R2 FASTQ files" >&2; exit 1; }
        r2=$path
    else
        echo "[dbitm] fastp: FASTQ filename must contain an independent R1 or R2 token: $filename" >&2
        exit 1
    fi
done
[[ -n "$r1" && -n "$r2" ]] || { echo "[dbitm] fastp: input directory must contain one R1 and one R2 FASTQ file" >&2; exit 1; }

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
