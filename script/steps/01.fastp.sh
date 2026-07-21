#!/usr/bin/env bash
set -euo pipefail

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ ${use_scratch:-false} == true && -n ${scratch_run:-} && -d $scratch_run ]]; then
        rm -rf -- "$scratch_run" || echo "fastp: warning: failed to clean scratch directory: $scratch_run" >&2
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
    echo "Usage: fastp.sh <assay> <raw_fastq_folder>" >&2; exit 1;
fi
assay=$1
raw_path=$2
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
use_scratch=false
run_input=$raw_abs
run_output=$final_dir/fastp

echo "====== dbitm fastp ======"
echo "[dbitm] input directory: $raw_abs"
echo "[dbitm] config: ${config_file}"
echo "[dbitm] threads: $FASTP_THREADS"
echo "[dbitm] output directory: $final_dir/fastp"

if [[ -n ${SCRATCH_ROOT:-} ]]; then
    if [[ "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] fastp: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    echo "[dbitm] scratch root: $SCRATCH_ROOT"
    mkdir -p "$SCRATCH_ROOT"
    scratch_root=$(realpath "$SCRATCH_ROOT")
    run_id=${SLURM_JOB_ID:-fastp_$$}
    scratch_run=$scratch_root/dbitm/$run_id
    run_input=$scratch_run/input
    run_output=$scratch_run/fastp
    use_scratch=true
    enable_cleanup
    mkdir -p "$run_input" "$run_output"
    echo "[dbitm] copying input files to scratch: $run_input"
    cp -a "$raw_abs/." "$run_input"/
else
    mkdir -p "$final_dir"
fi
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

if [[ "$use_scratch" == true ]]; then
    mkdir -p "$final_dir/fastp"
    rm -rf -- "$final_dir/fastp"
    echo "[dbitm] copying fastp result from scratch to $final_dir/fastp"
    cp -a "$run_output" "$final_dir/fastp"
    echo "[dbitm] fastp result copied from scratch to $final_dir/fastp"
fi
echo "[dbitm] fastp result: $final_dir/fastp"
echo "====== dbitm fastp finished ======"
