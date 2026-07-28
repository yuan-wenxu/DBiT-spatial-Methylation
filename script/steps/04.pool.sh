#!/usr/bin/env bash
set -euo pipefail

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ ${use_scratch:-false} == true && -n ${scratch_run:-} && -d $scratch_run ]]; then
        rm -rf -- "$scratch_run" || echo "[dbitm] pool: warning: failed to clean scratch directory: $scratch_run" >&2
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
    echo "Usage: 04.pool.sh <assay> <raw_fastq_folder>" >&2
    exit 1
fi
assay=$1
raw_path=$2
case "$assay" in
    taps|taps-v2|emseq) ;;
    *) echo "[dbitm] pool: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] pool: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] pool: config file not found: $config_file" >&2
    exit 1
fi
source "$config_file"

# derive per-thread sort memory: POOL_MEM / POOL_THREADS (fallback to POOL_SORT_MEM)
if [[ -n "${POOL_SORT_MEM:-}" ]]; then
    sort_mem=$POOL_SORT_MEM
else
    pool_mem_num=$(echo "${POOL_MEM:-64G}" | sed 's/[^0-9]//g')
    sort_mem=$((pool_mem_num / ${POOL_THREADS:-4}))G
fi

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
align_dir=$final_dir/align
if [[ ! -d "$align_dir" ]]; then
    echo "[dbitm] pool: align output directory not found: $align_dir" >&2
    exit 1
fi

shopt -s nullglob
cb_bams=("$align_dir"/*.cb.bam)
shopt -u nullglob
if (( ${#cb_bams[@]} == 0 )); then
    echo "[dbitm] pool: no .cb.bam files found under: $align_dir" >&2
    exit 1
fi

use_scratch=false
run_input=$align_dir
run_output=$final_dir/pooled

echo "====== dbitm pool ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] input chunks: ${#cb_bams[@]}"
echo "[dbitm] align directory: $align_dir"
echo "[dbitm] config: $config_file"
echo "[dbitm] output directory: $final_dir/pooled"

if [[ -n ${SCRATCH_ROOT:-} ]]; then
    if [[ "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] pool: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    echo "[dbitm] scratch root: $SCRATCH_ROOT"
    mkdir -p "$SCRATCH_ROOT"
    scratch_root=$(realpath "$SCRATCH_ROOT")
    run_id=${SLURM_JOB_ID:-pool_$$}
    scratch_run=$scratch_root/dbitm/$run_id
    run_input=$scratch_run/align
    run_output=$scratch_run/pooled
    use_scratch=true
    enable_cleanup
    mkdir -p "$run_input"
    echo "[dbitm] copying align BAM files to scratch: $run_input"
    cp -a "$align_dir/." "$run_input"/
fi

rm -rf -- "$run_output"
mkdir -p "$run_output"

cat_tmp=$run_output/unsorted.cb.bam
pooled_bam=$run_output/pooled.byCB.bam
pool_log=$run_output/pool.log
: > "$pool_log"

echo "[dbitm] concatenating BAM chunks..."
pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
    samtools cat -o "$cat_tmp" "$run_input"/*.cb.bam \
    >> "$pool_log" 2>&1

echo "[dbitm] sorting by genomic coordinate..."
pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
    samtools sort \
    -m "$sort_mem" \
    -@ "${POOL_THREADS:-4}" \
    -o "$pooled_bam" "$cat_tmp" \
    >> "$pool_log" 2>&1

rm -f -- "$cat_tmp"

echo "[dbitm] indexing pooled BAM..."
pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
    samtools index "$pooled_bam" \
    >> "$pool_log" 2>&1

echo "[dbitm] pool finished: $pooled_bam"

if [[ "$use_scratch" == true ]]; then
    rm -rf -- "$final_dir/pooled"
    mkdir -p "$final_dir"
    echo "[dbitm] copying pooled result from scratch to $final_dir/pooled"
    cp -a "$run_output" "$final_dir/pooled"
fi
echo "[dbitm] pool log: $final_dir/pooled/pool.log"
echo "[dbitm] pool result: $final_dir/pooled"
echo "====== dbitm pool finished ======"
