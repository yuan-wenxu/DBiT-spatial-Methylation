#!/usr/bin/env bash
set -euo pipefail

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ ${use_scratch:-false} == true && -n ${scratch_run:-} && -d $scratch_run ]]; then
        if (( status != 0 )) && [[ -n ${run_output:-} && -d $run_output && -n ${final_dir:-} ]]; then
            echo "[dbitm] pool: recovering scratch results after exit status $status: $final_dir/pooled" >&2
            if mkdir -p "$final_dir/pooled" && cp -a "$run_output/." "$final_dir/pooled/"; then
                echo "[dbitm] pool: scratch results recovered: $final_dir/pooled" >&2
            else
                echo "[dbitm] pool: warning: failed to recover scratch results: $run_output" >&2
            fi
        fi
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

if (( $# < 2 || $# > 3 )); then
    echo "Usage: 04.pool.sh <assay> <raw_fastq_folder> [--dry-run]" >&2
    exit 1
fi
assay=$1
raw_path=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] pool: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
case "$assay" in
    taps|taps-v2|emseq|cabernet|smc) ;;
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

# use configured per-thread sort memory, otherwise derive POOL_MEM / POOL_THREADS
if [[ -n "${POOL_SORT_MEM}" ]]; then
    sort_mem=$POOL_SORT_MEM
else
    pool_mem_num=$(echo "${POOL_MEM}" | sed 's/[^0-9]//g')
    sort_mem=$((pool_mem_num / ${POOL_THREADS}))G
fi

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
align_dir=$final_dir/align
spike_align_dir=$final_dir/spike_align
if [[ ! -d "$align_dir" ]]; then
    if [[ "$dry_run" == true ]]; then
        echo "[dbitm] pool: dry-run expects future align output directory: $align_dir"
    else
        echo "[dbitm] pool: align output directory not found: $align_dir" >&2
        exit 1
    fi
fi

shopt -s nullglob
cb_bams=("$align_dir"/*.cb.bam)
shopt -u nullglob
if (( ${#cb_bams[@]} == 0 )); then
    if [[ "$dry_run" == true ]]; then
        echo "[dbitm] pool: dry-run found no existing .cb.bam files under: $align_dir"
    else
        echo "[dbitm] pool: no .cb.bam files found under: $align_dir" >&2
        exit 1
    fi
fi

declare -a spike_names=(lambda puc19)
declare -a source_spike_bams=()
shopt -s nullglob
for spike_name in "${spike_names[@]}"; do
    source_spike_bams+=("$spike_align_dir"/*."$spike_name".bam)
done
shopt -u nullglob

use_scratch=false
run_input=$align_dir
run_spike_input=$spike_align_dir
run_output=$final_dir/pooled

echo "====== dbitm pool ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] input chunks: ${#cb_bams[@]}"
echo "[dbitm] align directory: $align_dir"
echo "[dbitm] spike-in BAM files found: ${#source_spike_bams[@]}"
echo "[dbitm] config: $config_file"
echo "[dbitm] output directory: $final_dir/pooled"

if [[ "$dry_run" == true ]]; then
    if [[ -n ${SCRATCH_ROOT:-} && "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] pool: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    echo "[dbitm] dry-run: no files will be written"
    echo "[dbitm] sort memory per thread: $sort_mem"
    echo "[dbitm] planned command: samtools cat/sort/index -> $final_dir/pooled"
    [[ -z ${SCRATCH_ROOT:-} ]] || echo "[dbitm] planned scratch root: $SCRATCH_ROOT"
    echo "====== dbitm pool dry-run finished ======"
    exit 0
fi

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
    run_spike_input=$scratch_run/spike_align
    run_output=$scratch_run/pooled
    use_scratch=true
    enable_cleanup
    mkdir -p "$run_input"
    echo "[dbitm] copying align BAM files to scratch: $run_input"
    cp -a "$align_dir/." "$run_input"/
    if (( ${#source_spike_bams[@]} > 0 )); then
        mkdir -p "$run_spike_input"
        echo "[dbitm] copying spike-in BAM files to scratch: $run_spike_input"
        cp -a "${source_spike_bams[@]}" "$run_spike_input/"
    fi
fi

rm -rf -- "$run_output"
mkdir -p "$run_output"

pooled_bam=$run_output/pooled.cb.bam
pool_log=$run_output/pool.log
: > "$pool_log"

pool_bam_set() {
    local pool_name=$1
    local output_bam=$2
    shift 2
    local input_bams=("$@")
    local cat_tmp=$run_output/unsorted.$pool_name.bam

    echo "[dbitm] concatenating $pool_name BAM chunks..."
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        samtools cat -o "$cat_tmp" "${input_bams[@]}" \
        >> "$pool_log" 2>&1

    echo "[dbitm] sorting $pool_name BAM by genomic coordinate..."
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        samtools sort \
        -m "$sort_mem" \
        -@ "${POOL_THREADS}" \
        -o "$output_bam" "$cat_tmp" \
        >> "$pool_log" 2>&1

    rm -f -- "$cat_tmp"

    echo "[dbitm] indexing $pool_name pooled BAM..."
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        samtools index "$output_bam" \
        >> "$pool_log" 2>&1
    echo "[dbitm] pool finished: $output_bam"
}

shopt -s nullglob
declare -a host_bams=("$run_input"/*.cb.bam)
shopt -u nullglob
pool_bam_set cb "$pooled_bam" "${host_bams[@]}"

for spike_name in "${spike_names[@]}"; do
    shopt -s nullglob
    spike_bams=("$run_spike_input"/*."$spike_name".bam)
    shopt -u nullglob
    if (( ${#spike_bams[@]} > 0 )); then
        pool_bam_set "$spike_name" "$run_output/pooled.$spike_name.bam" "${spike_bams[@]}"
    fi
done

if [[ "$use_scratch" == true ]]; then
    rm -rf -- "$final_dir/pooled"
    mkdir -p "$final_dir"
    echo "[dbitm] copying pooled result from scratch to $final_dir/pooled"
    cp -a "$run_output" "$final_dir/pooled"
fi
echo "[dbitm] pool log: $final_dir/pooled/pool.log"
echo "[dbitm] pool result: $final_dir/pooled"
echo "====== dbitm pool finished ======"
