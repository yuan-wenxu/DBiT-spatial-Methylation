#!/usr/bin/env bash
set -euo pipefail

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ ${use_scratch:-false} == true && -n ${scratch_run:-} && -d $scratch_run ]]; then
        rm -rf -- "$scratch_run" || echo "[dbitm] align: warning: failed to clean scratch directory: $scratch_run" >&2
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
    echo "Usage: 03.align.sh <assay> <raw_fastq_folder> [--dry-run]" >&2
    exit 1
fi
assay=$1
raw_path=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] align: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
case "$assay" in
    taps|taps-v2|emseq) ;;
    *) echo "[dbitm] align: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] align: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] align: config file not found: $config_file" >&2
    exit 1
fi
source "$config_file"

# Keep older local config files usable after this step is added.
BWA_INDEX=${BWA_INDEX:-}
BISCUIT_REFERENCE=${BISCUIT_REFERENCE:-}
BISCUIT_DIRECTIONAL_MODE=${BISCUIT_DIRECTIONAL_MODE:-1}
ALIGN_THREADS_PER_CHUNK=${ALIGN_THREADS_PER_CHUNK:-8}

if [[ ! "$ALIGN_THREADS_PER_CHUNK" =~ ^[1-9][0-9]*$ ]]; then
    echo "[dbitm] align: ALIGN_THREADS_PER_CHUNK must be >= 1" >&2
    exit 1
fi
if [[ "$BISCUIT_DIRECTIONAL_MODE" != 0 && "$BISCUIT_DIRECTIONAL_MODE" != 1 ]]; then
    echo "[dbitm] align: BISCUIT_DIRECTIONAL_MODE must be 0 or 1" >&2
    exit 1
fi

aligner=""
reference=""
case "$assay" in
    taps|taps-v2)
        aligner=bwa
        reference=$BWA_INDEX
        if [[ -z "$reference" ]]; then
            echo "[dbitm] align: BWA_INDEX is required for assay '$assay'" >&2
            exit 1
        fi
        ;;
    emseq)
        aligner=biscuit
        reference=$BISCUIT_REFERENCE
        if [[ -z "$reference" ]]; then
            echo "[dbitm] align: BISCUIT_REFERENCE is required for assay 'emseq'" >&2
            exit 1
        fi
        ;;
esac
if [[ "$reference" != /* ]]; then
    reference=$REPO_DIR/$reference
fi
reference=$(realpath -m "$reference")
if [[ "$aligner" == bwa ]]; then
    if [[ ! -e "$reference" && ! -e "$reference.bwt" && ! -e "$reference.bwt.2bit.64" ]]; then
        echo "[dbitm] align: bwa reference/index not found: $reference" >&2
        exit 1
    fi
elif [[ ! -f "$reference" ]]; then
    echo "[dbitm] align: biscuit reference FASTA not found: $reference" >&2
    exit 1
fi

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
barcode_dir=$final_dir/barcode
if [[ ! -d "$barcode_dir" ]]; then
    if [[ "$dry_run" == true ]]; then
        echo "[dbitm] align: dry-run expects future barcode output directory: $barcode_dir"
    else
        echo "[dbitm] align: barcode output directory not found: $barcode_dir" >&2
        exit 1
    fi
fi

use_scratch=false
run_input=$barcode_dir
run_output=$final_dir/align

echo "====== dbitm align ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] aligner: $aligner"
echo "[dbitm] reference: $reference"
echo "[dbitm] input directory: $barcode_dir"
echo "[dbitm] config: $config_file"
if [[ "$aligner" == biscuit ]]; then
    echo "[dbitm] biscuit directional mode: $BISCUIT_DIRECTIONAL_MODE"
fi
echo "[dbitm] output directory: $final_dir/align"

if [[ "$dry_run" == true ]]; then
    if [[ -n ${SCRATCH_ROOT:-} && "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] align: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    declare -a dry_r1_files=()
    if [[ -d "$barcode_dir" ]]; then
        shopt -s nullglob
        dry_r1_files=("$barcode_dir"/*.R1.demux.fastq.gz)
        shopt -u nullglob
    fi
    echo "[dbitm] dry-run: no files will be written"
    echo "[dbitm] discovered chunks: ${#dry_r1_files[@]}"
    echo "[dbitm] expected input: $barcode_dir/*.R1.demux.fastq.gz"
    echo "[dbitm] planned command: $aligner + sinto nametotag -> $final_dir/align"
    [[ -z ${SCRATCH_ROOT:-} ]] || echo "[dbitm] planned scratch root: $SCRATCH_ROOT"
    echo "====== dbitm align dry-run finished ======"
    exit 0
fi

if [[ -n ${SCRATCH_ROOT:-} ]]; then
    if [[ "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] align: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    echo "[dbitm] scratch root: $SCRATCH_ROOT"
    mkdir -p "$SCRATCH_ROOT"
    scratch_root=$(realpath "$SCRATCH_ROOT")
    run_id=${SLURM_JOB_ID:-align_$$}
    scratch_run=$scratch_root/dbitm/$run_id
    run_input=$scratch_run/barcode
    run_output=$scratch_run/align
    use_scratch=true
    enable_cleanup
    mkdir -p "$run_input"
    echo "[dbitm] copying barcode FASTQ files to scratch: $run_input"
    cp -a "$barcode_dir/." "$run_input"/
else
    mkdir -p "$final_dir"
fi

rm -rf -- "$run_output"
mkdir -p "$run_output/logs"
align_log=$run_output/align.log
: > "$align_log"

shopt -s nullglob
r1_files=("$run_input"/*.R1.demux.fastq.gz)
shopt -u nullglob
if (( ${#r1_files[@]} == 0 )); then
    echo "[dbitm] align: no R1 demux FASTQ chunks found: $run_input/*.R1.demux.fastq.gz" >&2
    exit 1
fi

chunk_count=${#r1_files[@]}
parallel_jobs=$chunk_count
align_threads_per_job=$ALIGN_THREADS_PER_CHUNK
total_align_threads=$((chunk_count * align_threads_per_job))

echo "[dbitm] chunks: $chunk_count"
echo "[dbitm] parallel jobs: $parallel_jobs"
echo "[dbitm] threads per chunk: $align_threads_per_job"
echo "[dbitm] total aligner threads: $total_align_threads"

declare -a chunks=()
declare -a r2_files=()
declare -a output_bams=()
declare -a chunk_logs=()
for r1 in "${r1_files[@]}"; do
    filename=$(basename "$r1")
    chunk=${filename%.R1.demux.fastq.gz}
    r2=$run_input/$chunk.R2.demux.fastq.gz
    output_bam=$run_output/$chunk.cb.bam
    if [[ ! -f "$r2" ]]; then
        echo "[dbitm] align: paired R2 FASTQ not found for chunk '$chunk': $r2" >&2
        exit 1
    fi

    chunks+=("$chunk")
    r2_files+=("$r2")
    output_bams+=("$output_bam")
    chunk_logs+=("$run_output/logs/$chunk.log")
done

align_chunk() {
    local chunk=$1
    local r1=$2
    local r2=$3
    local output_bam=$4
    local chunk_log=$5
    echo "[dbitm] aligning chunk: $chunk"
    echo "[$chunk] aligner=$aligner threads=$align_threads_per_job r1=$r1 r2=$r2 output=$output_bam" > "$chunk_log"
    if [[ "$aligner" == bwa ]]; then
        pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
            bwa mem -t "$align_threads_per_job" "$reference" "$r1" "$r2" \
            2>> "$chunk_log" \
            | pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
                sinto nametotag -b - -O b -o "$output_bam" \
                >> "$chunk_log" 2>&1
    else
        pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
            biscuit align -@ "$align_threads_per_job" -b "$BISCUIT_DIRECTIONAL_MODE" \
                "$reference" "$r1" "$r2" \
            2>> "$chunk_log" \
            | pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
                sinto nametotag -b - -O b -o "$output_bam" \
                >> "$chunk_log" 2>&1
    fi
    pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        samtools quickcheck "$output_bam" >> "$chunk_log" 2>&1
    echo "[dbitm] chunk finished: $chunk"
}

declare -a active_pids=()
declare -a active_chunks=()
alignment_failed=false

wait_for_active_jobs() {
    local job_index
    for job_index in "${!active_pids[@]}"; do
        if ! wait "${active_pids[$job_index]}"; then
            echo "[dbitm] align: chunk failed: ${active_chunks[$job_index]}" >&2
            alignment_failed=true
        fi
    done
    active_pids=()
    active_chunks=()
}

for chunk_index in "${!chunks[@]}"; do
    align_chunk \
        "${chunks[$chunk_index]}" \
        "${r1_files[$chunk_index]}" \
        "${r2_files[$chunk_index]}" \
        "${output_bams[$chunk_index]}" \
        "${chunk_logs[$chunk_index]}" &
    active_pids+=("$!")
    active_chunks+=("${chunks[$chunk_index]}")

    if (( ${#active_pids[@]} == parallel_jobs )); then
        wait_for_active_jobs
        [[ "$alignment_failed" == false ]] || break
    fi
done
if (( ${#active_pids[@]} > 0 )); then
    wait_for_active_jobs
fi

for chunk_log in "${chunk_logs[@]}"; do
    if [[ -f "$chunk_log" ]]; then
        printf '\n===== %s =====\n' "$(basename "$chunk_log")" >> "$align_log"
        cat "$chunk_log" >> "$align_log"
    fi
done
if [[ "$alignment_failed" == true ]]; then
    echo "[dbitm] align: one or more chunks failed; see per-chunk logs under: $run_output/logs" >&2
    exit 1
fi

if [[ "$use_scratch" == true ]]; then
    rm -rf -- "$final_dir/align"
    mkdir -p "$final_dir"
    echo "[dbitm] copying alignment result from scratch to $final_dir/align"
    cp -a "$run_output" "$final_dir/align"
fi
echo "[dbitm] align log: $final_dir/align/align.log"
echo "[dbitm] align result: $final_dir/align"
echo "====== dbitm align finished ======"
