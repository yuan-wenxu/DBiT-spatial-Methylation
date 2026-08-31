#!/usr/bin/env bash
set -euo pipefail

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ ${use_scratch:-false} == true && -n ${scratch_run:-} && -d $scratch_run ]]; then
        if (( status != 0 )) && [[ -n ${run_output:-} && -d $run_output && -n ${final_dir:-} ]]; then
            echo "[dbitm] spike-align: recovering scratch results after exit status $status: $final_dir/spike_align" >&2
            if mkdir -p "$final_dir/spike_align" && cp -a "$run_output/." "$final_dir/spike_align/"; then
                echo "[dbitm] spike-align: scratch results recovered: $final_dir/spike_align" >&2
            else
                echo "[dbitm] spike-align: warning: failed to recover scratch results: $run_output" >&2
            fi
        fi
        rm -rf -- "$scratch_run" || echo "[dbitm] spike-align: warning: failed to clean scratch directory: $scratch_run" >&2
    fi
    exit "$status"
}

enable_cleanup() {
    trap cleanup_scratch EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM
    trap 'exit 129' HUP
}

require_index_directory() {
    local spike_name=$1
    local index_path=$2
    local index_dir
    index_dir=$(dirname "$index_path")

    if [[ ! -d "$index_dir" ]]; then
        echo "[dbitm] spike-align: index directory not found ($spike_name): $index_dir" >&2
        exit 1
    fi
    if [[ -z $(find "$index_dir" -mindepth 1 -maxdepth 1 -print -quit) ]]; then
        echo "[dbitm] spike-align: index directory is empty ($spike_name): $index_dir" >&2
        exit 1
    fi
}

if (( $# < 2 || $# > 3 )); then
    echo "Usage: 03.spike_align.sh <assay> <raw_fastq_folder> [--dry-run]" >&2
    exit 1
fi
assay=$1
raw_path=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] spike-align: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
case "$assay" in
    taps|taps-v2|emseq|cabernet) ;;
    *) echo "[dbitm] spike-align: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] spike-align: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] spike-align: config file not found: $config_file" >&2
    exit 1
fi
source "$config_file"

if [[ ! "$SPIKE_ALIGN_THREADS_PER_CHUNK" =~ ^[1-9][0-9]*$ ]]; then
    echo "[dbitm] spike-align: SPIKE_ALIGN_THREADS_PER_CHUNK must be >= 1" >&2
    exit 1
fi
if [[ "$BISCUIT_DIRECTIONAL_MODE" != 0 && "$BISCUIT_DIRECTIONAL_MODE" != 1 ]]; then
    echo "[dbitm] spike-align: BISCUIT_DIRECTIONAL_MODE must be 0 or 1" >&2
    exit 1
fi

aligner=""
index_variable=""
case "$assay" in
    taps|taps-v2)
        aligner=bwa
        index_variable=BWA_SPIKE_IN_INDEXES
        ;;
    emseq|cabernet)
        aligner=biscuit
        index_variable=BISCUIT_SPIKE_IN_INDEXES
        ;;
esac

index_declaration=$(declare -p "$index_variable" 2>/dev/null || true)
if [[ "$index_declaration" != "declare -A "* ]]; then
    echo "[dbitm] spike-align: $index_variable must be defined with declare -A" >&2
    exit 1
fi
declare -n configured_indexes=$index_variable
mapfile -t configured_names < <(
    printf '%s\n' "${!configured_indexes[@]}" | LC_ALL=C sort
)

declare -a spike_names=()
declare -a spike_references=()
for spike_name in "${configured_names[@]}"; do
    spike_reference=${configured_indexes[$spike_name]}
    spike_reference=${spike_reference#"${spike_reference%%[![:space:]]*}"}
    spike_reference=${spike_reference%"${spike_reference##*[![:space:]]}"}
    if [[ -z "$spike_reference" ]]; then
        continue
    fi
    if [[ ! "$spike_name" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        echo "[dbitm] spike-align: invalid spike-in name: $spike_name" >&2
        exit 1
    fi
    if [[ "$spike_reference" != /* ]]; then
        echo "[dbitm] spike-align: spike-in index path must be absolute ($spike_name): $spike_reference" >&2
        exit 1
    fi
    spike_reference=$(realpath -m "$spike_reference")
    require_index_directory "$spike_name" "$spike_reference"
    spike_names+=("$spike_name")
    spike_references+=("$spike_reference")
done

if (( ${#spike_names[@]} == 0 )); then
    echo "[dbitm] spike-align: no indexes configured in $index_variable" >&2
    exit 1
fi

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
barcode_dir=$final_dir/barcode
if [[ ! -d "$barcode_dir" ]]; then
    if [[ "$dry_run" == true ]]; then
        echo "[dbitm] spike-align: dry-run expects future barcode output directory: $barcode_dir"
    else
        echo "[dbitm] spike-align: barcode output directory not found: $barcode_dir" >&2
        exit 1
    fi
fi

use_scratch=false
run_input=$barcode_dir
run_output=$final_dir/spike_align

echo "====== dbitm spike-align ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] aligner: $aligner"
echo "[dbitm] configured spike-ins: ${#spike_names[@]}"
for spike_index in "${!spike_names[@]}"; do
    echo "[dbitm] spike-in ${spike_names[$spike_index]}: ${spike_references[$spike_index]}"
done
echo "[dbitm] input directory: $barcode_dir"
echo "[dbitm] config: $config_file"
if [[ "$aligner" == biscuit ]]; then
    echo "[dbitm] biscuit directional mode: $BISCUIT_DIRECTIONAL_MODE"
fi
echo "[dbitm] output directory: $final_dir/spike_align"

if [[ "$dry_run" == true ]]; then
    if [[ -n ${SCRATCH_ROOT:-} && "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] spike-align: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    declare -a dry_r1_files=()
    if [[ -d "$barcode_dir" ]]; then
        shopt -s nullglob
        dry_r1_files=("$barcode_dir"/*.R1.spike-in.fastq.gz)
        shopt -u nullglob
    fi
    echo "[dbitm] dry-run: no files will be written"
    echo "[dbitm] discovered chunks: ${#dry_r1_files[@]}"
    echo "[dbitm] expected input: $barcode_dir/*.R1.spike-in.fastq.gz"
    echo "[dbitm] planned command: $aligner -> $final_dir/spike_align"
    [[ -z ${SCRATCH_ROOT:-} ]] || echo "[dbitm] planned scratch root: $SCRATCH_ROOT"
    echo "====== dbitm spike-align dry-run finished ======"
    exit 0
fi

if [[ -n ${SCRATCH_ROOT:-} ]]; then
    if [[ "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] spike-align: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    echo "[dbitm] scratch root: $SCRATCH_ROOT"
    mkdir -p "$SCRATCH_ROOT"
    scratch_root=$(realpath "$SCRATCH_ROOT")
    run_id=${SLURM_JOB_ID:-spike_align_$$}
    scratch_run=$scratch_root/dbitm/$run_id
    run_input=$scratch_run/barcode
    run_output=$scratch_run/spike_align
    use_scratch=true
    enable_cleanup
    mkdir -p "$run_input"
    shopt -s nullglob
    scratch_fastq_files=(
        "$barcode_dir"/*.R1.spike-in.fastq.gz
        "$barcode_dir"/*.R2.spike-in.fastq.gz
    )
    shopt -u nullglob
    echo "[dbitm] copying ${#scratch_fastq_files[@]} spike-in FASTQ files to scratch: $run_input"
    if (( ${#scratch_fastq_files[@]} > 0 )); then
        cp -a "${scratch_fastq_files[@]}" "$run_input"/
    fi
else
    mkdir -p "$final_dir"
fi

rm -rf -- "$run_output"
mkdir -p "$run_output/logs"
align_log=$run_output/spike_align.log
: > "$align_log"

shopt -s nullglob
r1_files=("$run_input"/*.R1.spike-in.fastq.gz)
shopt -u nullglob
if (( ${#r1_files[@]} == 0 )); then
    echo "[dbitm] spike-align: no R1 spike-in FASTQ chunks found: $run_input/*.R1.spike-in.fastq.gz" >&2
    exit 1
fi

chunk_count=${#r1_files[@]}
parallel_jobs=$chunk_count
align_threads_per_job=$SPIKE_ALIGN_THREADS_PER_CHUNK
total_align_threads=$((chunk_count * align_threads_per_job))

echo "[dbitm] chunks: $chunk_count"
echo "[dbitm] parallel jobs: $parallel_jobs"
echo "[dbitm] threads per chunk: $align_threads_per_job"
echo "[dbitm] total aligner threads: $total_align_threads"

declare -a chunks=()
declare -a r2_files=()
declare -a chunk_logs=()
for r1 in "${r1_files[@]}"; do
    filename=$(basename "$r1")
    chunk=${filename%.R1.spike-in.fastq.gz}
    r2=$run_input/$chunk.R2.spike-in.fastq.gz
    if [[ ! -f "$r2" ]]; then
        echo "[dbitm] spike-align: paired R2 FASTQ not found for chunk '$chunk': $r2" >&2
        exit 1
    fi
    chunks+=("$chunk")
    r2_files+=("$r2")
    chunk_logs+=("$run_output/logs/$chunk.log")
done

align_spike_chunk() {
    local chunk=$1
    local r1=$2
    local r2=$3
    local chunk_log=$4
    local spike_index
    local spike_name
    local spike_reference
    local output_bam
    local flagstat_report

    echo "[dbitm] aligning spike-in chunk: $chunk"
    echo "[$chunk] aligner=$aligner threads=$align_threads_per_job r1=$r1 r2=$r2" > "$chunk_log"
    for spike_index in "${!spike_names[@]}"; do
        spike_name=${spike_names[$spike_index]}
        spike_reference=${spike_references[$spike_index]}
        output_bam=$run_output/$chunk.$spike_name.bam
        flagstat_report=$run_output/$chunk.$spike_name.flagstat.txt
        echo "[dbitm] mapping chunk $chunk to: $spike_name"
        echo "[$chunk] spike_in=$spike_name reference=$spike_reference output=$output_bam" >> "$chunk_log"
        if [[ "$aligner" == bwa ]]; then
            pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
                bwa mem -t "$align_threads_per_job" "$spike_reference" "$r1" "$r2" \
                2>> "$chunk_log" \
                | pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
                    samtools view -b -o "$output_bam" - \
                    >> "$chunk_log" 2>&1
        else
            pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
                biscuit align -@ "$align_threads_per_job" -b "$BISCUIT_DIRECTIONAL_MODE" \
                    "$spike_reference" "$r1" "$r2" \
                2>> "$chunk_log" \
                | pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
                    samtools view -b -o "$output_bam" - \
                    >> "$chunk_log" 2>&1
        fi
        if ! pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
            samtools quickcheck -v "$output_bam" >> "$chunk_log" 2>&1; then
            echo "[dbitm] spike-align: BAM quickcheck failed: $output_bam" >&2
            return 1
        fi
        pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
            samtools flagstat -@ "$align_threads_per_job" "$output_bam" \
            > "$flagstat_report" 2>> "$chunk_log"
        echo "[$chunk] flagstat=$flagstat_report" >> "$chunk_log"
    done
    echo "[dbitm] spike-in chunk finished: $chunk"
}

declare -a active_pids=()
declare -a active_chunks=()
alignment_failed=false

wait_for_active_jobs() {
    local job_index
    for job_index in "${!active_pids[@]}"; do
        if ! wait "${active_pids[$job_index]}"; then
            echo "[dbitm] spike-align: chunk failed: ${active_chunks[$job_index]}" >&2
            alignment_failed=true
        fi
    done
    active_pids=()
    active_chunks=()
}

for chunk_index in "${!chunks[@]}"; do
    align_spike_chunk \
        "${chunks[$chunk_index]}" \
        "${r1_files[$chunk_index]}" \
        "${r2_files[$chunk_index]}" \
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
    echo "[dbitm] spike-align: one or more chunks failed; see logs under: $run_output/logs" >&2
    exit 1
fi

if [[ "$use_scratch" == true ]]; then
    rm -rf -- "$final_dir/spike_align"
    mkdir -p "$final_dir"
    echo "[dbitm] copying spike-in alignment result from scratch to $final_dir/spike_align"
    cp -a "$run_output" "$final_dir/spike_align"
fi
echo "[dbitm] spike-align log: $final_dir/spike_align/spike_align.log"
echo "[dbitm] spike-align result: $final_dir/spike_align"
echo "====== dbitm spike-align finished ======"
