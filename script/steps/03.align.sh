#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 )); then
    echo "Usage: 03.align.sh <assay> <raw_fastq_folder> [--chunk NNNN] [--dry-run]" >&2
    exit 1
fi
assay=$1
raw_path=$2
shift 2
dry_run=false
selected_chunk=""
while (( $# > 0 )); do
    case "$1" in
        --chunk)
            if (( $# < 2 )); then
                echo "[dbitm] align: --chunk requires a chunk number" >&2
                exit 1
            fi
            selected_chunk=$2
            shift 2
            ;;
        --dry-run)
            dry_run=true
            shift
            ;;
        *)
            echo "[dbitm] align: unknown argument: $1" >&2
            exit 1
            ;;
    esac
done
if [[ -n "$selected_chunk" ]]; then
    if [[ ! "$selected_chunk" =~ ^[0-9]+$ ]] || (( 10#$selected_chunk < 1 )); then
        echo "[dbitm] align: --chunk must be a positive integer" >&2
        exit 1
    fi
    printf -v selected_chunk '%04d' "$((10#$selected_chunk))"
fi
case "$assay" in
    taps|taps-v2|emseq|cabernet|smc) ;;
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
biscuit_directional_mode=$BISCUIT_DIRECTIONAL_MODE
case "$assay" in
    taps|taps-v2)
        aligner=bwa
        reference=$BWA_INDEX
        if [[ -z "$reference" ]]; then
            echo "[dbitm] align: BWA_INDEX is required for assay '$assay'" >&2
            exit 1
        fi
        ;;
    emseq|cabernet)
        aligner=biscuit
        reference=$BISCUIT_REFERENCE
        if [[ -z "$reference" ]]; then
            echo "[dbitm] align: BISCUIT_REFERENCE is required for assay '$assay'" >&2
            exit 1
        fi
        ;;
    smc)
        aligner=biscuit
        reference=$BISCUIT_REFERENCE
        biscuit_directional_mode=1
        if [[ -z "$reference" ]]; then
            echo "[dbitm] align: BISCUIT_REFERENCE is required for assay '$assay'" >&2
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

run_input=$barcode_dir
run_output=$final_dir/align

echo "====== dbitm align ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] aligner: $aligner"
echo "[dbitm] reference: $reference"
echo "[dbitm] input directory: $barcode_dir"
echo "[dbitm] config: $config_file"
if [[ "$aligner" == biscuit ]]; then
    echo "[dbitm] biscuit directional mode: $biscuit_directional_mode"
fi
echo "[dbitm] output directory: $final_dir/align"
if [[ -n "$selected_chunk" ]]; then
    echo "[dbitm] selected chunk: $selected_chunk"
fi

if [[ "$dry_run" == true ]]; then
    declare -a dry_chunk_files=()
    if [[ -d "$barcode_dir" ]]; then
        shopt -s nullglob
        if [[ "$assay" == smc ]]; then
            dry_chunk_files=("$barcode_dir"/*.watson.short-genomic.fastq.gz)
        else
            dry_chunk_files=("$barcode_dir"/*.R1.demux.fastq.gz)
        fi
        shopt -u nullglob
    fi
    echo "[dbitm] dry-run: no files will be written"
    echo "[dbitm] discovered chunks: ${#dry_chunk_files[@]}"
    if [[ "$assay" == smc ]]; then
        echo "[dbitm] expected Watson input: $barcode_dir/*.watson.{short-genomic,genomic}.fastq.gz"
        echo "[dbitm] expected Crick input: $barcode_dir/*.crick.{short-genomic,genomic}.fastq.gz"
        echo "[dbitm] Watson mate assignment: R1/parent=genomic R2/daughter=short-genomic"
        echo "[dbitm] Crick mate assignment: R1/parent=short-genomic R2/daughter=genomic"
    else
        echo "[dbitm] expected input: $barcode_dir/*.R1.demux.fastq.gz"
    fi
    echo "[dbitm] planned command: $aligner + sinto nametotag -> $final_dir/align"
    echo "====== dbitm align dry-run finished ======"
    exit 0
fi

mkdir -p "$final_dir"
if [[ -z "$selected_chunk" ]]; then
    rm -rf -- "$run_output"
fi
mkdir -p "$run_output/logs"
if [[ -n "$selected_chunk" ]]; then
    align_log=$run_output/align.$selected_chunk.log
else
    align_log=$run_output/align.log
fi
: > "$align_log"

declare -a chunks=()
declare -a r1_files=()
declare -a r2_files=()
declare -a output_bams=()
declare -a chunk_logs=()
if [[ "$assay" == smc ]]; then
    declare -a crick_chunks=()
    declare -a crick_r1_files=()
    declare -a crick_r2_files=()
    declare -a crick_output_bams=()
    declare -a crick_chunk_logs=()
    if [[ -n "$selected_chunk" ]]; then
        watson_short_files=("$run_input/$selected_chunk.watson.short-genomic.fastq.gz")
    else
        shopt -s nullglob
        watson_short_files=("$run_input"/*.watson.short-genomic.fastq.gz)
        shopt -u nullglob
    fi
    if (( ${#watson_short_files[@]} == 0 )); then
        echo "[dbitm] align: no Watson FASTQ chunks found: $run_input/*.watson.short-genomic.fastq.gz" >&2
        exit 1
    fi
    chunk_count=${#watson_short_files[@]}
    for watson_short in "${watson_short_files[@]}"; do
        if [[ ! -f "$watson_short" ]]; then
            echo "[dbitm] align: Watson FASTQ chunk not found: $watson_short" >&2
            exit 1
        fi
        filename=$(basename "$watson_short")
        chunk=${filename%.watson.short-genomic.fastq.gz}
        watson_long=$run_input/$chunk.watson.genomic.fastq.gz
        crick_short=$run_input/$chunk.crick.short-genomic.fastq.gz
        crick_long=$run_input/$chunk.crick.genomic.fastq.gz
        for paired_fastq in "$watson_long" "$crick_short" "$crick_long"; do
            if [[ ! -f "$paired_fastq" ]]; then
                echo "[dbitm] align: paired SmC FASTQ not found for chunk '$chunk': $paired_fastq" >&2
                exit 1
            fi
        done

        chunks+=("$chunk.watson")
        r1_files+=("$watson_long")
        r2_files+=("$watson_short")
        output_bams+=("$run_output/$chunk.watson.cb.bam")
        chunk_logs+=("$run_output/logs/$chunk.watson.log")

        crick_chunks+=("$chunk.crick")
        crick_r1_files+=("$crick_short")
        crick_r2_files+=("$crick_long")
        crick_output_bams+=("$run_output/$chunk.crick.cb.bam")
        crick_chunk_logs+=("$run_output/logs/$chunk.crick.log")
    done
    chunks+=("${crick_chunks[@]}")
    r1_files+=("${crick_r1_files[@]}")
    r2_files+=("${crick_r2_files[@]}")
    output_bams+=("${crick_output_bams[@]}")
    chunk_logs+=("${crick_chunk_logs[@]}")
else
    if [[ -n "$selected_chunk" ]]; then
        r1_files=("$run_input/$selected_chunk.R1.demux.fastq.gz")
    else
        shopt -s nullglob
        r1_files=("$run_input"/*.R1.demux.fastq.gz)
        shopt -u nullglob
    fi
    if (( ${#r1_files[@]} == 0 )); then
        echo "[dbitm] align: no R1 demux FASTQ chunks found: $run_input/*.R1.demux.fastq.gz" >&2
        exit 1
    fi
    chunk_count=${#r1_files[@]}
    for r1 in "${r1_files[@]}"; do
        if [[ ! -f "$r1" ]]; then
            echo "[dbitm] align: R1 demux FASTQ chunk not found: $r1" >&2
            exit 1
        fi
        filename=$(basename "$r1")
        chunk=${filename%.R1.demux.fastq.gz}
        r2=$run_input/$chunk.R2.demux.fastq.gz
        if [[ ! -f "$r2" ]]; then
            echo "[dbitm] align: paired R2 FASTQ not found for chunk '$chunk': $r2" >&2
            exit 1
        fi

        chunks+=("$chunk")
        r2_files+=("$r2")
        output_bams+=("$run_output/$chunk.cb.bam")
        chunk_logs+=("$run_output/logs/$chunk.log")
    done
fi

align_threads_per_job=$ALIGN_THREADS_PER_CHUNK

echo "[dbitm] chunks: $chunk_count"
echo "[dbitm] alignments: ${#chunks[@]}"
echo "[dbitm] threads per alignment: $align_threads_per_job"

align_chunk() {
    local chunk=$1
    local r1=$2
    local r2=$3
    local output_bam=$4
    local chunk_log=$5
    echo "[dbitm] aligning chunk: $chunk"
    echo "[$chunk] aligner=$aligner threads=$align_threads_per_job r1=$r1 r2=$r2 output=$output_bam" > "$chunk_log"
    if [[ "$assay" == smc ]]; then
        echo "[$chunk] r1-role=parent r2-role=daughter biscuit-directional-mode=1" >> "$chunk_log"
    fi
    if [[ "$aligner" == bwa ]]; then
        pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
            bwa mem -t "$align_threads_per_job" "$reference" "$r1" "$r2" \
            2>> "$chunk_log" \
            | pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
                sinto nametotag -b - -O b -o "$output_bam" \
                >> "$chunk_log" 2>&1
    else
        pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
            biscuit align -@ "$align_threads_per_job" -b "$biscuit_directional_mode" \
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

for chunk_index in "${!chunks[@]}"; do
    chunk_name=${chunks[$chunk_index]}
    trap 'job_status=$?; if (( job_status != 0 )); then echo "[dbitm] align: chunk failed: $chunk_name" >&2; fi' EXIT
    align_chunk \
        "$chunk_name" \
        "${r1_files[$chunk_index]}" \
        "${r2_files[$chunk_index]}" \
        "${output_bams[$chunk_index]}" \
        "${chunk_logs[$chunk_index]}"
    trap - EXIT
done

for chunk_log in "${chunk_logs[@]}"; do
    if [[ -f "$chunk_log" ]]; then
        printf '\n===== %s =====\n' "$(basename "$chunk_log")" >> "$align_log"
        cat "$chunk_log" >> "$align_log"
    fi
done

echo "[dbitm] align log: $align_log"
echo "[dbitm] align result: $final_dir/align"
echo "====== dbitm align finished ======"
