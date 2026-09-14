#!/usr/bin/env bash
set -euo pipefail

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

if (( $# < 2 )); then
    echo "Usage: 03.spike_align.sh <assay> <raw_fastq_folder> [--chunk NNNN] [--dry-run]" >&2
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
                echo "[dbitm] spike-align: --chunk requires a chunk number" >&2
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
            echo "[dbitm] spike-align: unknown argument: $1" >&2
            exit 1
            ;;
    esac
done
if [[ -n "$selected_chunk" ]]; then
    if [[ ! "$selected_chunk" =~ ^[0-9]+$ ]] || (( 10#$selected_chunk < 1 )); then
        echo "[dbitm] spike-align: --chunk must be a positive integer" >&2
        exit 1
    fi
    printf -v selected_chunk '%04d' "$((10#$selected_chunk))"
fi
case "$assay" in
    taps|taps-v2|emseq|cabernet|smc) ;;
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
biscuit_directional_mode=$BISCUIT_DIRECTIONAL_MODE
case "$assay" in
    taps|taps-v2)
        aligner=bwa
        index_variable=BWA_SPIKE_IN_INDEXES
        ;;
    emseq|cabernet)
        aligner=biscuit
        index_variable=BISCUIT_SPIKE_IN_INDEXES
        ;;
    smc)
        aligner=biscuit
        index_variable=BISCUIT_SPIKE_IN_INDEXES
        biscuit_directional_mode=1
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
    echo "[dbitm] biscuit directional mode: $biscuit_directional_mode"
fi
echo "[dbitm] output directory: $final_dir/spike_align"
if [[ -n "$selected_chunk" ]]; then
    echo "[dbitm] selected chunk: $selected_chunk"
fi

if [[ "$dry_run" == true ]]; then
    declare -a dry_r1_files=()
    if [[ -d "$barcode_dir" ]]; then
        shopt -s nullglob
        if [[ "$assay" == smc ]]; then
            dry_r1_files=("$barcode_dir"/*.watson.short-genomic.fastq.gz)
        else
            dry_r1_files=("$barcode_dir"/*.R1.spike-in.fastq.gz)
        fi
        shopt -u nullglob
    fi
    echo "[dbitm] dry-run: no files will be written"
    echo "[dbitm] discovered chunks: ${#dry_r1_files[@]}"
    if [[ "$assay" == smc ]]; then
        echo "[dbitm] expected Watson input: $barcode_dir/*.watson.{short-genomic,genomic}.fastq.gz"
        echo "[dbitm] expected Crick input: $barcode_dir/*.crick.{short-genomic,genomic}.fastq.gz"
        echo "[dbitm] Watson mate assignment: R1/parent=genomic R2/daughter=short-genomic"
        echo "[dbitm] Crick mate assignment: R1/parent=short-genomic R2/daughter=genomic"
    else
        echo "[dbitm] expected input: $barcode_dir/*.R1.spike-in.fastq.gz"
    fi
    echo "[dbitm] planned command: $aligner -> $final_dir/spike_align"
    echo "====== dbitm spike-align dry-run finished ======"
    exit 0
fi

mkdir -p "$final_dir"
if [[ -z "$selected_chunk" ]]; then
    rm -rf -- "$run_output"
fi
mkdir -p "$run_output/logs"
if [[ -n "$selected_chunk" ]]; then
    align_log=$run_output/spike_align.$selected_chunk.log
else
    align_log=$run_output/spike_align.log
fi
: > "$align_log"

declare -a chunks=()
declare -a r1_files=()
declare -a r2_files=()
declare -a chunk_logs=()
if [[ "$assay" == smc ]]; then
    declare -a crick_chunks=()
    declare -a crick_r1_files=()
    declare -a crick_r2_files=()
    declare -a crick_chunk_logs=()
    if [[ -n "$selected_chunk" ]]; then
        watson_short_files=("$run_input/$selected_chunk.watson.short-genomic.fastq.gz")
    else
        shopt -s nullglob
        watson_short_files=("$run_input"/*.watson.short-genomic.fastq.gz)
        shopt -u nullglob
    fi
    if (( ${#watson_short_files[@]} == 0 )); then
        echo "[dbitm] spike-align: no Watson FASTQ chunks found: $run_input/*.watson.short-genomic.fastq.gz" >&2
        exit 1
    fi
    chunk_count=${#watson_short_files[@]}
    for watson_short in "${watson_short_files[@]}"; do
        if [[ ! -f "$watson_short" ]]; then
            echo "[dbitm] spike-align: Watson FASTQ chunk not found: $watson_short" >&2
            exit 1
        fi
        filename=$(basename "$watson_short")
        chunk=${filename%.watson.short-genomic.fastq.gz}
        watson_long=$run_input/$chunk.watson.genomic.fastq.gz
        crick_short=$run_input/$chunk.crick.short-genomic.fastq.gz
        crick_long=$run_input/$chunk.crick.genomic.fastq.gz
        for paired_fastq in "$watson_long" "$crick_short" "$crick_long"; do
            if [[ ! -f "$paired_fastq" ]]; then
                echo "[dbitm] spike-align: paired SmC FASTQ not found for chunk '$chunk': $paired_fastq" >&2
                exit 1
            fi
        done
        chunks+=("$chunk.watson")
        r1_files+=("$watson_long")
        r2_files+=("$watson_short")
        chunk_logs+=("$run_output/logs/$chunk.watson.log")

        crick_chunks+=("$chunk.crick")
        crick_r1_files+=("$crick_short")
        crick_r2_files+=("$crick_long")
        crick_chunk_logs+=("$run_output/logs/$chunk.crick.log")
    done
    chunks+=("${crick_chunks[@]}")
    r1_files+=("${crick_r1_files[@]}")
    r2_files+=("${crick_r2_files[@]}")
    chunk_logs+=("${crick_chunk_logs[@]}")
else
    if [[ -n "$selected_chunk" ]]; then
        r1_files=("$run_input/$selected_chunk.R1.spike-in.fastq.gz")
    else
        shopt -s nullglob
        r1_files=("$run_input"/*.R1.spike-in.fastq.gz)
        shopt -u nullglob
    fi
    if (( ${#r1_files[@]} == 0 )); then
        echo "[dbitm] spike-align: no R1 spike-in FASTQ chunks found: $run_input/*.R1.spike-in.fastq.gz" >&2
        exit 1
    fi
    chunk_count=${#r1_files[@]}
    for r1 in "${r1_files[@]}"; do
        if [[ ! -f "$r1" ]]; then
            echo "[dbitm] spike-align: R1 spike-in FASTQ chunk not found: $r1" >&2
            exit 1
        fi
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
fi

align_threads_per_job=$SPIKE_ALIGN_THREADS_PER_CHUNK

echo "[dbitm] chunks: $chunk_count"
echo "[dbitm] alignments: ${#chunks[@]}"
echo "[dbitm] threads per alignment: $align_threads_per_job"

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
                biscuit align -@ "$align_threads_per_job" -b "$biscuit_directional_mode" \
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

for chunk_index in "${!chunks[@]}"; do
    chunk_name=${chunks[$chunk_index]}
    trap 'job_status=$?; if (( job_status != 0 )); then echo "[dbitm] spike-align: chunk failed: $chunk_name" >&2; fi' EXIT
    align_spike_chunk \
        "$chunk_name" \
        "${r1_files[$chunk_index]}" \
        "${r2_files[$chunk_index]}" \
        "${chunk_logs[$chunk_index]}"
    trap - EXIT
done

for chunk_log in "${chunk_logs[@]}"; do
    if [[ -f "$chunk_log" ]]; then
        printf '\n===== %s =====\n' "$(basename "$chunk_log")" >> "$align_log"
        cat "$chunk_log" >> "$align_log"
    fi
done

echo "[dbitm] spike-align log: $align_log"
echo "[dbitm] spike-align result: $final_dir/spike_align"
echo "====== dbitm spike-align finished ======"
