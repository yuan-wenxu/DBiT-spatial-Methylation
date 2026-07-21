#!/usr/bin/env bash
set -euo pipefail

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ ${use_scratch:-false} == true && -n ${scratch_run:-} && -d $scratch_run ]]; then
        rm -rf -- "$scratch_run" || echo "[dbitm] barcode: warning: failed to clean scratch directory: $scratch_run" >&2
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
    echo "Usage: 02.barcode.sh <assay> <raw_fastq_folder>" >&2
    exit 1
fi
assay=$1
raw_path=$2
case "$assay" in
    taps|taps-v2|emseq) ;;
    *) echo "[dbitm] barcode: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] barcode: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file="$REPO_DIR/config/dbitm.config.sh"
python_script="$REPO_DIR/script/steps/python/02.extract_bc.py"
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] barcode: config file not found: $config_file" >&2
    exit 1
fi
if [[ ! -f "$python_script" ]]; then
    echo "[dbitm] barcode: Python script not found: $python_script" >&2
    exit 1
fi
source "$config_file"

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
fastp_dir=$final_dir/fastp
input_r1=$fastp_dir/R1.filtered.fastq.gz
input_r2=$fastp_dir/R2.filtered.fastq.gz
if [[ ! -f "$input_r1" || ! -f "$input_r2" ]]; then
    echo "[dbitm] barcode: fastp outputs not found under: $fastp_dir" >&2
    echo "[dbitm] barcode: expected R1.filtered.fastq.gz and R2.filtered.fastq.gz" >&2
    exit 1
fi

barcode_whitelist=${BARCODE_WHITELIST:-$REPO_DIR/docs/barcodes/barcodes50.tsv}
if [[ "$barcode_whitelist" != /* ]]; then
    barcode_whitelist=$REPO_DIR/$barcode_whitelist
fi
if [[ ! -f "$barcode_whitelist" ]]; then
    echo "[dbitm] barcode: whitelist not found: $barcode_whitelist" >&2
    exit 1
fi
barcode_whitelist=$(realpath "$barcode_whitelist")

if [[ ! "$BARCODE_CHUNK" =~ ^[1-9][0-9]*$ ]] ; then
    echo "[dbitm] barcode: BARCODE_CHUNK must be >= 1" >&2
    exit 1
fi
if [[ ! "$BARCODE_BATCH_SIZE" =~ ^[1-9][0-9]*$ ]] ; then
    echo "[dbitm] barcode: BARCODE_BATCH_SIZE must be >= 1" >&2
    exit 1
fi
case "$BARCODE_COMPRESSION_STEP" in python|shell) ;;
    *) echo "[dbitm] barcode: BARCODE_COMPRESSION_STEP must be python or shell" >&2; exit 1 ;;
esac
if [[ ! "$BARCODE_GZIP_LEVEL" =~ ^[0-9]$ ]] ; then
    echo "[dbitm] barcode: BARCODE_GZIP_LEVEL must be between 0 and 9" >&2
    exit 1
fi
if [[ ! "$BARCODE_LINKER_EDIT_DISTANCE" =~ ^[0-9]+$ ]] ; then
    echo "[dbitm] barcode: BARCODE_LINKER_EDIT_DISTANCE must be >= 0" >&2
    exit 1
fi
if [[ ! "$BARCODE_HAMMING_DISTANCE" =~ ^[0-9]+$ ]] ; then
    echo "[dbitm] barcode: BARCODE_HAMMING_DISTANCE must be >= 0" >&2
    exit 1
fi
if [[ ! "$BARCODE_INSERT_LEFT_EDIT_DISTANCE" =~ ^[0-9]+$ ]] ; then
    echo "[dbitm] barcode: BARCODE_INSERT_LEFT_EDIT_DISTANCE must be >= 0" >&2
    exit 1
fi
if [[ ! "$BARCODE_PROGRESS_READS" =~ ^[0-9]+$ ]] ; then
    echo "[dbitm] barcode: BARCODE_PROGRESS_READS must be >= 0" >&2
    exit 1
fi

use_scratch=false
run_r1=$input_r1
run_r2=$input_r2
run_output=$final_dir/barcode

echo "====== dbitm barcode ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] input R1: $input_r1"
echo "[dbitm] input R2: $input_r2"
echo "[dbitm] config: $config_file"
echo "[dbitm] whitelist: $barcode_whitelist"
echo "[dbitm] chunks: $BARCODE_CHUNK"
echo "[dbitm] batch size: $BARCODE_BATCH_SIZE read pairs"
echo "[dbitm] progress interval: $BARCODE_PROGRESS_READS read pairs per chunk"
if [[ "$assay" == taps-v2 ]]; then
    echo "[dbitm] methylated C positions: $BARCODE_METHYLATED_C_POSITIONS"
fi
echo "[dbitm] compression step: $BARCODE_COMPRESSION_STEP"
echo "[dbitm] output directory: $final_dir/barcode"

if [[ -n ${SCRATCH_ROOT:-} ]]; then
    if [[ "$SCRATCH_ROOT" != /* ]]; then
        echo "[dbitm] barcode: SCRATCH_ROOT must be an absolute path or empty" >&2
        exit 1
    fi
    echo "[dbitm] scratch root: $SCRATCH_ROOT"
    mkdir -p "$SCRATCH_ROOT"
    scratch_root=$(realpath "$SCRATCH_ROOT")
    run_id=${SLURM_JOB_ID:-barcode_$$}
    scratch_run=$scratch_root/dbitm/$run_id
    scratch_input=$scratch_run/input
    run_output=$scratch_run/barcode
    run_r1=$scratch_input/R1.filtered.fastq.gz
    run_r2=$scratch_input/R2.filtered.fastq.gz
    use_scratch=true
    enable_cleanup
    mkdir -p "$scratch_input"
    echo "[dbitm] copying fastp FASTQ files to scratch: $scratch_input"
    cp -a "$input_r1" "$run_r1"
    cp -a "$input_r2" "$run_r2"
else
    mkdir -p "$final_dir"
fi

rm -rf -- "$run_output"
mkdir -p "$run_output"
echo "[dbitm] starting barcode extraction..."

python_compression=gzip
if [[ "$BARCODE_COMPRESSION_STEP" == shell ]]; then
    python_compression=none
fi

pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default python "$python_script" \
    --assay "$assay" \
    "$run_r1" "$run_r2" \
    -b1 "$barcode_whitelist" \
    -b2 "$barcode_whitelist" \
    -o "$run_output" \
    --chunks "$BARCODE_CHUNK" \
    --batch-size "$BARCODE_BATCH_SIZE" \
    --compression "$python_compression" \
    --linker-bc "$BARCODE_LINKER_BC" \
    --insert-left "$BARCODE_INSERT_LEFT" \
    --methylated-c-positions "$BARCODE_METHYLATED_C_POSITIONS" \
    --linker-edit-distance "$BARCODE_LINKER_EDIT_DISTANCE" \
    --barcode-hamming-distance "$BARCODE_HAMMING_DISTANCE" \
    --insert-left-edit-distance "$BARCODE_INSERT_LEFT_EDIT_DISTANCE" \
    --progress-reads "$BARCODE_PROGRESS_READS" \
    --gzip-level "$BARCODE_GZIP_LEVEL" \
    > "$run_output/barcode.log" 2>&1

echo "[dbitm] barcode extraction finished successfully"
if [[ "$BARCODE_COMPRESSION_STEP" == shell ]]; then
    fastq_files=("$run_output"/*.fastq)
    if [[ -e "${fastq_files[0]}" ]]; then
        echo "[dbitm] compressing barcode FASTQ files with pigz..."
        compression_pids=()
        for ((chunk_index = 1; chunk_index <= BARCODE_CHUNK; chunk_index++)); do
            chunk_prefix=$(printf '%04d' "$chunk_index")
            chunk_fastq_files=("$run_output/$chunk_prefix".*.fastq)
            pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
                pigz -p 1 -"$BARCODE_GZIP_LEVEL" \
                -- "${chunk_fastq_files[@]}" \
                >> "$run_output/barcode.log" 2>&1 &
            compression_pids+=("$!")
        done
        compression_failed=false
        for compression_pid in "${compression_pids[@]}"; do
            if ! wait "$compression_pid"; then
                compression_failed=true
            fi
        done
        if [[ "$compression_failed" == true ]]; then
            echo "[dbitm] barcode: pigz compression failed; see barcode.log" >&2
            exit 1
        fi
        sed -i \
            -e 's/\.fastq"/\.fastq.gz"/g' \
            -e 's/"compression": "none"/"compression": "gzip"/g' \
            "$run_output"/*.json
        echo "[dbitm] pigz compression finished successfully"
    else
        echo "[dbitm] no barcode FASTQ files to compress"
    fi
fi
if [[ "$use_scratch" == true ]]; then
    rm -rf -- "$final_dir/barcode"
    mkdir -p "$final_dir"
    echo "[dbitm] copying barcode result from scratch to $final_dir/barcode"
    cp -a "$run_output" "$final_dir/barcode"
fi
echo "[dbitm] barcode log: $final_dir/barcode/barcode.log"
echo "[dbitm] barcode result: $final_dir/barcode"
echo "====== dbitm barcode finished ======"
