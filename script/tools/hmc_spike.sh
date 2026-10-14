#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 || $# > 3 )); then
    echo "Usage: hmc_spike.sh <taps|taps-v2> <barcode_dir> [--dry-run]" >&2
    exit 1
fi

assay=$1
barcode_dir=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] hmc-spike: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
case "$assay" in
    taps|taps-v2) ;;
    *) echo "[dbitm] hmc-spike: assay must be taps or taps-v2" >&2; exit 1 ;;
esac
if [[ ! -d "$barcode_dir" ]]; then
    echo "[dbitm] hmc-spike: barcode directory not found: $barcode_dir" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
motif_script=$SCRIPT_DIR/python/find_r1_5hmc_motif.py
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] hmc-spike: config file not found: $config_file" >&2
    exit 1
fi
if [[ ! -f "$motif_script" ]]; then
    echo "[dbitm] hmc-spike: Python script not found: $motif_script" >&2
    exit 1
fi

source "$config_file"

left_flank=$HMC_SPIKE_LEFT_FLANK
middle=$HMC_SPIKE_MIDDLE
right_flank=$HMC_SPIKE_RIGHT_FLANK
c_positions=$HMC_SPIKE_C_POSITIONS
max_jobs=$HMC_SPIKE_THREADS
left_flank=${left_flank^^}
middle=${middle^^}
right_flank=${right_flank^^}
if [[ ! "$left_flank" =~ ^[ACGT]+$ || ! "$middle" =~ ^[ACGTN]+$ || ! "$right_flank" =~ ^[ACGT]+$ ]]; then
    echo "[dbitm] hmc-spike: invalid motif sequence in configuration" >&2
    exit 1
fi
if [[ ! "$c_positions" =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ]]; then
    echo "[dbitm] hmc-spike: HMC_SPIKE_C_POSITIONS must be comma-separated 1-based positions" >&2
    exit 1
fi
IFS=, read -r -a positions <<< "$c_positions"
declare -A seen_positions=()
for position in "${positions[@]}"; do
    if (( position > ${#middle} )) || [[ ${middle:position-1:1} != C || -n ${seen_positions[$position]:-} ]]; then
        echo "[dbitm] hmc-spike: each configured C position must be unique and point to C in the middle motif" >&2
        exit 1
    fi
    seen_positions[$position]=1
done
if [[ ! "$max_jobs" =~ ^[1-9][0-9]*$ ]]; then
    echo "[dbitm] hmc-spike: HMC_SPIKE_THREADS must be a positive integer" >&2
    exit 1
fi

barcode_dir=$(realpath "$barcode_dir")
output_dir=$(dirname "$barcode_dir")/hmc-spike
shopt -s nullglob
candidate_files=("$barcode_dir"/*.R1.spike-in.fastq.gz)
shopt -u nullglob
input_files=()
for candidate_file in "${candidate_files[@]}"; do
    filename=${candidate_file##*/}
    if [[ "$filename" =~ ^[0-9]{4,}\.R1\.spike-in\.fastq\.gz$ ]]; then
        input_files+=("$candidate_file")
    fi
done
if (( ${#input_files[@]} == 0 )); then
    echo "[dbitm] hmc-spike: no numbered R1 spike-in FASTQ files in $barcode_dir" >&2
    exit 1
fi
worker_count=${#input_files[@]}
if (( max_jobs < worker_count )); then
    worker_count=$max_jobs
fi

output_files=()
for input_file in "${input_files[@]}"; do
    if [[ ! -f "$input_file" ]]; then
        echo "[dbitm] hmc-spike: input is not a file: $input_file" >&2
        exit 1
    fi
    stem=${input_file##*/}
    stem=${stem%.fastq.gz}
    output_file=$output_dir/$stem.hmc-spike.txt
    if [[ -e "$output_file" || -L "$output_file" ]]; then
        echo "[dbitm] hmc-spike: output already exists: $output_file" >&2
        exit 1
    fi
    output_files+=("$output_file")
done

echo "[dbitm] hmc-spike: assay=$assay input=$barcode_dir"
echo "[dbitm] hmc-spike: output directory=$output_dir"
echo "[dbitm] hmc-spike: middle C positions (1-based)=$c_positions"
echo "[dbitm] hmc-spike: parallel files=$worker_count (HMC_SPIKE_THREADS=$max_jobs)"

run_motif_file() (
    local output_file=$1
    shift
    local input_file=$1
    local temp_output
    temp_output=$(mktemp "$output_dir/.${input_file##*/}.XXXXXX.tmp")
    trap 'rm -f -- "$temp_output"' EXIT
    if ! pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        python "$motif_script" "$@" > "$temp_output"; then
        echo "[dbitm] hmc-spike: failed while processing $input_file" >&2
        return 1
    fi
    mv -n -- "$temp_output" "$output_file"
    if [[ -e "$temp_output" ]]; then
        echo "[dbitm] hmc-spike: output already exists: $output_file" >&2
        return 1
    fi
)

if [[ "$dry_run" == false ]]; then
    mkdir -p "$output_dir"
fi
active_pids=()
for index in "${!input_files[@]}"; do
    input_file=${input_files[$index]}
    output_file=${output_files[$index]}
    motif_args=(
        "$input_file"
        --left-flank "$left_flank"
        --middle "$middle"
        --right-flank "$right_flank"
        --c-positions "$c_positions"
    )
    if [[ "$dry_run" == true ]]; then
        printf '[dbitm] planned command:'
        printf ' %q' pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
            python "$motif_script" "${motif_args[@]}"
        printf ' > %q\n' "$output_file"
        continue
    fi

    echo "[dbitm] hmc-spike: processing $input_file"
    run_motif_file "$output_file" "${motif_args[@]}" &
    active_pids+=("$!")
    if (( ${#active_pids[@]} >= max_jobs )); then
        failed=false
        for pid in "${active_pids[@]}"; do
            if ! wait "$pid"; then
                failed=true
            fi
        done
        if [[ "$failed" == true ]]; then
            echo "[dbitm] hmc-spike: one or more files failed" >&2
            exit 1
        fi
        active_pids=()
    fi
done
if [[ "$dry_run" == true ]]; then
    echo "[dbitm] hmc-spike: dry-run; no files written"
    exit 0
fi
failed=false
for pid in "${active_pids[@]}"; do
    if ! wait "$pid"; then
        failed=true
    fi
done
if [[ "$failed" == true ]]; then
    echo "[dbitm] hmc-spike: one or more files failed" >&2
    exit 1
fi
for output_file in "${output_files[@]}"; do
    echo "[dbitm] hmc-spike: saved $output_file"
    cat "$output_file"
done
