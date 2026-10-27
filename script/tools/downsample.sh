#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 )); then
    echo "Usage: downsample.sh <assay> <raw_fastq_folder> (--fraction <float> | --reads <int>) [--dry-run]" >&2
    exit 1
fi

assay=$1
raw_path=$2
shift 2

fraction=""
reads=""
dry_run=false
while (( $# > 0 )); do
    case "$1" in
        --fraction|--reads)
            if (( $# < 2 )) || [[ -z "$2" ]]; then
                echo "[dbitm] downsample: $1 requires a value" >&2
                exit 1
            fi
            case "$1" in
                --fraction) fraction=$2 ;;
                --reads) reads=$2 ;;
            esac
            shift 2
            ;;
        --dry-run)
            dry_run=true
            shift
            ;;
        *)
            echo "[dbitm] downsample: unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

case "$assay" in
    taps|taps-v2|emseq|cabernet|smc) ;;
    *) echo "[dbitm] downsample: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] downsample: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] downsample: config file not found: $config_file" >&2
    exit 1
fi
# shellcheck source=/dev/null
source "$config_file"

if [[ -n "$fraction" ]]; then
    if [[ -n "$reads" ]]; then
        echo "[dbitm] downsample: choose --fraction or --reads, not both" >&2
        exit 1
    fi
    if [[ ! "$fraction" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]] ||
        ! awk -v value="$fraction" 'BEGIN { exit !(value >= 0 && value <= 1) }'; then
        echo "[dbitm] downsample: --fraction must be in [0, 1]" >&2
        exit 1
    fi
elif [[ ! "$reads" =~ ^(0|[1-9][0-9]*)$ || ${#reads} -gt 18 ]]; then
    echo "[dbitm] downsample: --reads must be a non-negative integer (at most 18 digits)" >&2
    exit 1
fi
seed=${MBIAS_SAMPLING_SEED:-}
if [[ ! "$seed" =~ ^(0|[1-9][0-9]*)$ || ${#seed} -gt 18 ]]; then
    echo "[dbitm] downsample: MBIAS_SAMPLING_SEED must be a non-negative integer (at most 18 digits)" >&2
    exit 1
fi

raw_abs=$(realpath "$raw_path")
output_dir=$raw_abs/downsample
if [[ -e "$output_dir" && ! -d "$output_dir" ]]; then
    echo "[dbitm] downsample: output is not a directory: $output_dir" >&2
    exit 1
fi

fail() {
    echo "[dbitm] downsample: $*" >&2
    exit 1
}

declare -a r1_files=()
declare -a r2_files=()
shopt -s nullglob nocaseglob
declare -a candidates=(
    "$raw_abs"/*.fastq "$raw_abs"/*.fq
    "$raw_abs"/*.fastq.gz "$raw_abs"/*.fq.gz
)
shopt -u nullglob nocaseglob
declare -A first_mates=() second_mates=()
for path in "${candidates[@]}"; do
    [[ -f "$path" ]] || continue
    filename=${path##*/}
    if [[ ! "$filename" =~ ^((.*[^[:alnum:]])?)[Rr]([12])([^[:alnum:]].*|$)$ ]]; then
        fail "filename needs an independent R1/R2 token: $filename"
    fi
    key=${BASH_REMATCH[1]}R#${BASH_REMATCH[4]}
    mate=${BASH_REMATCH[3]}
    if [[ "$mate" == 1 ]]; then
        [[ -z ${first_mates[$key]:-} ]] || fail "duplicate R1: $filename"
        first_mates[$key]=$path
    else
        [[ -z ${second_mates[$key]:-} ]] || fail "duplicate R2: $filename"
        second_mates[$key]=$path
    fi
done
(( ${#first_mates[@]} > 0 )) || fail "no R1 FASTQ files found"
(( ${#first_mates[@]} == ${#second_mates[@]} )) || fail "missing mates"
# Sort the complete path list rather than associative-array keys.
while IFS= read -r -d '' path; do
    filename=${path##*/}
    [[ "$filename" =~ ^((.*[^[:alnum:]])?)[Rr]([12])([^[:alnum:]].*|$)$ ]]
    key=${BASH_REMATCH[1]}R#${BASH_REMATCH[4]}
    [[ -n ${second_mates[$key]:-} ]] || fail "missing R2: $filename"
    r1_files+=("$path")
    r2_files+=("${second_mates[$key]}")
done < <(printf '%s\0' "${first_mates[@]}" | sort -z)

declare -a input_files=("${r1_files[@]}" "${r2_files[@]}")
for path in "${input_files[@]}"; do
    target=$output_dir/${path##*/}
    [[ ! -e "$target" && ! -L "$target" ]] ||
        fail "output already exists: $target"
done

echo "====== dbitm downsample ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] input directory: $raw_abs"
echo "[dbitm] output directory: $output_dir"
echo "[dbitm] config: $config_file"
echo "[dbitm] input pairs: ${#r1_files[@]}"
echo "[dbitm] sampling seed: $seed"
if [[ -n "$fraction" ]]; then
    echo "[dbitm] sampling fraction: $fraction"
else
    echo "[dbitm] target read pairs per FASTQ pair: $reads"
fi

if [[ "$dry_run" == true ]]; then
    for index in "${!r1_files[@]}"; do
        echo "[dbitm] input pair $((index + 1)) R1: ${r1_files[$index]}"
        echo "[dbitm] input pair $((index + 1)) R2: ${r2_files[$index]}"
    done
    echo "[dbitm] planned command: seqtk sample -> $output_dir"
    echo "[dbitm] dry-run: records not scanned; no files will be written"
    echo "====== dbitm downsample dry-run finished ======"
    exit 0
fi

mkdir -p "$output_dir"
temp_dir=$(mktemp -d "$output_dir/.downsample.XXXXXX")
declare -a published=()
finished=false
cleanup() {
    if [[ "$finished" == false && ${#published[@]} -gt 0 ]]; then
        rm -f -- "${published[@]}"
    fi
    rm -rf -- "$temp_dir"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

declare -a seqtk_args=(sample "-s$seed")
if [[ -n "$fraction" ]]; then
    target_count=$fraction
    if awk -v value="$fraction" 'BEGIN { exit !(value == 1) }'; then
        # seqtk sample 1 means one record rather than all records.
        seqtk_args=(seq)
    fi
else
    target_count=$reads
    (( reads == 0 )) || seqtk_args+=(-2)
fi

for path in "${input_files[@]}"; do
    sample_command=(
        pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default
        seqtk "${seqtk_args[@]}" "$path"
    )
    [[ ${seqtk_args[0]} == seq ]] || sample_command+=("$target_count")
    printf '[dbitm] command:'
    printf ' %q' "${sample_command[@]}"
    printf '\n'
    target=$temp_dir/${path##*/}
    if [[ "${path,,}" == *.gz ]]; then
        "${sample_command[@]}" | \
            pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
                pigz -p 1 -c > "$target"
    else
        "${sample_command[@]}" > "$target"
    fi
done

for path in "${input_files[@]}"; do
    target=$output_dir/${path##*/}
    ln -- "$temp_dir/${path##*/}" "$target"
    published+=("$target")
done
finished=true
echo "[dbitm] downsample result: $output_dir"
echo "====== dbitm downsample finished ======"
