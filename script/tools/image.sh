#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 || $# > 3 )); then
    echo "Usage: image.sh <assay> <full_resolution_image> [--dry-run]" >&2
    exit 1
fi
assay=$1
image_path=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] image: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
case "$assay" in
    taps|taps-v2|emseq|cabernet|smc) ;;
    *) echo "[dbitm] image: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -f "$image_path" ]]; then
    echo "[dbitm] image: source image not found: $image_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
segment_script=$SCRIPT_DIR/python/image_segment.py
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] image: config file not found: $config_file" >&2
    exit 1
fi
if [[ ! -f "$segment_script" ]]; then
    echo "[dbitm] image: segmentation script not found: $segment_script" >&2
    exit 1
fi
# shellcheck source=/dev/null
source "$config_file"

image_path=$(realpath "$image_path")
image_dir=$(dirname "$image_path")
mask_path=$image_dir/meth-mask.png
output_dir=$image_dir/meth-image
methscan_dir=$(realpath -m "$(dirname "$image_dir")/meth/dbitm/methscan")
if [[ ! -f "$mask_path" ]]; then
    echo "[dbitm] image: adjacent registered mask not found: $mask_path" >&2
    exit 1
fi

if [[ "$assay" == smc ]]; then
    barcode_whitelist=${SMC_BARCODE_WHITELIST:-$REPO_DIR/docs/barcodes/barcodes-smc.tsv}
else
    barcode_whitelist=${BARCODE_WHITELIST:-$REPO_DIR/docs/barcodes/barcodes50.tsv}
fi
if [[ "$barcode_whitelist" != /* ]]; then
    barcode_whitelist=$REPO_DIR/$barcode_whitelist
fi
if [[ ! -f "$barcode_whitelist" ]]; then
    echo "[dbitm] image: barcode whitelist not found: $barcode_whitelist" >&2
    exit 1
fi
barcode_whitelist=$(realpath "$barcode_whitelist")
spot_count=$(awk 'NF && $1 !~ /^#/ { count++ } END { print count + 0 }' "$barcode_whitelist")
if [[ ! "$spot_count" =~ ^[1-9][0-9]*$ ]]; then
    echo "[dbitm] image: barcode whitelist has no usable entries" >&2
    exit 1
fi

spot_length=${SUMMARY_FRAME_SPOT_LENGTH:-50}
interval=${SUMMARY_FRAME_INTERVAL:-50}
pixel_length=${SUMMARY_FRAME_PIXEL_LENGTH:-0.294}
if [[ ! "$spot_length" =~ ^[1-9][0-9]*$ ]]; then
    echo "[dbitm] image: SUMMARY_FRAME_SPOT_LENGTH must be > 0" >&2
    exit 1
fi
if [[ ! "$interval" =~ ^[0-9]+$ ]]; then
    echo "[dbitm] image: SUMMARY_FRAME_INTERVAL must be >= 0" >&2
    exit 1
fi
if [[ ! "$pixel_length" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]] || \
    ! awk -v value="$pixel_length" 'BEGIN { exit !(value > 0) }'; then
    echo "[dbitm] image: SUMMARY_FRAME_PIXEL_LENGTH must be > 0" >&2
    exit 1
fi

declare -a segment_args=(
    --image_path "$image_path"
    --mask_path "$mask_path"
    --result_path "$output_dir"
    --barcodeA_whitelist "$barcode_whitelist"
    --barcodeB_whitelist "$barcode_whitelist"
    --x_spots_number "$spot_count"
    --y_spots_number "$spot_count"
    --length_spot "$spot_length"
    --interval "$interval"
    --pixel_length "$pixel_length"
)

declare -a matrix_dirs=()
if [[ -d "$methscan_dir" ]]; then
    while IFS= read -r -d '' matrix_dir; do
        if [[ -f "$matrix_dir/matrix.mtx.gz" && \
              -f "$matrix_dir/barcodes.tsv.gz" && \
              -f "$matrix_dir/features.tsv.gz" ]]; then
            matrix_dirs+=("$matrix_dir")
        else
            echo "[dbitm] image: skipping incomplete 10x directory: $matrix_dir" >&2
        fi
    done < <(
        find "$methscan_dir" -type d \
            \( -path '*/matrix/10x/methylation_fractions' -o \
               -path '*/matrix/10x/mean_shrunken_residuals' \) \
            -print0 | sort -z
    )
fi

echo "====== dbitm image segmentation tool ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] source image: $image_path"
echo "[dbitm] frame mask: $mask_path"
echo "[dbitm] output directory: $output_dir"
echo "[dbitm] MethSCAn directory: $methscan_dir"
echo "[dbitm] 10x copy targets: ${#matrix_dirs[@]}"
echo "[dbitm] barcode whitelist: $barcode_whitelist"
echo "[dbitm] chip size: ${spot_count}x${spot_count}"
echo "[dbitm] spot length: $spot_length"
echo "[dbitm] interval: $interval"
echo "[dbitm] pixel length: $pixel_length"

if [[ "$dry_run" == true ]]; then
    printf '[dbitm] planned command:'
    printf ' %q' pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        python "$segment_script" "${segment_args[@]}"
    printf '\n'
    for matrix_dir in "${matrix_dirs[@]}"; do
        echo "[dbitm] planned copy: $output_dir/meth-tissue_positions.tsv.gz -> $matrix_dir/tissue_positions.tsv.gz"
        echo "[dbitm] planned copy: $output_dir/meth-fullres_grayscale.png -> $matrix_dir/tissue_raw_image.png"
    done
    echo "[dbitm] dry-run: no files will be written"
    echo "====== dbitm image segmentation dry-run finished ======"
    exit 0
fi

mkdir -p "$output_dir"
pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
    python "$segment_script" "${segment_args[@]}" \
    2>&1 | tee "$output_dir/meth-seg.log"
positions_path=$output_dir/meth-tissue_positions.tsv.gz
grayscale_path=$output_dir/meth-fullres_grayscale.png
if [[ ! -f "$positions_path" || ! -f "$grayscale_path" ]]; then
    echo "[dbitm] image: required segmentation outputs are missing" >&2
    exit 1
fi
for matrix_dir in "${matrix_dirs[@]}"; do
    cp -f -- "$positions_path" "$matrix_dir/tissue_positions.tsv.gz"
    cp -f -- "$grayscale_path" "$matrix_dir/tissue_raw_image.png"
    echo "[dbitm] copied image outputs: $matrix_dir"
done
if (( ${#matrix_dirs[@]} == 0 )); then
    echo "[dbitm] image: warning: no complete MethSCAn 10x directories found; no outputs copied" >&2
fi
echo "[dbitm] tissue positions: $output_dir/meth-tissue_positions.tsv.gz"
echo "====== dbitm image segmentation finished ======"
