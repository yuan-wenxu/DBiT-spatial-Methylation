#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 || $# > 3 )); then
    echo "Usage: 03.align_prepare.sh <assay> <raw_fastq_folder> [--dry-run]" >&2
    exit 1
fi
assay=$1
raw_path=$2
dry_run=false
if (( $# == 3 )); then
    if [[ $3 != --dry-run ]]; then
        echo "[dbitm] align-prepare: unknown argument: $3" >&2
        exit 1
    fi
    dry_run=true
fi
case "$assay" in
    taps|taps-v2|emseq|cabernet|smc) ;;
    *) echo "[dbitm] align-prepare: unsupported assay: $assay" >&2; exit 1 ;;
esac
if [[ ! -d "$raw_path" ]]; then
    echo "[dbitm] align-prepare: FASTQ directory not found: $raw_path" >&2
    exit 1
fi

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] align-prepare: config file not found: $config_file" >&2
    exit 1
fi

raw_abs=$(realpath "$raw_path")
final_dir=$(dirname "$raw_abs")/dbitm
run_output=$final_dir/align

echo "====== dbitm align prepare ======"
echo "[dbitm] assay: $assay"
echo "[dbitm] config: $config_file"
echo "[dbitm] output directory: $run_output"
if [[ "$dry_run" == true ]]; then
    echo "[dbitm] dry-run: no files will be written"
    echo "[dbitm] planned action: recreate $run_output"
    echo "====== dbitm align prepare dry-run finished ======"
    exit 0
fi

mkdir -p "$final_dir"
rm -rf -- "$run_output"
mkdir -p "$run_output/logs"
: > "$run_output/align.log"

echo "[dbitm] align output prepared: $run_output"
echo "====== dbitm align prepare finished ======"
