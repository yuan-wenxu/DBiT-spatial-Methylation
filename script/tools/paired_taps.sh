#!/usr/bin/env bash
set -euo pipefail

if (( $# < 1 )); then
    echo "Usage: paired_taps.sh <assay> --taps-cov <file> --taps-beta-cov <file> --spot-map <file> [--dry-run]" >&2
    exit 1
fi

assay=$1
shift
taps_cov=""
taps_beta_cov=""
spot_map=""
dry_run=false
while (( $# > 0 )); do
    case "$1" in
        --taps-cov)
            taps_cov=$2
            shift 2
            ;;
        --taps-beta-cov)
            taps_beta_cov=$2
            shift 2
            ;;
        --spot-map)
            spot_map=$2
            shift 2
            ;;
        --dry-run)
            dry_run=true
            shift
            ;;
        *)
            echo "[dbitm] paired-taps: unknown argument: $1" >&2
            exit 1
            ;;
    esac
done
case "$assay" in
    taps|taps-v2) ;;
    *)
        echo "[dbitm] paired-taps: assay must be taps or taps-v2" >&2
        exit 1
        ;;
esac
for required_option in taps_cov taps_beta_cov spot_map; do
    if [[ -z ${!required_option} ]]; then
        echo "[dbitm] paired-taps: --${required_option//_/-} is required" >&2
        exit 1
    fi
done

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=${DBITM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}
config_file=${DBITM_CONFIG:-$REPO_DIR/config/dbitm.config.sh}
paired_taps_script=$SCRIPT_DIR/python/paired_taps.py
if [[ ! -f "$config_file" ]]; then
    echo "[dbitm] paired-taps: config file not found: $config_file" >&2
    exit 1
fi
if [[ ! -f "$paired_taps_script" ]]; then
    echo "[dbitm] paired-taps: Python script not found: $paired_taps_script" >&2
    exit 1
fi
# shellcheck source=/dev/null
source "$config_file"

for input_path in "$taps_cov" "$taps_beta_cov" "$spot_map"; do
    if [[ ! -f "$input_path" ]]; then
        echo "[dbitm] paired-taps: input file not found: $input_path" >&2
        exit 1
    fi
done
taps_cov=$(realpath "$taps_cov")
taps_beta_cov=$(realpath "$taps_beta_cov")
spot_map=$(realpath "$spot_map")
output_dir=$(dirname "$spot_map")/paired-taps

threads=${PAIRED_TAPS_THREADS:-8}
declare -a paired_taps_args=(
    --taps-cov "$taps_cov"
    --taps-beta-cov "$taps_beta_cov"
    --spot-map "$spot_map"
    --bin-size "${PAIRED_TAPS_BIN_SIZE}"
    --min-coverage-per-site "${PAIRED_TAPS_MIN_COVERAGE_PER_SITE}"
    --min-coverage-per-bin "${PAIRED_TAPS_MIN_COVERAGE_PER_BIN}"
    --min-site-per-bin "${PAIRED_TAPS_MIN_SITE_PER_BIN}"
    --min-valid-spots-per-bin "${PAIRED_TAPS_MIN_VALID_SPOTS_PER_BIN}"
    --variable-fraction "${PAIRED_TAPS_VARIABLE_FRACTION}"
    --shrinkage-lambda "${PAIRED_TAPS_SHRINKAGE_LAMBDA}"
    --chunksize "${PAIRED_TAPS_CHUNKSIZE}"
    --threads "$threads"
)
if [[ -n ${PAIRED_TAPS_MAX_COVERAGE_PER_SITE:-} ]]; then
    paired_taps_args+=(
        --max-coverage-per-site "$PAIRED_TAPS_MAX_COVERAGE_PER_SITE"
    )
fi

echo "====== dbitm paired TAPS/TAPS-beta tool ======"
echo "[dbitm] TAPS coverage: $taps_cov"
echo "[dbitm] TAPS-beta coverage: $taps_beta_cov"
echo "[dbitm] spot map: $spot_map"
echo "[dbitm] output directory: $output_dir"
echo "[dbitm] worker processes: $threads"

if [[ "$dry_run" == true ]]; then
    printf '[dbitm] planned command:'
    printf ' %q' pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
        python "$paired_taps_script" "${paired_taps_args[@]}"
    printf '\n'
    echo "[dbitm] dry-run: no files will be written"
    echo "====== dbitm paired TAPS/TAPS-beta dry-run finished ======"
    exit 0
fi

pixi run --manifest-path "$REPO_DIR/pixi.toml" -e default \
    python "$paired_taps_script" "${paired_taps_args[@]}"
echo "[dbitm] paired TAPS/TAPS-beta output: $output_dir"
echo "====== dbitm paired TAPS/TAPS-beta finished ======"
