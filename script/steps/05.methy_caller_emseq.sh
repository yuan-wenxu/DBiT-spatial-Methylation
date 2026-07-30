#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage:
  05.methy_caller_emseq.sh \
    --bam BAM \
    --reference FASTA \
    --out-dir DIR \
    --chromosomes CHR1,CHR2,... \
    [--context-mode cg|ch|both] \
    [--min-base-quality INT] \
    [--min-mapping-quality INT] \
    [--r1-left-trim INT] \
    [--r1-right-trim INT] \
    [--r2-left-trim INT] \
    [--r2-right-trim INT] \
    [--jobs INT] \
    [--max-spots INT]
EOF
}

bam=
reference=
out_dir=
chromosomes=
context_mode=cg
min_base_quality=30
min_mapping_quality=10
r1_left_trim=0
r1_right_trim=0
r2_left_trim=0
r2_right_trim=0
jobs=8
max_spots=10000

while (( $# > 0 )); do
    case "$1" in
        --bam) bam=${2:-}; shift 2 ;;
        --reference) reference=${2:-}; shift 2 ;;
        --out-dir) out_dir=${2:-}; shift 2 ;;
        --chromosomes) chromosomes=${2:-}; shift 2 ;;
        --context-mode) context_mode=${2:-}; shift 2 ;;
        --min-base-quality) min_base_quality=${2:-}; shift 2 ;;
        --min-mapping-quality) min_mapping_quality=${2:-}; shift 2 ;;
        --r1-left-trim) r1_left_trim=${2:-}; shift 2 ;;
        --r1-right-trim) r1_right_trim=${2:-}; shift 2 ;;
        --r2-left-trim) r2_left_trim=${2:-}; shift 2 ;;
        --r2-right-trim) r2_right_trim=${2:-}; shift 2 ;;
        --jobs) jobs=${2:-}; shift 2 ;;
        --max-spots) max_spots=${2:-}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *)
            echo "[methy-caller] unknown argument: $1" >&2
            usage
            exit 2
            ;;
    esac
done

for required_name in bam reference out_dir chromosomes; do
    if [[ -z ${!required_name} ]]; then
        echo "[methy-caller] --${required_name//_/-} is required" >&2
        usage
        exit 2
    fi
done

is_nonnegative_integer() {
    [[ $1 =~ ^[0-9]+$ ]]
}

is_positive_integer() {
    [[ $1 =~ ^[1-9][0-9]*$ ]]
}

if [[ ! -f "$bam" ]]; then
    echo "[methy-caller] pooled BAM not found: $bam" >&2
    exit 1
fi
if [[ ! -f "$reference" ]]; then
    echo "[methy-caller] reference FASTA not found: $reference" >&2
    exit 1
fi
if [[ ! -f "$reference.fai" ]]; then
    echo "[methy-caller] reference FASTA index not found: $reference.fai" >&2
    exit 1
fi
case "$context_mode" in
    cg|ch|both) ;;
    *)
        echo "[methy-caller] --context-mode must be one of: cg, ch, both" >&2
        exit 2
        ;;
esac
if ! is_positive_integer "$jobs"; then
    echo "[methy-caller] --jobs must be a positive integer: $jobs" >&2
    exit 2
fi
if ! is_positive_integer "$max_spots"; then
    echo "[methy-caller] --max-spots must be a positive integer: $max_spots" >&2
    exit 2
fi
for value_name in \
    min_base_quality min_mapping_quality \
    r1_left_trim r1_right_trim r2_left_trim r2_right_trim
do
    if ! is_nonnegative_integer "${!value_name}"; then
        echo "[methy-caller] --${value_name//_/-} must be a non-negative integer: ${!value_name}" >&2
        exit 2
    fi
done
if [[ "$r1_left_trim" != "$r2_left_trim" ]]; then
    echo "[methy-caller] BISCUIT requires equal R1/R2 left trimming" >&2
    exit 2
fi
if [[ "$r1_right_trim" != "$r2_right_trim" ]]; then
    echo "[methy-caller] BISCUIT requires equal R1/R2 right trimming" >&2
    exit 2
fi

for executable in samtools biscuit awk; do
    if ! command -v "$executable" >/dev/null 2>&1; then
        echo "[methy-caller] executable not found: $executable" >&2
        exit 1
    fi
done

mkdir -p "$out_dir"
work_parent=$(dirname "$(realpath -m "$out_dir")")
work_dir=$(mktemp -d "$work_parent/.emseq-call-work.XXXXXX")
split_dir=$work_dir/split
mkdir -p "$split_dir"

cleanup() {
    local status=$?
    trap - EXIT INT TERM HUP
    rm -rf -- "$work_dir"
    exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'exit 129' HUP

echo "[methy-caller] bam=$bam"
echo "[methy-caller] reference=$reference"
echo "[methy-caller] out-dir=$out_dir"
echo "[methy-caller] chromosomes=$chromosomes"
echo "[methy-caller] jobs=$jobs"
echo "[methy-caller] splitting pooled BAM by CB tag"

samtools split \
    -@ "$jobs" \
    -d CB \
    -M "$max_spots" \
    -u "$work_dir/reads-without-cb.bam" \
    -f "$split_dir/%!.bam" \
    "$bam"

mapfile -d '' spot_bams < <(
    find "$split_dir" -maxdepth 1 -type f -name '*.bam' -print0 | sort -z
)
if (( ${#spot_bams[@]} == 0 )); then
    echo "[methy-caller] no reads with a CB tag were found in: $bam" >&2
    exit 1
fi

missing_cb_reads=$(samtools view -c "$work_dir/reads-without-cb.bam")
if (( missing_cb_reads > 0 )); then
    echo "[methy-caller] warning: skipped $missing_cb_reads reads without a CB tag" >&2
fi
echo "[methy-caller] spots=${#spot_bams[@]}"

call_spot() {
    local spot_bam=$1
    local filename cb prefix cg_output ch_output cg_temporary ch_temporary
    local temporary_vcf need_cg=false need_ch=false

    filename=$(basename "$spot_bam")
    cb=${filename%.bam}
    prefix=${cb:0:2}
    cg_output=$out_dir/host/$prefix/$cb.CG.cov
    ch_output=$out_dir/host/$prefix/$cb.CH.cov
    cg_temporary=$cg_output.tmp
    ch_temporary=$ch_output.tmp
    temporary_vcf=$spot_bam.vcf
    mkdir -p "$(dirname "$cg_output")"

    if [[ "$context_mode" == cg || "$context_mode" == both ]]; then
        if [[ -f "$cg_output" ]]; then
            echo "[methy-caller] skip-existing-output=$cg_output"
        else
            need_cg=true
        fi
    fi
    if [[ "$context_mode" == ch || "$context_mode" == both ]]; then
        if [[ -f "$ch_output" ]]; then
            echo "[methy-caller] skip-existing-output=$ch_output"
        else
            need_ch=true
        fi
    fi
    if [[ "$need_cg" == false && "$need_ch" == false ]]; then
        return 0
    fi

    samtools index "$spot_bam"
    if ! biscuit pileup \
        -b "$min_base_quality" \
        -m "$min_mapping_quality" \
        -a 0 \
        -c \
        -u \
        -p \
        -5 "$r1_left_trim" \
        -3 "$r1_right_trim" \
        -@ 1 \
        "$reference" \
        "$spot_bam" \
        > "$temporary_vcf"
    then
        rm -f -- "$temporary_vcf"
        return 1
    fi

    if [[ "$need_cg" == true ]] && ! (
        set -o pipefail
        biscuit vcf2bed -t cg "$temporary_vcf" \
        | biscuit mergecg -c "$reference" - \
        | awk -v chromosome_csv="$chromosomes" '
            BEGIN {
                OFS = "\t"
                count = split(chromosome_csv, values, ",")
                for (i = 1; i <= count; i++) {
                    gsub(/^[[:space:]]+|[[:space:]]+$/, "", values[i])
                    if (values[i] != "") {
                        retained[values[i]] = 1
                    }
                }
            }
            NF >= 6 && ($1 in retained) {
                methylated = $5 + 0
                unmethylated = $6 + 0
                depth = methylated + unmethylated
                if (depth == 0) {
                    next
                }
                position = $2 + 1
                percentage = 100 * methylated / depth
                printf "%s\t%d\t%d\t%.2f\t%d\t%d\n",
                    $1, position, position, percentage, methylated, unmethylated
            }
        ' > "$cg_temporary"
    ); then
        rm -f -- "$temporary_vcf" "$cg_temporary" "$ch_temporary"
        return 1
    fi

    if [[ "$need_ch" == true ]] && ! (
        set -o pipefail
        biscuit vcf2bed -t ch -e -c "$temporary_vcf" \
        | awk -v chromosome_csv="$chromosomes" '
            BEGIN {
                OFS = "\t"
                count = split(chromosome_csv, values, ",")
                for (i = 1; i <= count; i++) {
                    gsub(/^[[:space:]]+|[[:space:]]+$/, "", values[i])
                    if (values[i] != "") {
                        retained[values[i]] = 1
                    }
                }
            }
            NF >= 10 && ($1 in retained) {
                methylated = $9 + 0
                unmethylated = $10 + 0
                depth = methylated + unmethylated
                if (depth == 0) {
                    next
                }
                position = $2 + 1
                percentage = 100 * methylated / depth
                strand = ($4 == "C") ? "+" : "-"
                printf "%s\t%d\t%d\t%.2f\t%d\t%d\t%s\t%s\n",
                    $1, position, position, percentage, methylated,
                    unmethylated, $6, strand
            }
        ' > "$ch_temporary"
    ); then
        rm -f -- "$temporary_vcf" "$cg_temporary" "$ch_temporary"
        return 1
    fi

    rm -f -- "$temporary_vcf"
    if [[ "$need_cg" == true ]]; then
        mv -- "$cg_temporary" "$cg_output"
        echo "[methy-caller] spot=$cb output=$cg_output"
    fi
    if [[ "$need_ch" == true ]]; then
        mv -- "$ch_temporary" "$ch_output"
        echo "[methy-caller] spot=$cb output=$ch_output"
    fi
}

declare -a active_pids=()
declare -a active_spots=()
calling_failed=false

wait_for_active_jobs() {
    local index
    for index in "${!active_pids[@]}"; do
        if ! wait "${active_pids[$index]}"; then
            echo "[methy-caller] spot failed: ${active_spots[$index]}" >&2
            calling_failed=true
        fi
    done
    active_pids=()
    active_spots=()
}

for spot_bam in "${spot_bams[@]}"; do
    call_spot "$spot_bam" &
    active_pids+=("$!")
    active_spots+=("$(basename "$spot_bam" .bam)")
    if (( ${#active_pids[@]} >= jobs )); then
        wait_for_active_jobs
    fi
done
if (( ${#active_pids[@]} > 0 )); then
    wait_for_active_jobs
fi
if [[ "$calling_failed" == true ]]; then
    echo "[methy-caller] one or more spots failed" >&2
    exit 1
fi

echo "[methy-caller] done spots=${#spot_bams[@]}"
