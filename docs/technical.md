# DBiT-spatial-Methylation Technical Reference

## Overview

DBiT-spatial-Methylation processes paired-end spatial methylation FASTQs for
`taps`, `taps-v2`, `emseq`, `cabernet`, and `smc`. It produces aligned BAMs,
per-spot methylation coverage, QC summaries, and optional MethSCAn VMR matrices.

```text
fastp -> barcode -> (spike-align + align) -> pool -> mbias -> [smc-filter]
    -> call -> saturation -> summary -> methscan
    -> spike-call ----------> summary
```

## Installation and usage

Run from the repository root in Linux with Pixi installed:

```bash
pixi run init
cp config/dbitm.config.example.sh config/dbitm.config.sh
```

The installer adds `dbitm` to `~/.local/bin`; reload the shell if needed.
`DBITM_INSTALL_DIR` changes the installation directory and `DBITM_BASHRC`
changes the shell startup file updated by the installer.

Place one paired R1/R2 FASTQ set in the input directory. Set the references
and execution mode in `config/dbitm.config.sh` before running. Build indexes
with the aligner used by the assay and index each calling FASTA:

```bash
# TAPS / TAPS-v2
pixi run -e default bwa index /path/to/genome.fa
# EM-seq / Cabernet / SmC
pixi run -e default biscuit index /path/to/genome.fa
```

Prepare spike-in references the same way when using them.

```bash
dbitm taps all --input /data/sample/fastq --dry-run
dbitm taps all --input /data/sample/fastq
```

| Argument | Meaning |
|---|---|
| `<assay>` | `taps`, `taps-v2`, `emseq`, `cabernet`, or `smc`. |
| `<step>` | One stage shown in the overview, `all`, `image`, or `paired-taps`. |
| `--input PATH` | Raw FASTQ directory; for `image`, the full-resolution image file. |
| `--config PATH` | Configuration file; default: `config/dbitm.config.sh`. |
| `--resume STEP` | With `all`, start at this stage and run all later stages. |
| `--dry-run` | Validate inputs/configuration and print the plan without writing results. |
| `-h`, `--help` | Show command help. |
| `--taps-cov FILE` | TAPS coverage input for `paired-taps`. |
| `--taps-beta-cov FILE` | TAPS-beta coverage input for `paired-taps`. |
| `--spot-map FILE` | Spot-pairing table for `paired-taps`. |

Standalone examples:

```bash
# Requires mask.png beside the image; results go into image-segmentation/.
dbitm taps image --input /data/sample/image/tissue.png

# Results go into paired-taps/ beside the spot-pairing table.
dbitm taps paired-taps \
    --taps-cov /path/to/taps/host.CG.cov \
    --taps-beta-cov /path/to/taps-beta/host.CG.cov \
    --spot-map /path/to/spot-pairs.tsv
```

`paired-taps` accepts `taps` and `taps-v2`. `image` uses the barcode whitelist
and `SUMMARY_FRAME_*` geometry settings.

## Library structure

The barcode-bearing R1 in TAPS/TAPS-v2/EM-seq/Cabernet has this layout:

```text
[barcode2][linker-BC][barcode1][insert-left][genomic insert]
```

The default whitelist has 50 barcodes of 8 bp. The default linker-BC is
`ATCCACGTGCTTGAGAGGCCAGAGCATTCG` (30 bp), and insert-left is
`CATCGGCGTACGACTAGATGTGTATAAGAGACAG` (34 bp). These anchors delimit the
barcodes and genomic insert; the non-genomic prefix is removed.

SmC uses 11 bp barcodes and different anchors:

```text
[fixed leader][barcode2][linker-BC][barcode1][linker2/insert-left][genomic insert]
```

Its default whitelist has 96 entries. SmC reads are classified as Watson or
Crick from linker conversion evidence and aligned separately. See
[SmC library structure](smc-technical.md) for the anchor sequences and mate
assignments.


| Assay | Converted C-to-T evidence | Retained C evidence |
|---|---|---|
| TAPS / TAPS-v2 | Methylated | Unmethylated |
| EM-seq / Cabernet / SmC | Unmethylated | Methylated |

## Configuration parameters

Defaults below come from [dbitm.config.example.sh](../config/dbitm.config.example.sh).
Edit the local configuration or select another file with `--config`. An empty
reference path must be set when its corresponding alignment/calling target is
used. Sequence defaults are listed in the library structure sections.

### Execution and references

| Parameter | Default | Meaning |
|---|---|---|
| `RUN_MODE` | `hpc` | `local`: run directly; `hpc`: submit Slurm jobs with dependencies. |
| `SCRATCH_ROOT` | Empty | Absolute scratch root for `barcode`, `pool`, `call`, and `methscan`; empty uses the sample directory. Results are copied back. |
| `BWA_INDEX` | Empty | Host BWA index prefix for TAPS/TAPS-v2. |
| `BWA_SPIKE_IN_INDEXES` | Empty `lambda`, `puc19` entries | Associative array mapping spike-in names to BWA index prefixes. |
| `BISCUIT_REFERENCE` | Empty | BISCUIT-indexed host FASTA for EM-seq/Cabernet/SmC. |
| `BISCUIT_SPIKE_IN_INDEXES` | Empty `lambda`, `puc19` entries | Associative array mapping spike-in names to BISCUIT-indexed FASTAs. |
| `CALL_REFERENCE` | Empty | Host calling FASTA; requires a `.fai` index. |
| `CALL_SPIKE_IN_REFERENCES` | Empty `lambda`, `puc19` entries | Associative array mapping spike-in names to calling FASTAs with `.fai` indexes. |


### Barcodes and SmC classification

| Parameter | Default | Meaning |
|---|---|---|
| `BARCODE_WHITELIST` | Empty | Standard-assay whitelist; empty selects `docs/barcodes/barcodes50.tsv`. |
| `BARCODE_CHUNK` | `10` | Number of barcode worker processes and output FASTQ chunks. |
| `BARCODE_BATCH_SIZE` | `50000` | Read pairs dispatched to a barcode worker at a time. |
| `BARCODE_COMPRESSION_STEP` | `python` | Compress during Python output, or use `shell` for subsequent gzip compression. |
| `BARCODE_GZIP_LEVEL` | `1` | gzip compression level, `0`–`9`. |
| `BARCODE_LINKER_BC` | Standard linker-BC | Anchor between barcode2 and barcode1. |
| `BARCODE_INSERT_LEFT` | Standard insert-left | Anchor defining the start of the genomic insert. |
| `BARCODE_METHYLATED_C_POSITIONS` | `3,6,10` | Zero-based methylated C positions in insert-left for TAPS-v2 conversion QC. |
| `BARCODE_LINKER_EDIT_DISTANCE` | `1` | Linker matching error allowance; also the standard TAPS insert-left allowance and SmC anchor allowance. |
| `BARCODE_HAMMING_DISTANCE` | `1` | Maximum substitutions for unique nearest-whitelist barcode correction. |
| `BARCODE_INSERT_LEFT_EDIT_DISTANCE` | `1` | Non-C mismatch allowance for C/T-aware insert-left matching in TAPS-v2/EM-seq/Cabernet. |
| `BARCODE_PROGRESS_READS` | `1000000` | Progress log interval in read pairs per worker; `0` disables reporting. |
| `SMC_BARCODE_WHITELIST` | Empty | SmC whitelist; empty selects `docs/barcodes/barcodes-smc.tsv`. |
| `SMC_LINKER_BC` | SmC linker-BC | SmC anchor between the two barcodes. |
| `SMC_INSERT_LEFT` | SmC insert-left | SmC anchor defining the genomic insert boundary. |
| `SMC_MAX_CONVERSION_MISMATCHES` | `4` | Maximum C/T-position mismatches allowed for one SmC conversion class. |
| `SMC_MINIMUM_SCORE_MARGIN` | `3` | Minimum difference between C-retained and C-to-T evidence for a confident class. |
| `SMC_SAVE_UNINFORMATIVE_FASTQ` | `0` | `1` saves ambiguous/discarded FASTQs; counts are always recorded. |

### Alignment and pooling

| Parameter | Default | Meaning |
|---|---|---|
| `ALIGN_THREADS_PER_CHUNK` | `24` | Aligner threads per host FASTQ chunk. |
| `SPIKE_ALIGN_THREADS_PER_CHUNK` | `8` | Aligner threads per spike-in FASTQ chunk. |
| `BISCUIT_DIRECTIONAL_MODE` | `1` | EM-seq/Cabernet library mode: `1` directional, `0` non-directional. SmC always uses `1`. |
| `POOL_SORT_MEM` | `12G` | Memory per SAMtools sort thread; empty derives it from `POOL_MEM / POOL_THREADS`. |

### M-bias and calling

| Parameter | Default | Meaning |
|---|---|---|
| `MBIAS_MODE` | `all` | Analyze `host`, configured `spike` targets, or `all`. |
| `MBIAS_HOST_SUBSAMPLE_FRACTION` | `0.1` | Host record sampling fraction; applied independently to SmC Watson/Crick BAMs. Spike-ins use all records. |
| `MBIAS_HOST_MAX_RECORDS` | `10000000` | Maximum sampled host records per target; `0` removes the cap. |
| `MBIAS_MAX_CYCLE` | `150` | Maximum read cycle included in M-bias analysis. |
| `MBIAS_MIN_CYCLE_COVERAGE` | `500` | A cycle needs more CpG observations than this to be included. |
| `MBIAS_R1_ORIGINAL_LENGTH` | `150` | Original read length used for barcode-mate cycle offsets and both cutoff axes; `0` infers axes from observed cycles. |
| `MBIAS_MIN_BASE_QUALITY` | `30` | Minimum base quality for M-bias evidence. |
| `MBIAS_MIN_MAPPING_QUALITY` | `10` | Minimum mapping quality for M-bias evidence. |
| `MBIAS_SAMPLING_SEED` | `42` | Seed for reproducible host sampling. |
| `MBIAS_CUTOFF_RATE_TOLERANCE` | `0.05` | Maximum absolute methylation-rate deviation from the central baseline for a stable cycle. |
| `CALL_CHROMOSOMES` | `chr1`–`chr19`, `chrX` | Comma-separated host chromosomes to call; also limits host SmC filtering. Match the reference names. |
| `CALL_MITO_CHROMOSOMES` | `chrM` | Comma-separated mitochondrial reference names. |
| `CALL_CONTEXT_MODE` | `both` | `cg`: CpG; `ch`: CA/CC/CT; `both`: all four. Also selects summary and MethSCAn contexts. |
| `SPIKE_CALL_MODE` | `all` | Call `mito`, configured `spike` targets, or `all`. |
| `CALL_MIN_BASE_QUALITY` | `MBIAS_MIN_BASE_QUALITY` | Minimum base quality for methylation calling. |
| `CALL_MIN_MAPPING_QUALITY` | `MBIAS_MIN_MAPPING_QUALITY` | Minimum mapping quality for calling and summary BAM metrics. |
| `CALL_MAX_DEPTH` | `1000000` | Maximum pileup depth per genomic position. |
| `CALL_BATCH_SIZE` | `5000000` | Maximum calling interval length in bp; reduced for short chromosomes to occupy workers. |
| `CALL_JOBS` | `16` | Worker processes within a host chromosome during calling. |

### Saturation and spatial geometry

| Parameter | Default | Meaning |
|---|---|---|
| `SATURATION_READS_THRESHOLD` | `1000000` | Fallback reads-per-spot cutoff when automatic threshold inference fails. |
| `SATURATION_PRED_FRACTION` | `2.0` | Sequencing-depth multiplier used to predict unique CpGs. |
| `SATURATION_LINEAR_R2_THRESHOLD` | `0.99` | R² threshold for treating the subsampling curve as linear. |
| `SUMMARY_FRAME_SPOT_LENGTH` | `50` | Spot side length in micrometers for the registration frame and image tool. |
| `SUMMARY_FRAME_INTERVAL` | `50` | Gap between neighboring spots in micrometers. |
| `SUMMARY_FRAME_PIXEL_LENGTH` | `0.294` | Image resolution in micrometers per pixel. |


### MethSCAn and paired TAPS

| Parameter | Default | Meaning |
|---|---|---|
| `METHSCAN_CHUNKSIZE` | `10000000` | Coverage rows read per preparation chunk. |
| `METHSCAN_CG_MIN_SITES` | Empty | Minimum observed CG sites per spot to retain it; empty disables CG processing. |
| `METHSCAN_CA_MIN_SITES` | Empty | Same setting for CA. |
| `METHSCAN_CC_MIN_SITES` | Empty | Same setting for CC. |
| `METHSCAN_CT_MIN_SITES` | Empty | Same setting for CT. |
| `PAIRED_TAPS_BIN_SIZE` | `2000` | Genomic bin length in bp. |
| `PAIRED_TAPS_MIN_COVERAGE_PER_SITE` | `1` | Minimum read coverage for an eligible site. |
| `PAIRED_TAPS_MAX_COVERAGE_PER_SITE` | Empty | Maximum eligible site coverage; empty applies no upper limit. |
| `PAIRED_TAPS_MIN_COVERAGE_PER_BIN` | `1` | Minimum total coverage per spot/bin for a usable observation. |
| `PAIRED_TAPS_MIN_SITE_PER_BIN` | `1` | Minimum contributing sites per spot/bin. |
| `PAIRED_TAPS_MIN_VALID_SPOTS_PER_BIN` | `6` | Minimum valid spots before testing bin variability. |
| `PAIRED_TAPS_VARIABLE_FRACTION` | `0.02` | Fraction of eligible bins retained by residual variance; `1` retains all eligible bins. |
| `PAIRED_TAPS_SHRINKAGE_LAMBDA` | `1.0` | Strength of coverage-aware residual shrinkage; `0` disables shrinkage. |
| `PAIRED_TAPS_CHUNKSIZE` | `1000000` | Coverage rows read per chunk. |


### Threads and Slurm resources

Each prefix in the table below has five parameters:

| Parameter pattern | Meaning |
|---|---|
| `<PREFIX>_NAME` | Slurm job name; defaults to the stage name in the first column below. |
| `<PREFIX>_THREADS` | CPUs requested per Slurm job; also controls processing threads/workers where supported. |
| `<PREFIX>_PARTITION` | Slurm partition; empty uses the cluster default. |
| `<PREFIX>_MEM` | Total memory requested per Slurm job, such as `16G`. |
| `<PREFIX>_TIME` | Slurm wall-time limit in `HH:MM:SS`. |

| Stage / default name | Prefix | Default CPUs | Default memory | Default time |
|---|---|---|---|---|
| `fastp` | `FASTP` | `8` | `16G` | `24:00:00` |
| `barcode` | `BARCODE` | `BARCODE_CHUNK` | `32G` | `48:00:00` |
| `align` | `ALIGN` | `ALIGN_THREADS_PER_CHUNK` | `64G` | `96:00:00` |
| `align-prepare` | `ALIGN_PREPARE` | `1` | `1G` | `00:30:00` |
| `spike-align` | `SPIKE_ALIGN` | `SPIKE_ALIGN_THREADS_PER_CHUNK` | `16G` | `48:00:00` |
| `spike-align-prepare` | `SPIKE_ALIGN_PREPARE` | `1` | `1G` | `00:30:00` |
| `pool` | `POOL` | `4` | `64G` | `48:00:00` |
| `mbias` | `MBIAS` | `1` | `16G` | `24:00:00` |
| `smc-filter` | `SMC_FILTER` | `8` | `32G` | `24:00:00` |
| `call` | `CALL` | `CALL_JOBS` | `96G` | `24:00:00` |
| `spike-call` | `SPIKE_CALL` | `1` | `32G` | `24:00:00` |
| `saturation` | `SATURATION` | `8` | `16G` | `24:00:00` |
| `summary` | `SUMMARY` | `8` | `16G` | `24:00:00` |
| `methscan` | `METHSCAN` | `10` | `64G` | `24:00:00` |
| `image` | `IMAGE` | `1` | `32G` | `12:00:00` |
| `paired-taps` | `PAIRED_TAPS` | `8` | `64G` | `24:00:00` |


| Parameter | Default | Meaning |
|---|---|---|
| `SBATCH_OUTPUT` | `%x_%j.out` | Slurm stdout filename pattern; `%x` is the job name and `%j` the job ID. |
| `SBATCH_ERROR` | `%x_%j.err` | Slurm stderr filename pattern. |
| `SBATCH_REQUEUE` | `true` | Allow Slurm requeueing when enabled. |
