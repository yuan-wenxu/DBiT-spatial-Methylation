# DBiT-spatial-Methylation Technical Reference

## 1. Scope

DBiT-spatial-Methylation is a paired-end sequencing workflow for spatial DNA
methylation assays. It extracts spatial barcodes from read 1, aligns host and
spike-in reads, measures cycle-dependent methylation bias, calls methylation at
CG, CA, CC, and CT sites, and creates per-spot and sample-level quality-control
summaries.

The workflow supports four assays:

- TAPS
- TAPS-v2
- EM-seq
- Cabernet

This document describes the complete pipeline implemented by the scripts in
`script/`, with `config/dbitm.config.example.sh` as the configuration reference.
When this document and the source code disagree, the source code is the final
authority.

## 2. Pipeline architecture

The complete execution graph is:

```text
fastp -> barcode -> +-> align ------+
                    |               |
                    +-> spike-align +-> pool -> mbias -> call -> spike-call -> saturation -> summary -> methscan
```

`align` and `spike-align` can run concurrently. All later stages are ordered by
data dependency.

| Stage | Shell entry point | Main purpose |
|---|---|---|
| `fastp` | `script/steps/01.fastp.sh` | Filter paired FASTQ reads |
| `barcode` | `script/steps/02.barcode.sh` | Extract and correct spatial barcodes |
| `align` | `script/steps/03.align.sh` | Align barcoded host reads and attach `CB` tags |
| `spike-align` | `script/steps/03.spike_align.sh` | Align candidate spike-in reads |
| `pool` | `script/steps/04.pool.sh` | Merge, coordinate-sort, and index BAM files |
| `mbias` | `script/steps/05.mbias.sh` | Estimate cycle-specific methylation bias and trim cutoffs |
| `call` | `script/steps/06.call.sh` | Call per-spot host CG, CA, CC, and CT methylation |
| `spike-call` | `script/steps/06.spike_call.sh` | Call aggregate mitochondrial and spike-in methylation |
| `saturation` | `script/steps/07.saturation.sh` | Estimate host CpG saturation across spots |
| `summary` | `script/steps/08.summary.sh` | Produce QC tables and plots |
| `methscan` | `script/steps/09.methscan.sh` | Detect context-specific VMRs and construct spot-by-VMR matrices |

## 3. Installation and runtime environment

The project uses Pixi for dependency and environment management. Important
runtime packages include Python 3.11, BWA, BISCUIT, SAMtools, Sinto, MethSCAn,
pysam, fastp, pigz, fuzzysearch, and Matplotlib.

Initialize the locked environment and install the command-line entry point with:

```bash
./install-cli.sh
```

By default, this creates a `dbitm` symlink under `~/.local/bin`. Set
`DBITM_INSTALL_DIR` before running the installer to select another destination.
The pipeline itself invokes tools through the Pixi environment.

## 4. Command-line interface

The general syntax is:

```bash
dbitm <assay> <step> --input /path/to/sample/fastq [options]
```

Valid assays are `taps`, `taps-v2`, `emseq`, and `cabernet`. A stage name from
the table above or `all` can be supplied as `<step>`.

Examples:

```bash
dbitm taps all --input /data/sample/fastq
dbitm emseq mbias --input /data/sample/fastq
dbitm cabernet call --input /data/sample/fastq --dry-run
dbitm taps summary --input /data/sample/fastq --config config/custom.sh
dbitm taps methscan --input /data/sample/fastq --config config/custom.sh
```

The input directory must contain exactly one paired FASTQ set whose filenames
can be distinguished by independent R1 and R2 tokens. For an input directory
`/data/sample/fastq`, the workflow writes results beneath:

```text
/data/sample/dbitm
```

### 4.1 Local and HPC execution

`RUN_MODE=local` runs stages directly. `RUN_MODE=hpc` submits stages through
Slurm using `sbatch`; the example configuration selects HPC mode. In HPC mode,
the controller builds job dependencies so that host and spike alignment can run
in parallel and downstream stages wait for all required jobs.

Each stage has configurable Slurm job name, CPU count, memory, time, and an
optional partition. The main controller validates the required resource fields
before submission.

### 4.2 Dry-run behavior

All stages and the main controller accept `--dry-run`. In local mode, dry-run is
passed to the selected shell stage. In HPC mode, the controller prints the
`sbatch` commands and dependencies without submitting jobs. Dry-run is intended
to validate resolved paths, configuration, and command construction; it does not
validate biological inputs by executing the underlying tools.

## 5. Configuration model

Copy and edit the example configuration:

```bash
cp config/dbitm.config.example.sh config/dbitm.config.sh
```

The major configuration groups are described below.

### 5.1 Reference indexes

TAPS and TAPS-v2 alignment use:

- `BWA_INDEX` for the host genome
- `BWA_SPIKE_IN_INDEXES` for named spike-in genomes

EM-seq and Cabernet alignment use:

- `BISCUIT_REFERENCE` for the host genome
- `BISCUIT_SPIKE_IN_INDEXES` for named spike-in genomes
- `BISCUIT_DIRECTIONAL_MODE=1` for directional libraries or `0` for
  non-directional libraries

Methylation calling and M-bias analysis use the unconverted FASTA references:

- `CALL_REFERENCE` for host and mitochondrial calling
- `CALL_SPIKE_IN_REFERENCES` for named spike-ins

Each calling FASTA must have a SAMtools `.fai` index. A single FASTA can be used
as both the BISCUIT reference and calling reference when its required indexes
have been built.

### 5.2 Barcode configuration

Barcode configuration defines:

- the barcode whitelist
- barcode length
- expected linker sequence
- expected insert-left sequence
- linker and insert matching thresholds
- barcode correction distance
- chunk size and worker count
- TAPS-v2 conversion-QC positions

The supplied default whitelist is `docs/barcodes/barcodes50.tsv`. When whitelist
indices are absent, the extractor assigns zero-based, zero-padded indices such
as `00` through `49`.

### 5.3 Calling and quality thresholds

The main calling controls include:

- `CALL_CONTEXT_MODE`: `cg`, `ch`, or `both`
- `CALL_CHROMOSOMES`: host chromosomes to call; host M-bias is sampled from these contigs only
- `CALL_MITO_CHROMOSOMES`: mitochondrial contigs to call
- `SPIKE_CALL_MODE`: `all`, `mito`, or `spike`
- minimum mapping and base qualities
- pileup maximum depth
- maximum genomic interval size; shorter chromosomes are divided into at
  least `CALL_JOBS` intervals when their length permits
- parallel caller job count

### 5.4 Saturation configuration

`SATURATION_READS_THRESHOLD` selects high-read spots. The stage predicts at
`SATURATION_PRED_FRACTION` and uses `SATURATION_LINEAR_R2_THRESHOLD` to choose
between linear and saturation-curve extrapolation. `SATURATION_*` Slurm fields
control its standalone job.

### 5.5 MethSCAn configuration

The independent MethSCAn stage uses:

- `METHSCAN_CHUNKSIZE`: chromosome chunk size used by `methscan prepare`
- `METHSCAN_CG_MIN_SITES`: minimum observed CG sites required per spot
- `METHSCAN_CA_MIN_SITES`: minimum observed CA sites required per spot
- `METHSCAN_CC_MIN_SITES`: minimum observed CC sites required per spot
- `METHSCAN_CT_MIN_SITES`: minimum observed CT sites required per spot
- `METHSCAN_THREADS`: threads used by `methscan scan` and `methscan matrix`

`METHSCAN_NAME`, `METHSCAN_THREADS`, `METHSCAN_PARTITION`, `METHSCAN_MEM`, and
`METHSCAN_TIME` control the separate Slurm job. The configured CPU count is
also passed to the two MethSCAn commands that support parallel execution.

## 6. Read and coordinate data model

### 6.1 R1 structure and spatial barcode

The expected R1 layout is:

```text
[barcode2][linker][barcode1][insert-left][biological insert]
```

The barcode stage identifies the linker and insert-left sequence, extracts the
two adjacent barcodes, corrects them against the whitelist, and retains only the
biological insert in R1. R2 is not structurally trimmed by this stage.

The spatial identifier is represented as:

```text
barcode1+barcode2
```

The barcode stage places this identifier before the original read name. Host
alignment transfers it into the BAM `CB` tag. The methylation caller accepts the
delimiter form and the legacy concatenated form.

### 6.2 Read orientation and sequencing cycle

BAM query sequences are stored in the orientation used by the SAM/BAM record.
M-bias converts the aligned query position into a cycle measured from the
original molecule's 5-prime end:

- forward alignments use the query position directly
- reverse alignments invert the query position
- R1 receives an offset for the barcode/linker portion removed before alignment
- R2 has no barcode-removal offset

Cycles are one-based. `MBIAS_R1_ORIGINAL_LENGTH` supplies the original R1 length
used to reconstruct the removed prefix. Setting it to zero makes R1 cycles
relative to the trimmed biological sequence instead. Calling applies the
inferred end trims with reverse-read orientation accounted for.

### 6.3 Genomic coordinates in coverage files

Coverage output uses zero-based genomic coordinates. The `start` and `end`
columns both contain the selected cytosine-representation coordinate; they are
not a BED half-open interval.

For CpG calls, evidence from both strands is merged at the forward-strand
reference C.

CH calls retain strand: plus-strand sites use the reference C coordinate and
minus-strand sites use the corresponding reference G coordinate.

## 7. Assay-specific methylation chemistry

The callers inspect aligned dinucleotides to classify CpG evidence consistently
on both strands.

Methylation interpretation is assay dependent:

| Assay | `TG` or `CA` | `CG` |
|---|---|---|
| TAPS / TAPS-v2 | methylated | unmethylated |
| EM-seq / Cabernet | unmethylated | methylated |

Only the primary paired-end orientation flags `83`, `99`, `147`, and `163` are
accepted for methylation classification. For an observed `CG`, flags `99` and
`147` represent the top/plus family, whereas flags `83` and `163` represent the
bottom/minus family. CpG output does not retain that family label because both
families update the same forward-C counter.

For non-CpG CH contexts (`CA`, `CC`, and `CT` in the reference):

- plus-strand evidence uses observed C or T
- minus-strand evidence uses observed G or A
- retained C/G is methylated for EM-seq and Cabernet
- converted T/A is methylated for TAPS and TAPS-v2

## 8. Stage 1: FASTQ filtering

Entry point: `script/steps/01.fastp.sh`

The stage verifies that a single R1/R2 pair can be identified and runs fastp on
the paired reads. Adapter trimming is disabled because the downstream barcode
structure is assay-specific.

Outputs under `dbitm/fastp/` are:

```text
R1.filtered.fastq.gz
R2.filtered.fastq.gz
fastp.html
fastp.json
fastp.log
```

When scratch storage is enabled, raw inputs are copied to the stage scratch
directory. Successful output is installed at completion; if the stage exits
with an error or a handled signal, any partial output is recovered first.

## 9. Stage 2: barcode extraction

Entry points:

- `script/steps/02.barcode.sh`
- `script/steps/python/02.extract_bc.py`

The extractor searches R1 for the expected library structure, corrects barcode1
and barcode2 against the whitelist using the configured Hamming-distance limit,
adds the corrected barcode pair to the read name, and removes the non-biological
R1 prefix.

Matching differs by assay:

- TAPS uses standard linker and insert-left matching.
- TAPS-v2 uses standard linker matching and C/T-aware insert-left matching. It
  also evaluates conversion at configured zero-based cytosine positions.
- EM-seq and Cabernet use C/T-aware matching for both linker and insert-left;
  whitelist sequences are converted from C to T for matching.

Reads that do not satisfy the DBiT structure or barcode criteria are written to
the `spike-in` FASTQ streams. These files contain spike-in candidates; the
barcode stage does not itself prove their reference of origin.

Representative outputs under `dbitm/barcode/` are:

```text
0001.R1.demux.fastq.gz
0001.R2.demux.fastq.gz
0001.R1.spike-in.fastq.gz
0001.R2.spike-in.fastq.gz
0001.stats.json
stats.json
barcode.log
```

The input is processed in chunks and batches. Configured worker processes handle
batches in parallel, and output compression can use either Python gzip support
or an external compressor.

## 10. Stage 3A: host alignment

Entry point: `script/steps/03.align.sh`

The stage aligns each demultiplexed host FASTQ chunk:

- TAPS and TAPS-v2 use `bwa mem`.
- EM-seq and Cabernet use `biscuit align` with
  `BISCUIT_DIRECTIONAL_MODE`.

The alignment stream passes through `sinto nametotag`, which moves the
`barcode1+barcode2` read-name prefix into the `CB` BAM tag.

Outputs under `dbitm/align/` include one BAM and associated logs per chunk:

```text
0001.cb.bam
0001.align.log
```

Each BAM is checked with `samtools quickcheck`. Chunk BAMs are not required to
be coordinate sorted or indexed at this stage. Chunks are parallelized; when
planning HPC resources, the effective CPU demand may approach the number of
chunks running simultaneously multiplied by threads per chunk.

## 11. Stage 3B: spike-in alignment

Entry point: `script/steps/03.spike_align.sh`

Each `spike-in` FASTQ pair is aligned independently to every configured spike-in
reference. The aligner follows the same assay rule as host alignment. Spatial
`CB` tags are not required because downstream spike-in calling is aggregate,
not per spot.

Outputs under `dbitm/spike_align/` include:

```text
0001.lambda.bam
0001.lambda.flagstat.txt
0001.lambda.align.log
0001.puc19.bam
0001.puc19.flagstat.txt
0001.puc19.align.log
```

The spike alignment stage accepts filename-safe configured spike-in names. A
full pipeline run requires the relevant spike reference arrays to be populated.

## 12. Stage 4: BAM pooling

Entry point: `script/steps/04.pool.sh`

Host chunk BAMs are concatenated with `samtools cat`, coordinate sorted with
`samtools sort`, and indexed. Spike-in BAMs are pooled independently using the
same operations.

Outputs under `dbitm/pooled/` are:

```text
pooled.cb.bam
pooled.cb.bam.bai
pooled.lambda.bam
pooled.lambda.bam.bai
pooled.puc19.bam
pooled.puc19.bam.bai
pool.log
```

`POOL_SORT_MEM` is memory per SAMtools sort thread, not total stage memory.

Current limitation: spike alignment can accept arbitrary configured spike names,
but the pooling stage currently enumerates `lambda` and `puc19`. Supporting
additional names requires extending the pooling list.

## 13. Stage 5: M-bias analysis and cutoff inference

Entry points:

- `script/steps/05.mbias.sh`
- `script/steps/python/05.mbias.py`
- `script/steps/python/05.mbias_cutoff.py`

### 13.1 Target selection and filtering

`MBIAS_MODE` selects `host`, `spike`, or `all`. Host analysis is restricted to
the `CALL_CHROMOSOMES` contigs, so reads mapped outside that list — including
mitochondrial reads present in the pooled host BAM — never enter host M-bias
statistics; the stage fails fast when `CALL_CHROMOSOMES` is unset or empty in
`host` or `all` mode. Spike-in targets scan all contigs of their own references
without this restriction. Host analysis uses a deterministic pseudo-random read
fraction and an optional record cap; the example configuration uses a 0.1 host
fraction and a maximum of 10,000,000 records.
Spike-in targets are analyzed without host-style subsampling.

M-bias excludes unmapped, secondary, supplementary, unsupported paired-end
orientations, and observations below the configured mapping or base quality.
The same assay-specific CpG dinucleotide rules described in Section 7 are used
to avoid a mismatch between M-bias and methylation calling. M-bias only evaluates
CG and does not generate separate CA, CC, or CT profiles.

### 13.2 Per-cycle statistics

For each R1/R2 cycle, the program records methylated, unmethylated, and total
CpG observations. Plus- and minus-family CpG evidence is combined in the same
counter, using the caller-aligned query position without a one-base shift.
Cycles whose total coverage does not exceed
`MBIAS_MIN_CYCLE_COVERAGE` (default 500) are skipped and never reach the TSV,
PNG, or cutoff inference. TSV rates are fractions from 0 to 1; plots use
percentages from 0 to 100.

Each selected target produces:

```text
<label>.mbias.tsv
<label>.mbias.png
<label>.mbias.cutoffs.tsv
```

The combined stage log is `dbitm/mbias/mbias.log`.

### 13.3 Cutoff inference

The cutoff program smooths the methylation-rate profile, estimates a central
baseline, rejects poorly covered cycles, and searches from both ends for a
stable run within the configured tolerance. The default controls are:

- methylation-rate tolerance: `MBIAS_CUTOFF_RATE_TOLERANCE` (default 0.05)
- minimum coverage relative to the central region: 0.10
- absolute minimum coverage: 100
- smoothing window: 5 cycles
- stable run length: 5 cycles

The resulting table contains left and right trim values for R1 and R2. Existing
complete TSV/PNG outputs can be reused, while incomplete output sets are treated
as errors to prevent mixing stale and current results. When a target has no
usable cycles, or no stable region within tolerance can be found for either
read, inference is skipped instead of failing and the target receives a
zero-byte cutoff marker; host, mitochondrial, and spike-in calling proceed
without trimming.

## 14. Stage 6A: per-spot host methylation calling

Entry points:

- `script/steps/06.call.sh`
- `script/steps/python/06.methy_caller.py`

The caller reads the coordinate-sorted host BAM, reference FASTA, barcode
whitelist, and the host M-bias cutoff (`host.mbias.cutoffs.tsv`). It processes
configured host chromosomes in genomic intervals. Interval workers can run in
parallel, while chromosome dispatch remains controlled by the shell stage.

### 14.1 Trim selection

The host shell stage reads only `host.mbias.cutoffs.tsv` and passes its four
trim values to the per-spot host caller. Mitochondrial and spike-in cutoff
tables are consumed separately by the Stage 6B aggregate caller. A zero-byte
host cutoff table means no suitable cutoff was inferred; the host caller then
runs with zero trimming on all read ends.

Trim decisions use molecular 5-prime and 3-prime positions and therefore account
for reverse-aligned reads. They do not use raw BAM query positions as if every
read were forward aligned.

### 14.2 Spot manifest

The caller builds `dbitm/coverage/spot_manifest.tsv` from every whitelist pair,
including spots with no called sites. Spot coordinates are zero-based barcode
indices.

### 14.3 CpG strand handling

The caller combines plus- and minus-strand CpG evidence at the forward reference
C coordinate for every assay. It does not emit a second row at the symmetric G
coordinate and has no separate CpG-strand mode.

### 14.4 Per-spot output

Host coverage files are grouped by the first spatial coordinate:

```text
coverage/host/<X>/<X>_<Y>.CG.cov
coverage/host/<X>/<X>_<Y>.CA.cov
coverage/host/<X>/<X>_<Y>.CC.cov
coverage/host/<X>/<X>_<Y>.CT.cov
```

CG files contain:

```text
chrom  start  end  methylation_percent  methylated_count  unmethylated_count
```

CA, CC, and CT files add `context` and `strand` columns and do not combine the
two strands. A context file is created only when that spot has at least one row;
all whitelist spots remain in `spot_manifest.tsv`.

## 15. Stage 6B: mitochondrial and spike-in calling

Entry points:

- `script/steps/06.spike_call.sh`
- `script/steps/python/06.spike_caller.py`

This is a separate stage and does not modify the per-spot host caller. It creates
aggregate, non-spatial calls for mitochondrial and spike-in references.

`SPIKE_CALL_MODE` controls target selection:

- `mito`: call configured mitochondrial contigs from the host BAM
- `spike`: call configured spike-in BAMs
- `all`: call both groups

Mitochondrial calling uses `pooled.cb.bam`, `CALL_REFERENCE`, and
`CALL_MITO_CHROMOSOMES`. Spike-in calling uses each `pooled.<name>.bam`, its
configured FASTA, and the contigs listed in that FASTA's `.fai` index.

The host caller uses the host M-bias cutoff. Each aggregate target uses its own
M-bias cutoff: mitochondrial calls use the host cutoff and each spike-in uses
its matching cutoff.

Outputs under `dbitm/coverage/` include:

```text
host_mito.CG.cov
host_mito.CA.cov
host_mito.CC.cov
host_mito.CT.cov
lambda.CG.cov
lambda.CA.cov
lambda.CC.cov
lambda.CT.cov
puc19.CG.cov
puc19.CA.cov
puc19.CC.cov
puc19.CT.cov
spike_call.log
```

These files use the host calling rules: CG is strand-merged, while CA, CC, and
CT retain context and strand.

## 16. Stage 7: CpG saturation

Entry points:

- `script/steps/07.saturation.sh`
- `script/steps/python/07.saturation.py`

The stage counts reads per spot from `pooled/pooled.cb.bam`, selects spots above
`SATURATION_READS_THRESHOLD`, and reads their
`coverage/host/**/*.CG.cov` files. For each subsampling fraction it estimates
the expected unique CpGs from per-site depth, summarizes spots by median and
IQR, then fits a through-origin linear model and an exponential saturation
curve. A strong linear fit is reported as unsaturated with `saturation_rate=NA`.

Outputs are:

```text
saturation/
├── saturation_curve.png
├── saturation_summary.tsv
└── saturation.log
```

The summary table reports observed and predicted median unique CpGs, the fitted
model, saturation rate, and number of high-read spots. If no usable spots exist,
the stage writes valid outputs containing `NA` values.

## 17. Stage 8: summary and visualization

Entry points:

- `script/steps/08.summary.sh`
- `script/steps/python/08.summary.py`

The summary stage combines barcode statistics, fastp metrics, pooled BAM
metrics, per-spot host coverage, and aggregate mitochondrial/spike-in coverage.
It follows `CALL_CONTEXT_MODE`: `cg` summarizes CG files, `ch` summarizes CA,
CC, and CT files separately, and `both` summarizes all four contexts. It does
not require or emit a sample-ID field.

### 17.1 Inputs

The principal inputs are:

- `fastp/fastp.json`
- `barcode/stats.json`
- `pooled/pooled.cb.bam`
- `coverage/spot_manifest.tsv`
- `coverage/host/*/*.<context>.cov`
- `coverage/host_mito.<context>.cov`, when produced
- `coverage/<spike>.<context>.cov`, when produced
- corresponding pooled spike-in BAMs
- `saturation/saturation_summary.tsv`

Here, `<context>` is `CG`, `CA`, `CC`, or `CT` according to
`CALL_CONTEXT_MODE`.

### 17.2 Per-spot metrics

`per_spot_summary.tsv` contains one row for every manifest spot, including spots
with zero reads or zero called sites. It always reports spatial indices and read
counts. In CG mode it adds `mean_methylation` and `cpg_site_count`. In CH mode it
adds separate `mean_<context>_methylation` and `<context>_site_count` fields for
CA, CC, and CT, where `<context>` is lowercase. Both mode includes all eight
context-specific fields.

The per-spot `reads` field counts BAM records assigned by `CB`. Each
`*_site_count` is the number of nonzero coverage rows, not read depth.

The mean is:

```text
mean methylation = sum of per-site methylation percentages / number of context rows
```

Each CG, CA, CC, or CT row therefore has equal weight within its context. This
is not a read-depth-weighted ratio of total methylated observations to total
observations.

### 17.3 Sample-level metrics

`sample_summary.tsv` reports workflow-wide values including:

- raw reads from the fastp report
- barcode-retained reads and their fraction of raw reads
- host and spike-in mapped alignment records
- host alignment records assigned to valid paired-end orientations
- host mean methylation and median per-spot site count for each selected context
- mitochondrial and spike-in mean methylation for each selected context
- `saturation_rate` from Stage 7, or `NA` when unavailable

Within each context, the host-wide mean weights each spot mean by its site count,
which is equivalent to giving every spot-by-site output row equal weight.
Missing mitochondrial or configured spike-in coverage files produce `NA` for
that target and context instead of aborting the summary.

Barcode-kept records are paired reads, so the reported kept-read count is twice
the retained-pair count. BAM records are counted as individual alignments.

For mapped-read metrics, the summary requires a mapped record at or above the
mapping-quality threshold and accepts `NH` when it is absent or no greater than
one. This metric does not independently remove secondary or supplementary
records; interpret it as the implementation's QC count rather than a guaranteed
count of primary fragments.

### 17.4 Plots

The reads heatmap is always generated:

```text
reads_heatmap.png
```

CG mode additionally generates:

```text
cpg_site_count_heatmap.png
mean_methylation_heatmap.png
mean_methylation_violin.png
```

CH mode additionally generates three plot groups:

```text
ca_site_count_heatmap.png
mean_ca_methylation_heatmap.png
mean_ca_methylation_violin.png
cc_site_count_heatmap.png
mean_cc_methylation_heatmap.png
mean_cc_methylation_violin.png
ct_site_count_heatmap.png
mean_ct_methylation_heatmap.png
mean_ct_methylation_violin.png
```

Both mode generates all thirteen plots. Heatmap axes use the zero-based barcode
indices from the spot manifest. Each violin plot shows the distribution of mean
methylation across spots with at least one called site in that context. The
reads heatmap uses `Reds`; violin plots use an automatic y-axis range.

The complete summary output directory is:

```text
summary/
├── per_spot_summary.tsv
├── sample_summary.tsv
├── reads_heatmap.png
├── cpg_site_count_heatmap.png
├── mean_methylation_heatmap.png
├── mean_methylation_violin.png
├── ca_site_count_heatmap.png
├── mean_ca_methylation_heatmap.png
├── mean_ca_methylation_violin.png
├── cc_site_count_heatmap.png
├── mean_cc_methylation_heatmap.png
├── mean_cc_methylation_violin.png
├── ct_site_count_heatmap.png
├── mean_ct_methylation_heatmap.png
├── mean_ct_methylation_violin.png
└── summary.log
```

The CG or CA/CC/CT plot groups are omitted when that context is not selected.

## 18. Stage 9: CG, CA, CC, and CT VMR matrices

Entry point: `script/steps/09.methscan.sh`

This independent stage follows `CALL_CONTEXT_MODE`. It processes CG in `cg`
mode, processes CA, CC, and CT independently in `ch` mode, and runs all four
workflows sequentially in `both` mode after the summary stage:

```text
prepare -> filter -> smooth -> scan -> matrix
```

For each selected context, `prepare` reads every matching
`coverage/host/**/*.<CONTEXT>.cov` file as one spatial spot. `filter` retains
spots using the context-specific `METHSCAN_<CONTEXT>_MIN_SITES` threshold.
`smooth` constructs a separate pseudo-bulk methylation background for every
context, `scan` identifies VMRs, and `matrix` quantifies those regions across
spots.
MethSCAn defaults are used for smoothing and VMR detection parameters that are
not exposed in the configuration. The four contexts are never combined in one
MethSCAn data set.

The outputs are stored in a separate top-level directory:

```text
methscan/
├── CG/
├── CA/
├── CC/
└── CT/
```

Each selected context directory contains `compact/`, `filter/smoothed/`,
`VMRs.bed`, `matrix/`, and `methscan.log`. The matrix directory contains
`methylated_sites.csv.gz`, `total_sites.csv.gz`,
`methylation_fractions.csv.gz`, and `mean_shrunken_residuals.csv.gz`.

Only the selected context directories are generated. The matrices written by
MethSCAn are spot-by-VMR tables. Missing values mean that a spot has no usable
observation in that VMR; they are not zero methylation. This stage does not
perform the iterative-PCA imputation used in some downstream analyses.

## 19. Complete output layout

A successful full run normally has the following structure:

```text
dbitm/
├── fastp/
├── barcode/
├── align/
├── spike_align/
├── pooled/
├── mbias/
├── coverage/
├── saturation/
├── summary/
└── methscan/
```

Exact files depend on `MBIAS_MODE`, `CALL_CONTEXT_MODE`, `SPIKE_CALL_MODE`, and
the configured references.

## 20. Scratch storage and restart behavior

When `SCRATCH_ROOT` is configured, stages create isolated temporary workspaces
under a `dbitm` subdirectory, copy or stage the required inputs, and remove the
stage workspace through shell traps. The exact files staged depend on the
operation; reference FASTA and `.fai` files are included where random access is
required.

On successful completion, stages install their outputs and remove the scratch
workspace. On a nonzero exit or a handled `INT`, `TERM`, or `HUP` signal, every
scratch-enabled stage first copies its current output directory back to the
corresponding directory under `dbitm/`, preserves the original exit status, and
then removes the scratch workspace. If recovery copying fails, the stage reports
a warning and still removes the scratch workspace. An uncatchable `KILL` signal
or an abrupt host/container failure cannot run this recovery handler.

Recovered output may be incomplete. M-bias has explicit restart checks: a
complete existing TSV/PNG pair can be reused, whereas partial output is
rejected. Before rerunning a failed production stage, inspect its log and
distinguish incomplete output from a valid completed result.
