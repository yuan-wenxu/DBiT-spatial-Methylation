# DBiT-spatial-Methylation Technical Reference

## 1. Scope

DBiT-spatial-Methylation is a paired-end sequencing workflow for spatial DNA
methylation assays. It extracts spatial barcodes from read 1, aligns host and
spike-in reads, measures cycle-dependent methylation bias, calls methylation at
CG and CH sites, and creates per-spot and sample-level quality-control summaries.

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
                    +-> spike-align +-> pool -> mbias -> call -> spike-call -> summary
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
| `call` | `script/steps/06.call.sh` | Call per-spot host CG and CH methylation |
| `spike-call` | `script/steps/06.spike_call.sh` | Call aggregate mitochondrial and spike-in methylation |
| `summary` | `script/steps/07.summary.sh` | Produce QC tables and plots |

## 3. Installation and runtime environment

The project uses Pixi for dependency and environment management. Important
runtime packages include Python 3.11, BWA, BISCUIT, SAMtools, Sinto, pysam,
fastp, pigz, fuzzysearch, and Matplotlib.

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
- `CALL_CG_STRAND_MODE`: `separate` or `merged`
- `CALL_CHROMOSOMES`: host chromosomes to call
- `CALL_MITO_CHROMOSOMES`: mitochondrial contigs to call
- `SPIKE_CALL_MODE`: `all`, `mito`, or `spike`
- minimum mapping and base qualities
- pileup maximum depth
- genomic interval size
- parallel caller job count

Cabernet always uses `separate` CG strand mode, regardless of
`CALL_CG_STRAND_MODE`.

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
For M-bias and trimming, the scripts convert the aligned query position into a
cycle measured from the original molecule's 5-prime end:

- forward alignments use the query position directly
- reverse alignments invert the query position
- R1 receives an offset for the barcode/linker portion removed before alignment
- R2 has no barcode-removal offset

Cycles are one-based. `MBIAS_ORIGINAL_R1_LENGTH` supplies the original R1 length
used to reconstruct the removed prefix. Setting it to zero makes R1 cycles
relative to the trimmed biological sequence instead.

### 6.3 Genomic coordinates in coverage files

Coverage output uses zero-based genomic coordinates. The `start` and `end`
columns both contain the selected cytosine-representation coordinate; they are
not a BED half-open interval.

For merged CpG calls, both strands are placed at the forward-strand reference C.
For separate CpG calls, plus-strand evidence is placed at that C, while
minus-strand evidence is placed at the paired reference G and marked with its
strand. The G coordinate therefore represents the complementary-strand C.

## 7. Assay-specific methylation chemistry

The callers inspect aligned dinucleotides to classify CpG evidence consistently
on both strands.

| Observed dinucleotide | Represented strand | Cytosine-derived query base |
|---|---|---|
| `TG` | plus | first base |
| `CA` | minus | second base |
| `CG` | determined from alignment flag | C for top-family reads or G for bottom-family reads |

Methylation interpretation is assay dependent:

| Assay | `TG` or `CA` | `CG` |
|---|---|---|
| TAPS / TAPS-v2 | methylated | unmethylated |
| EM-seq / Cabernet | unmethylated | methylated |

Only the primary paired-end orientation flags `83`, `99`, `147`, and `163` are
accepted for methylation classification. For an observed `CG`, flags `99` and
`147` represent the top/plus family, whereas flags `83` and `163` represent the
bottom/minus family.

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
directory. The final output directory is installed only after successful
completion.

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

`MBIAS_MODE` selects `host`, `spike`, or `all`. Host analysis uses a
deterministic pseudo-random read fraction and an optional record cap; the example
configuration uses a 0.1 host fraction and a maximum of 10,000,000 records.
Spike-in targets are analyzed without host-style subsampling.

M-bias excludes unmapped, secondary, supplementary, unsupported paired-end
orientations, and observations below the configured mapping or base quality.
The same assay-specific CpG dinucleotide rules described in Section 7 are used
to avoid a mismatch between M-bias and methylation calling.

### 13.2 Per-cycle statistics

For each R1/R2 cycle, the program records methylated, unmethylated, and total
CpG observations. Output rates in TSV files are fractions from 0 to 1, while
plots display percentages from 0 to 100.

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

- methylation-rate tolerance: 0.02
- minimum coverage relative to the central region: 0.10
- absolute minimum coverage: 100
- smoothing window: 5 cycles
- stable run length: 5 cycles

The resulting table contains left and right trim values for R1 and R2. Existing
complete TSV/PNG outputs can be reused, while incomplete output sets are treated
as errors to prevent mixing stale and current results.

## 14. Stage 6A: per-spot host methylation calling

Entry points:

- `script/steps/06.call.sh`
- `script/steps/python/06.methy_caller.py`

The caller reads the coordinate-sorted host BAM, reference FASTA, barcode
whitelist, and M-bias cutoffs. It processes configured host chromosomes in
genomic intervals. Interval workers can run in parallel, while chromosome
dispatch remains controlled by the shell stage.

### 14.1 Trim selection

When multiple M-bias targets are selected, the host shell stage reads every
selected cutoff table and takes the maximum R1-left, R1-right, R2-left, and
R2-right trim. This conservative rule prevents a cycle considered biased in any
selected target from contributing to the host call.

Trim decisions use molecular 5-prime and 3-prime positions and therefore account
for reverse-aligned reads. They do not use raw BAM query positions as if every
read were forward aligned.

### 14.2 Spot manifest

The caller builds `dbitm/coverage/spot_manifest.tsv` from every whitelist pair,
including spots with no called sites. Spot coordinates are zero-based barcode
indices.

### 14.3 CpG strand modes

`CALL_CG_STRAND_MODE=merged` combines plus- and minus-strand evidence at the
forward reference C coordinate. `separate` retains two strand-specific records:
plus evidence at the forward C and minus evidence at the paired G. Cabernet
always uses `separate` mode.

### 14.4 Per-spot output

Host coverage files are grouped by the first spatial coordinate:

```text
coverage/host/<X>/<X>_<Y>.CG.cov
coverage/host/<X>/<X>_<Y>.CH.cov
```

CG files contain:

```text
chrom  start  end  methylation_percent  methylated_count  unmethylated_count  [strand]
```

The optional strand field is present in separate mode. CH files additionally
record context and strand. Methylation percentage is calculated from the
methylated and unmethylated observations for that output row.

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

Unlike the host caller's conservative combined cutoff, each aggregate target
uses its own M-bias cutoff: mitochondrial calls use the host cutoff and each
spike-in uses its matching cutoff.

Outputs under `dbitm/coverage/` include:

```text
host_mito.CG.cov
host_mito.CH.cov
lambda.CG.cov
lambda.CH.cov
puc19.CG.cov
puc19.CH.cov
spike_call.log
```

These files use the same assay chemistry, context selection, strand mode,
quality filters, coordinate convention, and trim logic as host calling.

## 16. Stage 7: summary and visualization

Entry points:

- `script/steps/07.summary.sh`
- `script/steps/python/07.summary.py`

The summary stage combines barcode statistics, fastp metrics, pooled BAM
metrics, per-spot host CG coverage, and aggregate mitochondrial/spike-in CG
coverage. It does not require or emit a sample-ID field.

### 16.1 Inputs

The principal inputs are:

- `fastp/fastp.json`
- `barcode/stats.json`
- `pooled/pooled.cb.bam`
- `coverage/spot_manifest.tsv`
- `coverage/host/*/*.CG.cov`
- `coverage/host_mito.CG.cov`, when produced
- `coverage/<spike>.CG.cov`, when produced
- corresponding pooled spike-in BAMs

CH files are not currently included in the summary tables or plots.

### 16.2 Per-spot metrics

`per_spot_summary.tsv` contains one row for every manifest spot, including spots
with zero reads or zero CpG sites. It reports spatial indices, read counts,
valid-orientation read counts, CpG site count, methylation sum, and mean CpG
methylation.

The mean is:

```text
mean methylation = sum of per-site methylation fractions / number of CpG rows
```

Each CpG row therefore has equal weight. This is not a read-depth-weighted ratio
of total methylated observations to total observations.

### 16.3 Sample-level metrics

`sample_summary.tsv` reports workflow-wide values including:

- fastp input and output reads
- barcode-kept read pairs and reads
- host BAM records and mapped reads
- reads assigned to valid paired-end orientations
- total host CpG rows and mean methylation
- mitochondrial and spike-in reads, CpG rows, and mean methylation when present
- saturation fields where supported

The host-wide mean weights each spot mean by its CpG row count, which is
equivalent to giving every spot-by-CpG output row equal weight.

Barcode-kept records are paired reads, so the reported kept-read count is twice
the retained-pair count. BAM records are counted as individual alignments.

For mapped-read metrics, the summary requires a mapped record at or above the
mapping-quality threshold and accepts `NH` when it is absent or no greater than
one. This metric does not independently remove secondary or supplementary
records; interpret it as the implementation's QC count rather than a guaranteed
count of primary fragments.

### 16.4 Plots

The stage creates four plots:

```text
reads_per_spot.heatmap.png
cpg_sites_per_spot.heatmap.png
mean_methylation_per_spot.heatmap.png
mean_methylation_per_spot.violin.png
```

Heatmap axes use the zero-based barcode indices from the spot manifest. The
violin plot shows the distribution of mean CpG methylation across spots with at
least one called CpG row.

The complete summary output directory is:

```text
summary/
├── per_spot_summary.tsv
├── sample_summary.tsv
├── reads_per_spot.heatmap.png
├── cpg_sites_per_spot.heatmap.png
├── mean_methylation_per_spot.heatmap.png
├── mean_methylation_per_spot.violin.png
└── summary.log
```

## 17. Complete output layout

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
└── summary/
```

Exact files depend on `MBIAS_MODE`, `CALL_CONTEXT_MODE`, `SPIKE_CALL_MODE`, and
the configured references.

## 18. Scratch storage and restart behavior

When `SCRATCH_ROOT` is configured, stages create isolated temporary workspaces
under a `dbitm` subdirectory, copy or stage the required inputs, and remove the
stage workspace through shell traps. The exact files staged depend on the
operation; reference FASTA and `.fai` files are included where random access is
required.

Most stages construct outputs in temporary locations and install them only after
successful completion. M-bias has explicit restart checks: a complete existing
TSV/PNG pair can be reused, whereas partial output is rejected. Before rerunning
a failed production stage, inspect its log and distinguish incomplete output
from a valid completed result.

## 19. Quality-control interpretation

- M-bias should be inspected separately for R1 and R2 because library structure,
  barcode removal, and end effects differ between mates.
- An inferred trim cutoff is a cycle filter, not genomic soft clipping. It is
  applied while accepting methylation observations.
- Merged CpG mode assumes complementary CpG methylation can be summarized at one
  forward-C coordinate. Separate mode preserves strand asymmetry and is required
  for Cabernet.
- Spike-in FASTQs are generated from reads rejected by DBiT barcode structure;
  their actual identity is established only by alignment to spike references.
- Per-site mean methylation and read-depth-weighted methylation answer different
  questions. The summary currently reports the former.
- Zero-based output coordinates must be converted explicitly before comparison
  with one-based formats or tools.

## 20. Known implementation constraints

1. The pool stage currently names `lambda` and `puc19` explicitly even though
   the spike alignment configuration can describe other spike-ins.
2. Summary statistics currently use CG coverage only; CH calling results are
   retained as coverage files but are not summarized.
3. Coverage `start` and `end` are equal point coordinates, not BED intervals.
4. The summary mapped-read metric can include secondary or supplementary records
   if they meet its mapped, MAPQ, and `NH` conditions.
5. Saturation may be reported as unavailable when the required calculation is
   not implemented for an input category.

These constraints should be considered when extending the workflow or comparing
its results with external methylation pipelines.
