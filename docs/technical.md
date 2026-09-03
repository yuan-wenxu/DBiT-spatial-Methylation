# DBiT-spatial-Methylation Technical Reference

## 1. Overview

DBiT-spatial-Methylation processes paired-end spatial DNA methylation data. It
filters FASTQs, extracts spatial barcodes, aligns host and spike-in reads,
estimates M-bias, calls methylation, and generates summary and MethSCAn outputs.

Supported assay names are `taps`, `taps-v2`, `emseq`, `cabernet`, and `smc`.
SmC library structure and strand splitting are documented separately in
`docs/smc-technical.md`.

```text
fastp -> barcode -> +-> spike-align -+
                    |                 |
                    +-> align --------+-> pool -> mbias -+-> spike-call ---------+
                                                         |                       |
                                                         +-> call -> saturation -+-> summary -> methscan
```

`spike-align` and `align` can run concurrently. After M-bias, `spike-call` and
`call` can also run concurrently.

| Stage | Entry point | Purpose |
|---|---|---|
| `fastp` | `script/steps/01.fastp.sh` | Filter paired FASTQs |
| `barcode` | `script/steps/02.barcode.sh` | Extract and correct spatial barcodes |
| `spike-align` | `script/steps/03.spike_align.sh` | Align spike-in candidates |
| `align` | `script/steps/03.align.sh` | Align host reads and add `CB` tags |
| `pool` | `script/steps/04.pool.sh` | Pool, sort, and index BAMs |
| `mbias` | `script/steps/05.mbias.sh` | Infer cycle trimming cutoffs |
| `spike-call` | `script/steps/06.spike_call.sh` | Call mitochondrial and spike-in methylation |
| `call` | `script/steps/06.call.sh` | Call per-spot host methylation |
| `saturation` | `script/steps/07.saturation.sh` | Estimate CpG saturation |
| `summary` | `script/steps/08.summary.sh` | Generate QC tables and plots |
| `methscan` | `script/steps/09.methscan.sh` | Generate context-specific VMR matrices |

The source code is authoritative if it differs from this document.

## 2. Installation and use

The project uses Pixi. Install the command-line entry point with:

```bash
./install-cli.sh
```

Run one stage or the complete pipeline with:

```bash
dbitm <assay> <step> --input /path/to/sample/fastq
dbitm taps all --input /data/sample/fastq
```

Useful options are:

- `--config PATH`: select a configuration file.
- `--resume STEP`: start an `all` run from a particular stage.
- `--dry-run`: print or validate the execution plan without running tools.

For input `/data/sample/fastq`, results are written under
`/data/sample/dbitm`. `RUN_MODE=local` executes stages directly;
`RUN_MODE=hpc` submits Slurm jobs with stage dependencies.

Copy the example configuration before a production run:

```bash
cp config/dbitm.config.example.sh config/dbitm.config.sh
```

The main reference settings are:

| Purpose | Configuration |
|---|---|
| TAPS host alignment | `BWA_INDEX` |
| TAPS spike-in alignment | `BWA_SPIKE_IN_INDEXES` |
| EM-seq/Cabernet/SmC host alignment | `BISCUIT_REFERENCE` |
| EM-seq/Cabernet/SmC spike-in alignment | `BISCUIT_SPIKE_IN_INDEXES` |
| Host calling | `CALL_REFERENCE` |
| Spike-in calling | `CALL_SPIKE_IN_REFERENCES` |

Calling FASTAs require SAMtools `.fai` indexes. Alignment references require
the indexes expected by BWA or BISCUIT.

Frequently adjusted processing settings include:

- barcode whitelist, linker, insert-left, correction distance, chunks, and
  compression;
- `CALL_CONTEXT_MODE` (`cg`, `ch`, or `both`);
- mapping/base-quality thresholds and caller parallelism;
- M-bias sampling and stability thresholds;
- Slurm CPU, memory, time, and partition values.

See `config/dbitm.config.example.sh` for all settings.

## 3. Data conventions

### 3.1 Spatial barcode

For the original assays, the barcode-bearing mate has the conceptual layout:

```text
[barcode2][linker][barcode1][insert-left][biological insert]
```

The barcode stage locates this structure, corrects both barcodes against the
whitelist, removes the non-biological prefix, and annotates the read name with:

```text
barcode1+barcode2
```

Host alignment moves this value into the BAM `CB` tag. Structure or barcode
failures are written to spike-in candidate FASTQs.

SmC uses different anchors and additionally splits informative reads into
Watson and Crick groups. See `docs/smc-technical.md` for the exact layout.

### 3.2 Coordinates and read cycles

Coverage output uses zero-based genomic coordinates. The `start` and `end`
fields both contain the called cytosine-representation coordinate; these files
are not BED half-open intervals.

CpG evidence from both reference strands is merged at the forward-strand C. CH
output retains strand: plus sites use the reference C coordinate, and minus
sites use the corresponding reference G coordinate.

M-bias cycles are one-based and measured from the original molecule's 5-prime
end. Reverse alignments invert the query position. R1 also receives the offset
for sequence removed during barcode extraction; R2 does not.

### 3.3 Methylation evidence

| Assay | Converted `TG`/`CA` | Retained `CG` |
|---|---|---|
| TAPS / TAPS-v2 | methylated | unmethylated |
| EM-seq / Cabernet / SmC | unmethylated | methylated |

Only paired-end orientation flags `83`, `99`, `147`, and `163` are used for
methylation classification.

For CpG, plus- and minus-family evidence updates one forward-C counter. For CH
(`CA`, `CC`, and `CT`):

- forward-aligned reads (`99`, `163`) contribute to plus-strand C sites;
- reverse-aligned reads (`83`, `147`) contribute to minus-strand G sites;
- plus evidence is observed C/T and minus evidence is observed G/A;
- the neighboring reference/read base must support the requested context.

CH sites on the two strands are therefore not merged.

## 4. Pipeline stages

### 4.1 FASTQ filtering

`01.fastp.sh` identifies one R1/R2 pair and runs fastp. Adapter trimming is
disabled because barcode structure is handled by the next stage.

```text
dbitm/fastp/
├── R1.filtered.fastq.gz
├── R2.filtered.fastq.gz
├── fastp.json
├── fastp.html
└── fastp.log
```

### 4.2 Barcode extraction

`02.barcode.sh` dispatches the assay-specific Python extractor. Barcode
correction uses the configured whitelist and Hamming-distance limit.

Standard output chunks contain:

```text
0001.R1.demux.fastq.gz
0001.R2.demux.fastq.gz
0001.R1.spike-in.fastq.gz
0001.R2.spike-in.fastq.gz
0001.stats.json
stats.json
barcode.log
```

SmC writes Watson, Crick, ambiguous, and discarded pairs. Its stats follow:

```text
total_reads = kept_reads + discarded_reads
kept_reads = informative_reads + ambiguous_reads
informative_reads = watson_reads + crick_reads
```

`kept_reads` means barcode extraction succeeded; `informative_reads` can enter
strand-specific SmC alignment.

### 4.3 Host and spike-in alignment

`03.align.sh` aligns demultiplexed host chunks:

- TAPS/TAPS-v2: `bwa mem`.
- EM-seq/Cabernet: `biscuit align` using `BISCUIT_DIRECTIONAL_MODE`.
- SmC: Watson and Crick are aligned separately with `biscuit align -b 1`.

For SmC the mate assignments are:

| Group | BISCUIT R1/parent | BISCUIT R2/daughter |
|---|---|---|
| Watson | `NNNN.watson.genomic.fastq.gz` | `NNNN.watson.short-genomic.fastq.gz` |
| Crick | `NNNN.crick.short-genomic.fastq.gz` | `NNNN.crick.genomic.fastq.gz` |

No sequence is reverse-complemented; only the R1/R2 arguments are reordered.
The alignment stream passes through `sinto nametotag` to create the `CB` tag.

Standard assays produce `NNNN.cb.bam`; SmC produces
`NNNN.watson.cb.bam` and `NNNN.crick.cb.bam`. Logs are stored in
`dbitm/align/logs/`. BAMs are checked but are not sorted or indexed until pool.

`03.spike_align.sh` maps each spike-in FASTQ pair independently to every
configured spike-in reference and writes one BAM and flagstat report per
chunk/reference combination. For SmC, spike-in molecules are already present
in the Watson and Crick FASTQs, so both groups are aligned with the same mate
assignments shown above and `biscuit align -b 1`. SmC outputs retain the group
in their names, for example `0001.watson.lambda.bam` and
`0001.crick.lambda.bam`; discarded reads are not used as spike-in input.

### 4.4 BAM pooling

`04.pool.sh` concatenates host chunk BAMs, coordinate-sorts the result, and
creates its index. Spike-in BAMs are pooled independently.

```text
dbitm/pooled/
├── pooled.cb.bam
├── pooled.cb.bam.bai
├── pooled.<spike>.bam
├── pooled.<spike>.bam.bai
└── pool.log
```

`POOL_SORT_MEM` is memory per SAMtools sort thread. Spike-in names are read
from the assay-specific configured alignment-index array.

### 4.5 M-bias

`05.mbias.sh` measures CpG methylation by R1/R2 cycle and infers end-trimming
cutoffs. It excludes unmapped, secondary, supplementary, unsupported-orientation,
and low-quality evidence. Host reads may be deterministically subsampled;
spike-in targets are processed in full.

Each target produces:

```text
<label>.mbias.tsv
<label>.mbias.png
<label>.mbias.cutoffs.tsv
```

If no stable region can be inferred, a zero-byte cutoff marker means subsequent
calling proceeds without trimming for that target.

### 4.6 Methylation calling

`06.call.sh` creates per-spot host calls from the pooled BAM, reference FASTA,
barcode whitelist, and host M-bias cutoff. All whitelist spots are recorded in
`coverage/spot_manifest.tsv`, including empty spots.

```text
coverage/host/<X>/<X>_<Y>.CG.cov
coverage/host/<X>/<X>_<Y>.CA.cov
coverage/host/<X>/<X>_<Y>.CC.cov
coverage/host/<X>/<X>_<Y>.CT.cov
```

CG output contains genomic position, methylation percentage, methylated count,
and unmethylated count. CH output additionally records context and strand.

`06.spike_call.sh` creates aggregate, non-spatial mitochondrial and spike-in
calls. `SPIKE_CALL_MODE` selects `mito`, `spike`, or `all`.

### 4.7 CpG saturation

`07.saturation.sh` selects sufficiently covered spots, estimates unique CpGs at
subsampling fractions, and fits linear and saturation-curve models.

```text
dbitm/saturation/
├── saturation_curve.png
├── saturation_summary.tsv
└── saturation.log
```

The summary records the configured multiplier in `prediction_fraction` and
the corresponding estimate in `predicted_median_unique_cpgs`.

### 4.8 Summary

`08.summary.sh` combines fastp, barcode, BAM, coverage, spike-in, and saturation
statistics. It writes per-spot and sample-level tables plus context-specific
heatmaps and violin plots.

Per-spot methylation gives every called site equal weight:

```text
spot mean = sum(per-site methylation rate) / number of called sites
```

It is not calculated by pooling methylated and unmethylated read counts across
sites. Host-wide context means weight each spot mean by its called-site count,
which is equivalent to equal weighting of spot-by-site rows.

Primary outputs are:

```text
dbitm/summary/
├── per_spot_summary.tsv
├── sample_summary.tsv
├── reads_heatmap.png
├── <context>_site_count_heatmap.png
├── mean_<context>_methylation_heatmap.png
├── mean_<context>_methylation_violin.png
└── summary.log
```

Only plots for contexts selected by `CALL_CONTEXT_MODE` are generated.

### 4.9 MethSCAn

`09.methscan.sh` runs the following workflow independently for each selected
context:

```text
prepare -> filter -> smooth -> scan -> matrix
```

Outputs are stored under `dbitm/methscan/CG`, `CA`, `CC`, and `CT`. Each selected
directory contains VMRs and spot-by-VMR matrices. Missing matrix values mean no
usable observation, not zero methylation.

## 5. Output and scratch behavior

A complete run uses this top-level layout:

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

Exact files depend on assay and configuration.

When `SCRATCH_ROOT` is set, only `barcode`, `pool`, `call`, and `methscan`
create isolated scratch workspaces and copy their results back to the sample
directory. These stages have substantial intermediate or random I/O. `fastp`,
`align`, `spike-align`, `mbias`, `spike-call`, `saturation`, and `summary`
always use the persistent sample directory so long-running or completed
outputs are not held only on ephemeral storage. Handled failures in the four
scratch-enabled stages attempt to recover partial output before removing their
scratch workspace; recovered output may be incomplete and should be inspected
before rerunning.
