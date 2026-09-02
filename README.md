# DBiT-spatial-Methylation

DBiT-spatial-Methylation processes paired-end DBiT methylation sequencing data from TAPS, TAPS-v2, EM-seq, Cabernet or SmC experiments. The command-line workflow supports local execution and SLURM submission.

## Input organization

Place exactly one R1/R2 FASTQ pair in an input directory. The filenames must
contain an independent `R1` or `R2` token.

```text
sample/
└── fastq/
    ├── sample_R1.fastq.gz
    └── sample_R2.fastq.gz
```

Pipeline results are written to `sample/dbitm/`, next to the input FASTQ
directory.

## Installation

Install [Pixi](https://pixi.prefix.dev/latest/installation/) first. Then clone
the repository and install the locked environment and the `dbitm` command:

```bash
cd /path/to/DBiT-spatial-Methylation
pixi run init
source ~/.bashrc
dbitm --help
```

By default, the command is installed at `~/.local/bin/dbitm`. Set
`DBITM_INSTALL_DIR` before installation to use another directory.

## Reference preparation

Build the reference indexes required by the selected assay.

For TAPS and TAPS-v2:

```bash
pixi run -e default bwa index /path/to/genome.fa
pixi run -e default samtools faidx /path/to/genome.fa
```

For EM-seq, Cabernet or SmC:

```bash
pixi run -e default biscuit index /path/to/genome.fa
pixi run -e default samtools faidx /path/to/genome.fa
```

Prepare lambda and pUC19 FASTAs in the same way when spike-in alignment is required.

## Configuration

Create the local configuration:

```bash
cp config/dbitm.config.example.sh config/dbitm.config.sh
```

Edit `config/dbitm.config.sh` and set at least the following values.

For TAPS or TAPS-v2:

```bash
BWA_INDEX=/path/to/genome.fa
CALL_REFERENCE=/path/to/genome.fa

declare -A BWA_SPIKE_IN_INDEXES=(
    [lambda]="/path/to/lambda.fa"
    [puc19]="/path/to/puc19.fa"
)
```

For EM-seq, Cabernet or SmC:

```bash
BISCUIT_REFERENCE=/path/to/genome.fa
CALL_REFERENCE=/path/to/genome.fa

declare -A BISCUIT_SPIKE_IN_INDEXES=(
    [lambda]="/path/to/lambda.fa"
    [puc19]="/path/to/puc19.fa"
)
```

Other commonly changed settings are:

```bash
RUN_MODE=hpc                 # hpc or local
SCRATCH_ROOT=/path/to/scratch
BARCODE_WHITELIST=           # empty uses docs/barcodes/barcodes50.tsv
SMC_BARCODE_WHITELIST=       # empty uses docs/barcodes/barcodes-smc.tsv
CALL_CONTEXT_MODE=both       # cg, ch, or both
```

Set `SCRATCH_ROOT` to a directory on a fast SSD. Adjust the SLURM CPU, memory,
time, and partition settings near the end of the configuration file.

## Usage

The command format is:

```text
dbitm <assay> <step> --input <FASTQ_directory> [--config <config_file>] [--resume <step>] [--dry-run]
```

Supported assays:

```text
taps
taps-v2
emseq
cabernet
smc
```

Run the complete workflow:

```bash
dbitm taps all --input /path/to/sample/fastq
```

With `RUN_MODE=local`, the steps run sequentially in the current process. With
`RUN_MODE=hpc`, they are submitted to SLURM with the required dependencies.

Resume at any step and run or submit that step plus all later steps with
`step=all` and `--resume`. Outputs required by earlier steps must already be
complete. For example:

```bash
dbitm taps all --resume mbias --input /path/to/sample/fastq
```


Use `--dry-run` to validate configuration and print the execution plan without
writing pipeline outputs. In HPC mode, this prints the `sbatch` commands and
dependencies without submitting jobs:

```bash
dbitm taps all --input /path/to/sample/fastq --dry-run
```

The required order is:

```text
fastp -> barcode -> (spike-align + align) -> pool -> mbias -+-> spike-call ---------+
                                                           +-> call -> saturation -+-> summary -> methscan
```

Use another configuration file with:

```bash
dbitm taps all \
    --input /path/to/sample/fastq \
    --config /path/to/dbitm.config.sh
```

## Output organization

```text
sample/
├── fastq/
└── dbitm/
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

Run `dbitm --help` to view the current command options.
