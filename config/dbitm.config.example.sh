# execution mode: local (run directly) or hpc (submit via sbatch)
RUN_MODE=${RUN_MODE:-local}

# reference configuration
SCRATCH_ROOT=${SCRATCH_ROOT:-}
# TAPS/TAPS-v2 bwa index prefix or indexed reference FASTA
BWA_INDEX=${BWA_INDEX:-}
# Absolute paths to spike-in BWA indexes used by TAPS/TAPS-v2 align.
# Bash equivalent of: {"lambda": "/path/to/index", "puc19": "/path/to/index"}
declare -A BWA_SPIKE_IN_INDEXES=(
    [lambda]=""
    [puc19]=""
)
# EM-seq biscuit-indexed reference FASTA
BISCUIT_REFERENCE=${BISCUIT_REFERENCE:-}
# Absolute paths to spike-in biscuit indexes used by EM-seq align.
# Bash equivalent of: {"lambda": "/path/to/index", "puc19": "/path/to/index"}
declare -A BISCUIT_SPIKE_IN_INDEXES=(
    [lambda]=""
    [puc19]=""
)
# Host reference FASTA used for methylation calling.
CALL_REFERENCE=${CALL_REFERENCE:-}
# Absolute paths to spike-in FASTAs used for methylation calling; each FASTA needs a .fai index.
# Bash equivalent of: {"lambda": "/path/to/FASTA", "puc19": "/path/to/FASTA"}
declare -A CALL_SPIKE_IN_REFERENCES=(
    [lambda]=""
    [puc19]=""
)

# barcode extraction configuration
# set BARCODE_WHITELIST to a custom whitelist file path if needed
# empty BARCODE_WHITELIST will use the default whitelist in docs/barcodes/barcodes50.tsv
BARCODE_WHITELIST=${BARCODE_WHITELIST:-}
# batches are assigned to these chunks in round-robin order
BARCODE_CHUNK=${BARCODE_CHUNK:-4}
# paired reads in each dispatched batch
BARCODE_BATCH_SIZE=${BARCODE_BATCH_SIZE:-50000}
# compression step for barcode extraction: shell (gzip) or python (zlib)
BARCODE_COMPRESSION_STEP=${BARCODE_COMPRESSION_STEP:-shell}
BARCODE_GZIP_LEVEL=${BARCODE_GZIP_LEVEL:-6}
BARCODE_LINKER_BC=${BARCODE_LINKER_BC:-ATCCACGTGCTTGAGAGGCCAGAGCATTCG}
BARCODE_INSERT_LEFT=${BARCODE_INSERT_LEFT:-CATCGGCGTACGACTAGATGTGTATAAGAGACAG}
# methylated C positions in BARCODE_INSERT_LEFT, using 0-based coordinates
BARCODE_METHYLATED_C_POSITIONS=${BARCODE_METHYLATED_C_POSITIONS:-3,6,10}
BARCODE_LINKER_EDIT_DISTANCE=${BARCODE_LINKER_EDIT_DISTANCE:-1}
BARCODE_HAMMING_DISTANCE=${BARCODE_HAMMING_DISTANCE:-1}
BARCODE_INSERT_LEFT_EDIT_DISTANCE=${BARCODE_INSERT_LEFT_EDIT_DISTANCE:-1}
BARCODE_PROGRESS_READS=${BARCODE_PROGRESS_READS:-1000000}

# genome alignment configuration
# aligner threads used by each barcode chunk
ALIGN_THREADS_PER_CHUNK=${ALIGN_THREADS_PER_CHUNK:-8}
# aligner threads used by each spike-in FASTQ chunk
SPIKE_ALIGN_THREADS_PER_CHUNK=${SPIKE_ALIGN_THREADS_PER_CHUNK:-8}
# paired-end library mode: 1 for directional, 0 for non-directional
BISCUIT_DIRECTIONAL_MODE=${BISCUIT_DIRECTIONAL_MODE:-1}

# BAM pooling configuration
# per-thread memory for samtools sort (defaults to POOL_MEM/POOL_THREADS if empty)
POOL_SORT_MEM=${POOL_SORT_MEM:-}

# methylation calling configuration
CALL_CHROMOSOMES=${CALL_CHROMOSOMES:-chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX}
CALL_MITO_CHROMOSOMES=${CALL_MITO_CHROMOSOMES:-chrM}
CALL_CONTEXT_MODE=${CALL_CONTEXT_MODE:-both}
CALL_MIN_BASE_QUALITY=${CALL_MIN_BASE_QUALITY:-30}
CALL_MIN_MAPPING_QUALITY=${CALL_MIN_MAPPING_QUALITY:-10}
CALL_MAX_DEPTH=${CALL_MAX_DEPTH:-250}
CALL_BATCH_SIZE=${CALL_BATCH_SIZE:-10000000}
CALL_R1_LEFT_TRIMMING=${CALL_R1_LEFT_TRIMMING:-0}
CALL_R1_RIGHT_TRIMMING=${CALL_R1_RIGHT_TRIMMING:-0}
CALL_R2_LEFT_TRIMMING=${CALL_R2_LEFT_TRIMMING:-0}
CALL_R2_RIGHT_TRIMMING=${CALL_R2_RIGHT_TRIMMING:-0}
CALL_JOBS=${CALL_JOBS:-8}

# SLURM resource configuration
SBATCH_OUTPUT=${SBATCH_OUTPUT:-%x_%j.out}
SBATCH_ERROR=${SBATCH_ERROR:-%x_%j.err}
SBATCH_REQUEUE=${SBATCH_REQUEUE:-true}

FASTP_NAME=${FASTP_NAME:-fastp}
FASTP_THREADS=${FASTP_THREADS:-8}
FASTP_PARTITION=${FASTP_PARTITION:-}
FASTP_MEM=${FASTP_MEM:-16G}
FASTP_TIME=${FASTP_TIME:-24:00:00}

BARCODE_NAME=${BARCODE_NAME:-barcode}
BARCODE_THREADS=${BARCODE_THREADS:-$BARCODE_CHUNK}
BARCODE_PARTITION=${BARCODE_PARTITION:-}
BARCODE_MEM=${BARCODE_MEM:-64G}
BARCODE_TIME=${BARCODE_TIME:-24:00:00}

ALIGN_NAME=${ALIGN_NAME:-align}
ALIGN_THREADS=${ALIGN_THREADS:-$((BARCODE_CHUNK * ALIGN_THREADS_PER_CHUNK))}
ALIGN_PARTITION=${ALIGN_PARTITION:-}
ALIGN_MEM=${ALIGN_MEM:-128G}
ALIGN_TIME=${ALIGN_TIME:-24:00:00}

SPIKE_ALIGN_NAME=${SPIKE_ALIGN_NAME:-spike-align}
SPIKE_ALIGN_THREADS=${SPIKE_ALIGN_THREADS:-$((BARCODE_CHUNK * SPIKE_ALIGN_THREADS_PER_CHUNK))}
SPIKE_ALIGN_PARTITION=${SPIKE_ALIGN_PARTITION:-}
SPIKE_ALIGN_MEM=${SPIKE_ALIGN_MEM:-64G}
SPIKE_ALIGN_TIME=${SPIKE_ALIGN_TIME:-24:00:00}

POOL_NAME=${POOL_NAME:-pool}
POOL_THREADS=${POOL_THREADS:-4}
POOL_PARTITION=${POOL_PARTITION:-}
POOL_MEM=${POOL_MEM:-64G}
POOL_TIME=${POOL_TIME:-24:00:00}

CALL_NAME=${CALL_NAME:-call}
CALL_THREADS=${CALL_THREADS:-$CALL_JOBS}
CALL_PARTITION=${CALL_PARTITION:-}
CALL_MEM=${CALL_MEM:-256G}
CALL_TIME=${CALL_TIME:-48:00:00}
