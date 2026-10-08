# SmC Library Structure

SmC uses the same installation and command-line interface as the other assays.
See the [technical reference](technical.md) for usage and configuration defaults.

## Barcode-bearing read

```text
[fixed leader][barcode B][linker-BC][barcode A][linker2/insert-left][genomic insert]
```

| Element | Length | Configuration / name |
|---|---:|---|
| Barcode B | 11 bp | `barcode2`; `SMC_BARCODE_WHITELIST` |
| Linker-BC | 30 bp | `SMC_LINKER_BC` |
| Barcode A | 11 bp | `barcode1`; `SMC_BARCODE_WHITELIST` |
| Insert-left | 34 bp | `SMC_INSERT_LEFT` |

```text
linker-BC:   ATCCACGTGCTTGAGCGCGCTGCATACTTG
insert-left: TGCAGTCGTGCCATGAGATGTGTATAAGAGACAG
```


## Conversion classes and mates

The anchors tolerate expected C/T and G/A conversions during location.
Classification uses the original C positions: retained C supports Watson,
and converted T supports Crick. Conflicting or insufficient anchor evidence
is ambiguous. These classes describe conversion lineage, not genomic
alignment strand.

| Class | Meaning | BISCUIT R1 | BISCUIT R2 |
|---|---|---|---|
| Watson | Valid barcodes; C-retained anchor evidence | Long genomic parent | Short barcode-derived daughter |
| Crick | Valid barcodes; C-to-T anchor evidence | Short barcode-derived daughter | Long genomic parent |
| Ambiguous | Valid barcodes; uncertain boundary or conversion class | Not aligned downstream | Not aligned downstream |
| Discarded | Structure or barcode correction failed | Not aligned downstream | Not aligned downstream |

SmC alignment uses directional mode `1`; mate assignment changes without
reverse-complementing the sequences. Watson and Crick BAMs stay separate;
called counts are merged for downstream coverage analysis.

## Parameter meanings

| Parameter | Meaning |
|---|---|
| `SMC_BARCODE_WHITELIST` | Barcode list; empty selects the bundled 96-entry, 11 bp whitelist. |
| `SMC_LINKER_BC` | Anchor between the two barcodes. |
| `SMC_INSERT_LEFT` | Anchor defining the genomic insert boundary. |
| `BARCODE_LINKER_EDIT_DISTANCE` | Allowed anchor errors outside expected conversion states. |
| `BARCODE_HAMMING_DISTANCE` | Maximum substitutions for unique nearest-whitelist correction. |
| `SMC_MAX_CONVERSION_MISMATCHES` | Maximum C/T-position mismatches for a conversion class. |
| `SMC_MINIMUM_SCORE_MARGIN` | Required difference between C-retained and C-to-T evidence. |
| `SMC_SAVE_UNINFORMATIVE_FASTQ` | Save ambiguous/discarded FASTQs when `1`; their counts are always recorded. |
| `SMC_FILTER_THREADS` | Worker processes for the post-M-bias XXX non-conversion filter. |

The XXX filter removes alignments containing three successive unconverted
cytosine observations after applying the matching M-bias trims. Host filtering
uses `CALL_CHROMOSOMES`; spike-in filtering covers all contigs. The lambda-rich
control spot `00_01` is excluded from saturation fitting and per-spot summaries.
