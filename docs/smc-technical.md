# SmC Barcode and Linker Handling

This document describes only the barcode-bearing read structure and the rules
used by `script/steps/python/02.smc_extract_bc.py` to locate linkers, extract
spatial barcodes, and assign Watson or Crick conversion classes.

## Barcode-bearing read structure

The relevant part of the read has the following layout:

```text
[fixed leader][barcode B][linker-BC][barcode A][linker2/insert-left][genomic insert]
```

The extractor uses the following names:

| Library element | Length | Extractor name | Position |
|---|---:|---|---|
| Barcode B | 11 bp | `barcode2` | Immediately before linker-BC |
| Linker-BC | 30 bp | `linker_bc` | Between the two barcodes |
| Barcode A | 11 bp | `barcode1` | Immediately after linker-BC |
| Insert-left | 34 bp | `insert_left` | Immediately before the genomic insert |

The fixed sequences are:

```text
linker-BC:   ATCCACGTGCTTGAGCGCGCTGCATACTTG
insert-left: TGCAGTCGTGCCATGAGATGTGTATAAGAGACAG
```

The order written to FASTQ headers is `barcode1+barcode2`, or barcode A
followed by barcode B. This is intentionally the reverse of their physical
left-to-right order around linker-BC:

```text
physical read:  barcode2 | linker-BC | barcode1
FASTQ header:   barcode1+barcode2
```

## Why linker2 is not searched separately

The library annotation contains a 30 bp barcode A linker, referred to here as
linker2:

```text
linker2:     CCCATGATCGTCCGA TGCAGTCGTGCCATG
                               |||||||||||||||
insert-left:                  TGCAGTCGTGCCATG AGATGTGTATAAGAGACAG
```

The last 15 bp of linker2 are also the first 15 bp of insert-left. Therefore,
linker2 and insert-left are overlapping annotations rather than two
consecutive adapters. From the end of barcode1 to the start of the genomic
insert, the structure is:

```text
15 bp unique to linker2 + 34 bp insert-left = 49 bp
```

Searching for linker2 and then searching for insert-left as independent,
non-overlapping elements would place the genomic boundary incorrectly. The
extractor instead uses linker-BC to delimit the barcodes and the right edge of
insert-left to define the start of the short genomic sequence.

## Conversion-compatible anchor location

SmC reads contain different conversion states, so an exact search for the
unconverted linker sequence would preferentially find only one class. Both
linker-BC and insert-left are located with a conversion-compatible mask:

```text
Expected C: observed C or T
Expected G: observed G or A
Expected A/T: exact match
```

`--anchor-edit-distance` controls the number of additional mismatches allowed
outside these expected conversion states. A match is accepted only when the
best position is unique.

Anchor location and strand classification are separate operations:

1. The conversion-compatible mask locates the anchor and determines its
   boundaries.
2. The original C positions within the located anchor provide C-retained or
   C-to-T evidence for Watson/Crick classification.

The reference G positions are tolerated during anchor location but are not
used in the C/T strand score.

## Barcode extraction and correction

For each read pair, the extractor performs these steps:

1. Locate linker-BC near the beginning of the barcode-bearing read.
2. Extract the 11 bases immediately before it as `barcode2`.
3. Extract the 11 bases immediately after it as `barcode1`.
4. Match both sequences against `--barcode-whitelist`.
5. Search for insert-left after barcode1.
6. Use the right edge of insert-left as the genomic insert boundary.
7. Classify the read using the C/T evidence in linker-BC and insert-left.

A barcode is accepted when it is an exact whitelist entry or has a unique
nearest whitelist entry within `--barcode-hamming-distance`. The corrected
whitelist sequences, rather than the observed error-containing sequences, are
written to the read header. A tied nearest match or a distance above the limit
is rejected.

## Watson and Crick classification

Linker-BC contains nine reference C positions and insert-left contains five.
Each anchor is scored independently:

```text
Watson mismatches = observed T + bases other than C/T
Crick mismatches  = observed C + bases other than C/T
```

A C-retained anchor supports Watson, while a C-to-T anchor supports Crick.
The maximum mismatch count is set by `--max-conversion-mismatches`, and the
required difference between C and T evidence is set by
`--minimum-score-margin`.

The two anchor calls are combined as follows:

| Linker-BC | Insert-left | Final class |
|---|---|---|
| Watson | Watson | Watson |
| Crick | Crick | Crick |
| Watson or Crick | Ambiguous | Use the informative anchor |
| Ambiguous | Watson or Crick | Use the informative anchor |
| Watson | Crick | Ambiguous |
| Crick | Watson | Ambiguous |
| Ambiguous | Ambiguous | Ambiguous |

These labels describe the linker conversion class determined before genomic
alignment. They are not inferred from the strand to which the read later
aligns.

## Output classes

Read pairs are assigned to four groups:

| Group | Meaning |
|---|---|
| `watson` | Both barcodes are valid and the conversion class is Watson |
| `crick` | Both barcodes are valid and the conversion class is Crick |
| `ambiguous` | Barcodes are valid, but the insert boundary or conversion class is uncertain |
| `discarded` | Linker/barcode extraction failed or a barcode could not be corrected uniquely |

When both barcodes are available, the read name is annotated as:

```text
@barcode1+barcode2:W|C|U|D:original_read_name
```

`W`, `C`, `U`, and `D` mean Watson, Crick, ambiguous, and discarded,
respectively. If linker-BC cannot be located, the barcodes are unknown and the
original header is retained.

The `discarded` group is a technical rejection category. It is not a
biological spike-in category and must not be used to identify lambda reads.
Barcoded lambda molecules pass through the same linker, barcode, and
Watson/Crick classification rules as other molecules.
