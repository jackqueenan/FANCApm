# Sum ABEd reads from CRISPResso2 output

This repository provides a Python script for quantifying **total adenine base edited (ABE) reads** from CRISPResso2 allele frequency tables for a given point mutation type within a given editing window. The script is designed for batch analysis of multiple amplicons and samples and supports both forward- and reverse-complemented protospacers with correct editing-window handling.

The workflow is intended for **base editing experiments** 

---

## Features

- Batch processing of multiple amplicons and samples
- Correct handling of reverse-complemented protospacers (`RC = TRUE`)
- Automatic conversion of editing windows to reverse-complement coordinates
- Orientation-agnostic matching of allele sequences
- Avoids double counting of allele rows
- Produces per-amplicon and combined summary CSV files
- Compatible with Python ≥ 3.10

---

## Inputs
- ***batch-file***: .csv/.tsv batch file containing the following variables 
  - required: ***amplicon***: Name of the sequencing amplicon (e.g. ***"FANCA"***)
  - required: ***input_seq***: Protospacer sequence 5'-3' exluding PAM (e.g. "GATTATATGTGCTAATGCTC")
  - optional: ***window*** Tuple representing the base editing window over which to calculate base edited reads (syntax is "2-10", "2:10", "2..10", or "2,10"). Alternatively, use window_start + window_end integer variables
  - optional: ***RC***: Boolean variable indicating the protospacer's sense relative to the sequencing amplicon (TRUE indicates the ***input_seq*** is reverse complement of the HTS amplicon read1
  - optional: ***print-variants***: print all base editing variants over which to sum read counts? 
- ***base-crispresso-dir***: directory containing CRISPResso_batch_on_<amplicon> folders

## Outputs
- Per-amplicon CSV in each CRISPResso_batch_on_<amplicon> directory: "ABE_reads__<amplicon>__<input_seq_sanitized>.csv"
- Compiled CSV in --base-crispresso-dir: "ABE_reads__combined.csv"

## Requirements

- Python **3.10 or higher**
- `pandas`

Install dependencies with:

```bash
pip install pandas
