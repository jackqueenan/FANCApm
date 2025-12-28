# ABE VLP Off-Target Quantification

This repository provides a Python script for quantifying **adenine base editor (ABE) off-target editing** from CRISPResso2 allele frequency tables for a given point mutation type within a given editing window. The script is designed for batch analysis of multiple amplicons and samples and supports both forward- and reverse-complemented protospacers with correct editing-window handling.

The workflow is optimized for **Vbase editing experiments**.

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

## Requirements

- Python **3.10 or higher**
- `pandas`

Install dependencies with:

```bash
pip install pandas
