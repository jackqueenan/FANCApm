#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python3
"""
abe_offtarget_vlp.py

Quantify ABE off-target editing from CRISPResso2 Alleles_frequency_table_around_sgRNA_*.txt
files using a batch file specifying amplicon name, protospacer, editing window, and whether
the protospacer is reverse-complemented (RC).

This version:
- Takes REQUIRED command-line inputs for batch file and CRISPResso base directory.
- Is compatible with Python 3.7+ (uses typing.Tuple/List/Set instead of tuple[int,int]).
- Handles RC window conversion correctly.
- Avoids double-counting allele rows when multiple variants match.
- Optionally prints variants in 3 aligned columns.

Inputs
- --batch-file: CSV/TSV with columns:
    required: amplicon, input_seq
    optional: window (e.g. "2-10", "2:10", "2..10", "2,10") OR window_start + window_end
    optional: RC (truthy/falsey values)
- --base-crispresso-dir: directory containing CRISPResso_batch_on_<amplicon> folders

Outputs
- Per-amplicon CSV in each CRISPResso_batch_on_<amplicon> directory:
    ABE_reads__<amplicon>__<input_seq_sanitized>.csv
- Combined CSV in --base-crispresso-dir:
    ABE_reads__combined.csv
"""

from __future__ import print_function

import argparse
import os
import re
import sys
from typing import Iterable, List, Set, Tuple, Optional

import pandas as pd


# --------------------
# Constants
# --------------------
ALLELE_TABLE_PREFIX = "Alleles_frequency_table_around_sgRNA_"  # file is {prefix}{sgRNA}.txt

_RC_MAP = str.maketrans({
    "A": "T", "T": "A", "C": "G", "G": "C",
    "N": "N", "n": "n",
    "-": "-",
})


# --------------------
# CLI
# --------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="abe_offtarget_vlp.py",
        description="Quantify ABE off-target editing from CRISPResso allele tables using a batch file."
    )
    p.add_argument(
        "--batch-file",
        required=True,
        help="CSV/TSV batch file with columns: amplicon,input_seq,(window or window_start+window_end),RC"
    )
    p.add_argument(
        "--base-crispresso-dir",
        required=True,
        help="Directory containing CRISPResso_batch_on_<amplicon> folders"
    )
    p.add_argument(
        "--variant-regex-chunk",
        type=int,
        default=500,
        help="Variants per regex chunk to avoid huge regexes (default: 500)"
    )
    p.add_argument(
        "--print-variants",
        action="store_true",
        help="Print generated variant strings (3 columns). Disable for large runs."
    )
    p.add_argument(
        "--print-orientation-equivalents",
        action="store_true",
        help="When RC=TRUE, print the orientation-equivalent variants (reverse-complements) instead of biological variants."
    )
    return p.parse_args()


# --------------------
# Utilities
# --------------------
def reverse_complement(seq: str) -> str:
    seq = (seq or "").strip()
    return seq.translate(_RC_MAP)[::-1]


def rc_window(window: Tuple[int, int], seq_len: int) -> Tuple[int, int]:
    """
    Convert a 1-based inclusive window on a sequence to the corresponding window
    on its reverse complement.

    If window = (w_start, w_end) on forward sequence of length L, then on RC:
      rc_start = L - w_end + 1
      rc_end   = L - w_start + 1
    """
    w_start, w_end = window
    return (seq_len - w_end + 1, seq_len - w_start + 1)


def sanitize_for_filename(s: str) -> str:
    return ''.join(ch if (ch.isalnum() or ch in ('-', '_')) else '_' for ch in (s or ""))


def read_batch_file(path: str) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.tsv', '.tab'):
        return pd.read_csv(path, sep='\t')
    return pd.read_csv(path)


def parse_window(row) -> Tuple[int, int]:
    """
    Accepts:
      - window column: '2-10' or '2,10' or '2..10' or '2:10'
      - or columns window_start/window_end
    Returns (start, end) 1-based inclusive tuple.
    """
    if ('window_start' in row) and ('window_end' in row) and pd.notna(row['window_start']) and pd.notna(row['window_end']):
        return int(row['window_start']), int(row['window_end'])

    w = str(row.get('window', '')).strip()
    if not w:
        raise ValueError("Missing window info (need 'window' or 'window_start'/'window_end').")

    for sep in ('-', ',', '..', ':'):
        if sep in w:
            a, b = w.split(sep, 1)
            return int(a.strip()), int(b.strip())

    raise ValueError("Unrecognized window format: {}".format(w))


def parse_rc(row) -> bool:
    """
    Parse per-row RC setting. Accepts truthy strings/numbers.
    Defaults to False if missing/empty.
    """
    val = row.get('RC', '')
    if pd.isna(val):
        return False
    s = str(val).strip().lower()
    if s in ('1', 'true', 't', 'yes', 'y'):
        return True
    if s in ('0', 'false', 'f', 'no', 'n', ''):
        return False
    try:
        return bool(int(s))
    except Exception:
        return True


def validate_window(window: Tuple[int, int], seq_len: int, label: str) -> None:
    a, b = window
    if a < 1 or b < 1:
        raise ValueError("{} window must be 1-based positive. Got {}.".format(label, window))
    if a > b:
        raise ValueError("{} window start must be <= end. Got {}.".format(label, window))
    if a > seq_len or b > seq_len:
        raise ValueError("{} window {} out of bounds for sequence length {}.".format(label, window, seq_len))


def print_variants_in_columns(variants: Iterable[str], n_cols: int = 3, col_width: Optional[int] = None, indent: str = "  "):
    variants = sorted(set(variants))
    if not variants:
        print("{}(no variants)".format(indent))
        return

    if col_width is None:
        col_width = max(len(v) for v in variants) + 2

    for i in range(0, len(variants), n_cols):
        row = variants[i:i + n_cols]
        line = "".join(v.ljust(col_width) for v in row)
        print("{}{}".format(indent, line))


# --------------------
# Variant generation
# --------------------
def generate_variants(seq: str, window: Tuple[int, int], from_base: str, to_base: str) -> Set[str]:
    """
    Generate all sequences formed by converting one or more occurrences of `from_base` to `to_base`
    only within the 1-based inclusive window.
    """
    seq = (seq or "").strip()
    if not seq:
        return set()

    start_1b, end_1b = window
    if start_1b > end_1b:
        return set()

    start = max(0, start_1b - 1)
    end = min(len(seq) - 1, end_1b - 1)
    if start > end:
        return set()

    idxs = [i for i in range(start, end + 1) if seq[i] == from_base]
    if not idxs:
        return set()

    s = list(seq)
    out = set()

    def rec(j: int, changed: bool):
        if j == len(idxs):
            if changed:
                out.add(''.join(s))
            return
        i = idxs[j]
        orig = s[i]
        s[i] = to_base
        rec(j + 1, True)
        s[i] = orig
        rec(j + 1, changed)

    rec(0, False)
    return out


# --------------------
# Fast matching (no double counting)
# --------------------
def _build_regex_from_variants(variants: List[str]) -> str:
    return "(?:" + "|".join(re.escape(v) for v in variants) + ")"


def sum_reads_for_variants_fast(df: pd.DataFrame, variants: Set[str], variant_regex_chunk: int) -> float:
    """
    Sum %Reads for rows whose Aligned_Sequence contains ANY variant.
    Uses chunked regex OR to avoid per-variant scans and avoids double counting rows.
    """
    if not variants:
        return 0.0

    seqs = df['Aligned_Sequence'].astype(str)
    mask = None

    vlist = list(variants)
    for i in range(0, len(vlist), variant_regex_chunk):
        chunk = vlist[i:i + variant_regex_chunk]
        rx = _build_regex_from_variants(chunk)
        hit = seqs.str.contains(rx, na=False, regex=True)
        mask = hit if mask is None else (mask | hit)

    if mask is None:
        return 0.0
    return float(df.loc[mask, '%Reads'].sum())


def load_allele_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', usecols=['Aligned_Sequence', '%Reads'])
    df['%Reads'] = pd.to_numeric(df['%Reads'], errors='coerce').fillna(0.0)
    return df


# --------------------
# Per-row processing
# --------------------
def process_row(
    amplicon: str,
    input_seq: str,
    window_tuple: Tuple[int, int],
    rc_flag: bool,
    base_dir: str,
    variant_regex_chunk: int,
    print_variants: bool,
    print_orientation_equivalents: bool,
) -> Tuple[pd.DataFrame, str]:
    """
    Off-target quantification.

    RC=False:
      - generate A>G variants on input_seq within window (1-based inclusive)

    RC=True:
      - convert window_tuple to the corresponding window on rc_input
      - define biology as T>C on rc_input within that RC-adjusted window
      - to be orientation-agnostic, also include RC of those edited rc_input variants
    """
    amplicon = str(amplicon).strip()
    input_seq = str(input_seq).strip()
    if not amplicon:
        raise ValueError("Empty amplicon name.")
    if not input_seq:
        raise ValueError("Empty input_seq.")

    batch_dir = os.path.join(base_dir, "CRISPResso_batch_on_{}".format(amplicon))
    if not os.path.isdir(batch_dir):
        raise FileNotFoundError("Missing CRISPResso batch directory: {}".format(batch_dir))

    L = len(input_seq)
    validate_window(window_tuple, L, label="Input")

    if rc_flag:
        rc_input = reverse_complement(input_seq)
        rc_window_tuple = rc_window(window_tuple, L)
        validate_window(rc_window_tuple, L, label="RC")

        rc_variants = generate_variants(rc_input, rc_window_tuple, from_base='T', to_base='C')
        fwd_equiv = set(reverse_complement(v) for v in rc_variants)
        variants = rc_variants | fwd_equiv
        seq_for_filename_candidates = [rc_input, input_seq]

        print("\n[{}] RC=TRUE".format(amplicon))
        print("  protospacer:    {}".format(input_seq))
        print("  window (input): {} (inclusive)".format(window_tuple))
        print("  window (RC):    {} (inclusive)".format(rc_window_tuple))
        print("  biological variants (T>C on rc_input): {}".format(len(rc_variants)))
        print("  strings searched (incl. orientation equivalents): {}".format(len(variants)))

        if print_variants:
            print("  variant strings (3 columns):")
            to_print = fwd_equiv if print_orientation_equivalents else rc_variants
            print_variants_in_columns(to_print, n_cols=3, indent="    ")

    else:
        variants = generate_variants(input_seq, window_tuple, from_base='A', to_base='G')
        seq_for_filename_candidates = [input_seq]

        print("\n[{}] RC=FALSE".format(amplicon))
        print("  protospacer: {}".format(input_seq))
        print("  window:      {} (inclusive)".format(window_tuple))
        print("  variants (A>G on input_seq): {}".format(len(variants)))

        if print_variants:
            print("  variant strings (3 columns):")
            print_variants_in_columns(variants, n_cols=3, indent="    ")

    # Enumerate sample subfolders
    subfolders = [f.name for f in os.scandir(batch_dir) if f.is_dir()]
    agg_name = "CRISPRessoAggregate_on_{}".format(amplicon)
    if agg_name in subfolders:
        subfolders.remove(agg_name)

    per_folder_reads = []
    folder_names = []

    for folder in sorted(subfolders):
        total_reads = float('nan')
        last_err = None

        for token in seq_for_filename_candidates:
            allele_path = os.path.join(batch_dir, folder, "{}{}.txt".format(ALLELE_TABLE_PREFIX, token))
            if not os.path.exists(allele_path):
                continue
            try:
                df = load_allele_table(allele_path)
                total_reads = sum_reads_for_variants_fast(df, variants, variant_regex_chunk)
                last_err = None
                break
            except Exception as e:
                last_err = e

        if last_err is not None and pd.isna(total_reads):
            print("  Skipping {}: {}".format(folder, last_err))

        per_folder_reads.append(total_reads)
        folder_names.append(folder)
        print("  Folder: {:40s} Total %Reads: {}".format(folder, total_reads))

    out_df = pd.DataFrame({
        'File': folder_names,
        '% aligned reads containing ABE mutation within spacer': per_folder_reads
    }).sort_values('File')

    safe_in = sanitize_for_filename(input_seq)
    per_row_out = os.path.join(batch_dir, "ABE_reads__{}__{}.csv".format(amplicon, safe_in))
    out_df.to_csv(per_row_out, index=False)
    print("  Wrote per-row output: {}".format(per_row_out))

    # Transposed single-row DF for combined output
    if out_df.empty:
        row_df = pd.DataFrame(index=[amplicon])
    else:
        row_series = out_df.set_index('File')['% aligned reads containing ABE mutation within spacer']
        row_df = row_series.T.to_frame().T
        row_df.index = [amplicon]

    return row_df, per_row_out


# --------------------
# Main
# --------------------
def main() -> int:
    args = parse_args()

    batch_file = args.batch_file
    base_dir = args.base_crispresso_dir
    variant_regex_chunk = args.variant_regex_chunk

    if variant_regex_chunk < 1:
        print("ERROR: --variant-regex-chunk must be >= 1", file=sys.stderr)
        return 2
    if not os.path.exists(batch_file):
        print("ERROR: batch file not found: {}".format(batch_file), file=sys.stderr)
        return 2
    if not os.path.isdir(base_dir):
        print("ERROR: base CRISPResso dir not found: {}".format(base_dir), file=sys.stderr)
        return 2

    batch_df = read_batch_file(batch_file)

    required = {'amplicon', 'input_seq'}
    if not required.issubset(batch_df.columns):
        raise ValueError("Batch file must contain columns: {}. Optional: window OR window_start+window_end, RC."
                         .format(required))

    combined_rows = []
    for idx, row in batch_df.iterrows():
        try:
            amplicon = str(row['amplicon']).strip()
            input_seq = str(row['input_seq']).strip()
            window_tuple = parse_window(row)
            rc_flag = parse_rc(row)

            row_df, _ = process_row(
                amplicon=amplicon,
                input_seq=input_seq,
                window_tuple=window_tuple,
                rc_flag=rc_flag,
                base_dir=base_dir,
                variant_regex_chunk=variant_regex_chunk,
                print_variants=args.print_variants,
                print_orientation_equivalents=args.print_orientation_equivalents,
            )
            combined_rows.append(row_df)

        except Exception as e:
            print("\n[Row {}] Skipping due to error: {}".format(idx, e))

    combined_out = os.path.join(base_dir, "ABE_reads__combined.csv")
    if combined_rows:
        combined_df = pd.concat(combined_rows, axis=0, sort=True)
        combined_df = combined_df.reindex(sorted(combined_df.columns), axis=1)
        combined_df.to_csv(combined_out, index=True)
        print("\nWrote combined (transposed) output: {}".format(combined_out))
    else:
        print("\nNo successful rows to combine.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

