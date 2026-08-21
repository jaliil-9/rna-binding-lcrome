#!/usr/bin/env python3
"""Convert a UniProt-style FASTA file into a reference XLSX.

Output columns: uniprot_accession, uniprot_sequence, uniprot_length
Output sheet:   "Combined"  (what merge_redundant_lcr.py expects as REF_SHEET)

Accession extraction: for headers like >sp|A0A087X1C5|CP2D7_HUMAN ... the
second pipe-delimited field is used. For non-piped headers, the first
whitespace-delimited token is used instead.

Usage:
    python fasta_to_reference_xlsx.py proteins.fasta
    python fasta_to_reference_xlsx.py proteins.fasta -o reference.xlsx
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_fasta(path):
    """Yield (header, sequence) tuples; sequence lines joined, whitespace stripped."""
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def accession_from(header):
    """UniProt header >sp|ACC|NAME ... -> ACC; otherwise first token."""
    parts = header.split("|")
    if len(parts) >= 2 and parts[0] in {"sp", "tr"}:
        return parts[1]
    return header.split()[0]


def main():
    parser = argparse.ArgumentParser(
        description="Convert FASTA to a reference XLSX (sheet 'Combined') "
                    "compatible with merge_redundant_lcr.py."
    )
    parser.add_argument("fasta", type=Path, help="Input FASTA file")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="Output XLSX path (default: <fasta_stem>_reference.xlsx)")
    args = parser.parse_args()

    out_path = args.output or args.fasta.with_name(args.fasta.stem + "_reference.xlsx")

    rows = []
    for header, seq in parse_fasta(args.fasta):
        seq = seq.upper()
        rows.append({
            "uniprot_accession": accession_from(header),
            "uniprot_sequence": seq,
            "uniprot_length": len(seq),
        })

    df = pd.DataFrame(rows, columns=["uniprot_accession", "uniprot_sequence", "uniprot_length"])

    n_dupes = df["uniprot_accession"].duplicated().sum()
    if n_dupes:
        # merge_redundant_lcr.py drops duplicates keeping the first anyway;
        # doing it here makes the behavior explicit and warns the user.
        dupes = df.loc[df["uniprot_accession"].duplicated(), "uniprot_accession"].tolist()
        print(f"WARNING: {n_dupes} duplicate accession(s), keeping first occurrence: "
              f"{dupes[:10]}{'...' if n_dupes > 10 else ''}")
        df = df.drop_duplicates(subset="uniprot_accession", keep="first")

    df.to_excel(out_path, sheet_name="Combined", index=False)
    print(f"Wrote {out_path} ({len(df)} entries, sheet 'Combined').")


if __name__ == "__main__":
    main()
