#!/usr/bin/env python3
"""Join a FASTA file to RNA-class annotations and plot sequence composition.

Usage:
    python plot_rna_binder_sequence_composition.py proteins.fasta annotations.xlsx

Outputs are written to ./rna_binder_plots/
"""

from argparse import ArgumentParser
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio import SeqIO

EXCLUDED_CLASSES = {"unknown"}


def read_fasta_metrics(fasta_file):
    """Return one row per UniProt accession in FASTA headers."""
    rows = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        accession = record.description.split("|")[0].lstrip(">")
        sequence = str(record.seq).upper().replace("*", "")
        length = len(sequence)
        if length == 0:
            continue

        counts = Counter(sequence)
        dominant_aa, dominant_count = counts.most_common(1)[0]
        rows.append(
            {
                "uniprot_accession": accession,
                "sequence_length": length,
                "dominant_amino_acid": dominant_aa,
                "dominant_amino_acid_fraction": dominant_count / length,
            }
        )
    return pd.DataFrame(rows).drop_duplicates("uniprot_accession")


def scatter_by_class(data, class_column, title, output_file):
    plot_data = data.dropna(subset=[class_column]).copy()
    plot_data[class_column] = plot_data[class_column].astype(str).str.strip()
    plot_data = plot_data[plot_data[class_column] != ""]
    plot_data = plot_data[~plot_data[class_column].isin(EXCLUDED_CLASSES)]

    classes = sorted(plot_data[class_column].unique())
    colors = plt.cm.tab20.colors # type: ignore

    fig, ax = plt.subplots(figsize=(10, 7))
    for i, class_name in enumerate(classes):
        subset = plot_data[plot_data[class_column] == class_name]
        ax.scatter(
            np.log10(subset["sequence_length"]),
            subset["dominant_amino_acid_fraction"],
            s=28,
            alpha=0.7,
            color=colors[i % len(colors)],
            label=class_name,
            edgecolors="none",
        )

    ax.set_title(title)
    ax.set_xlabel("log10(Protein length in amino acids)")
    ax.set_ylabel("Most dominant amino-acid fraction")
    ax.grid(alpha=0.25)
    ax.legend(title=class_column, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = ArgumentParser()
    parser.add_argument("fasta", help="Input FASTA with UniProt accession first in each header")
    parser.add_argument("xlsx", help="Excel file containing the all_combined sheet")
    parser.add_argument("--outdir", default="rna_binder_plots", help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sequence_metrics = read_fasta_metrics(args.fasta)
    annotations = pd.read_excel(args.xlsx, sheet_name="all_combined")

    required_columns = {
        "uniprot_accession",
        "rna_primary_class",
        "rna_secondary_class",
        "confidence",
        "evidence_source",
        "evidence_matched",
    }
    missing = required_columns - set(annotations.columns)
    if missing:
        raise ValueError(f"Missing columns in all_combined: {', '.join(sorted(missing))}")

    annotations["uniprot_accession"] = annotations["uniprot_accession"].astype(str).str.strip()
    merged = annotations.merge(sequence_metrics, on="uniprot_accession", how="inner")
    merged.to_csv(outdir / "rna_binder_sequence_metrics.csv", index=False)

    scatter_by_class(
        merged,
        "rna_primary_class",
        "Protein length vs dominant amino-acid fraction by primary RNA class",
        outdir / "scatter_primary_class.png",
    )
    scatter_by_class(
        merged,
        "rna_secondary_class",
        "Protein length vs dominant amino-acid fraction by secondary RNA class",
        outdir / "scatter_secondary_class.png",
    )

    merged["primary_secondary_class"] = (
        merged["rna_primary_class"].astype(str).str.strip()
        + " | "
        + merged["rna_secondary_class"].astype(str).str.strip()
    )
    scatter_by_class(
        merged,
        "primary_secondary_class",
        "Protein length vs dominant amino-acid fraction by primary × secondary class",
        outdir / "scatter_primary_secondary_intersection.png",
    )

    print(f"FASTA entries with valid sequences: {len(sequence_metrics):,}")
    print(f"Entries matched to Excel annotations: {len(merged):,}")
    print(f"Saved CSV and three plots in: {outdir.resolve()}")


if __name__ == "__main__":
    main()
