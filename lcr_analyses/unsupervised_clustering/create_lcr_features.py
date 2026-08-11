#!/usr/bin/env python3
"""
Create LCR feature tables for method-specific downstream clustering.

Inputs:
  lcr_methods_combined.xlsx
  rbp_rna_classification-2.xlsx

Outputs:
  lcr_sequence_features.csv
  lcr_length_filter_summary.csv
"""

from collections import Counter
from pathlib import Path 
import math
import pandas as pd


LCR_FILE = r"rbp_lcrs\lcr_methods_combined.xlsx"
RNA_FILE = r"rbp_superclasses\rbp_rna_classification.xlsx"
OUTPUT_FILE = Path(r"lcr_analyses\sequence_level\lcr_sequence_features.xlsx")
MIN_LCR_LENGTH = 10

def get_uniprot_accession(protein_id: str) -> str:
    """Extract leading UniProt accession from a composite protein identifier."""
    return str(protein_id).split("|")[0]


def shannon_entropy(sequence: str) -> float:
    """Return Shannon entropy in bits."""
    counts = Counter(sequence)
    length = len(sequence)

    return -sum(
        (count / length) * math.log2(count / length)
        for count in counts.values()
    )


def dominant_residue_features(sequence: str) -> tuple[str, float]:
    """Return the dominant amino acid and its sequence fraction."""
    dominant_aa, dominant_count = Counter(sequence).most_common(1)[0]
    return dominant_aa, dominant_count / len(sequence)


def mmpr(sequence: str) -> int:
    """
    Minimum mutations needed to transform a sequence into a perfect repeat.

    Tests all repeat periods from 1 to floor(L / 2).
    """
    length = len(sequence)
    best_mutations = length - 1

    for period in range(1, length // 2 + 1):
        mutations = 0

        for phase in range(period):
            residues = sequence[phase::period]
            modal_count = Counter(residues).most_common(1)[0][1]
            mutations += len(residues) - modal_count

        best_mutations = min(best_mutations, mutations)

    return best_mutations


def dipeptide_dominance(sequence: str) -> float:
    """Return the fraction of the most common overlapping dipeptide."""
    if len(sequence) < 2:
        return 0.0

    dipeptides = [sequence[i:i + 2] for i in range(len(sequence) - 1)]
    modal_count = Counter(dipeptides).most_common(1)[0][1]

    return modal_count / len(dipeptides)


def calculate_features(sequence: str) -> pd.Series:
    """Calculate all selected LCR sequence-level features."""
    sequence = str(sequence).upper().strip()
    dominant_aa, dominant_fraction = dominant_residue_features(sequence)

    return pd.Series(
        {
            "dominant_aa": dominant_aa,
            "dominant_aa_fraction": dominant_fraction,
            "mmpr": mmpr(sequence),
            "shannon_entropy": shannon_entropy(sequence),
            "aa_richness": len(set(sequence)),
            "dipeptide_dominance": dipeptide_dominance(sequence),
        }
    )


def format_excel_sheet(writer, sheet_name: str, dataframe: pd.DataFrame) -> None:
    """Apply basic readable formatting to an Excel worksheet."""
    worksheet = writer.sheets[sheet_name]
    worksheet.freeze_panes(1, 0)
    worksheet.autofilter(0, 0, len(dataframe), len(dataframe.columns) - 1)

    for column_index, column_name in enumerate(dataframe.columns):
        max_length = max(
            len(str(column_name)),
            dataframe[column_name].astype(str).str.len().max(),
        )
        worksheet.set_column(column_index, column_index, min(max_length + 2, 45))


def main():
    lcrs = pd.read_excel(LCR_FILE, sheet_name="all_results")
    rna_classes = pd.read_excel(RNA_FILE, sheet_name="all_combined")

    rna_classes = rna_classes[
        [
            "uniprot_accession",
            "rna_primary_class",
            "rna_secondary_class",
        ]
    ].drop_duplicates(subset="uniprot_accession")

    lcrs["uniprot_accession"] = lcrs["protein_id"].map(
        get_uniprot_accession
    )

    summary = (
        lcrs.groupby("method", as_index=False)
        .size()
        .rename(columns={"size": "lcrs_before_length_filter"})
    )

    filtered_lcrs = lcrs[lcrs["length"] >= MIN_LCR_LENGTH].copy()

    after_filter = (
        filtered_lcrs.groupby("method", as_index=False)
        .size()
        .rename(columns={"size": "lcrs_after_length_filter"})
    )

    summary = summary.merge(after_filter, on="method", how="left")
    summary["lcrs_after_length_filter"] = (
        summary["lcrs_after_length_filter"]
        .fillna(0)
        .astype(int)
    )
    summary["retained_percent"] = (
        100
        * summary["lcrs_after_length_filter"]
        / summary["lcrs_before_length_filter"]
    ).round(2)

    filtered_lcrs = filtered_lcrs.merge(
        rna_classes,
        on="uniprot_accession",
        how="left",
    )

    features = filtered_lcrs["sequence"].apply(calculate_features)

    output = pd.concat(
        [
            filtered_lcrs[
                [
                    "method",
                    "protein_id",
                    "uniprot_accession",
                    "start",
                    "end",
                    "length",
                    "sequence",
                    "rna_primary_class",
                    "rna_secondary_class",
                ]
            ].reset_index(drop=True),
            features.reset_index(drop=True),
        ],
        axis=1,
    ).rename(
        columns={
            "start": "lcr_start",
            "end": "lcr_end",
            "length": "lcr_length",
        }
    )

    with pd.ExcelWriter(OUTPUT_FILE, engine="xlsxwriter") as writer:
        output.to_excel(writer, sheet_name="all_LCR_features", index=False)
        summary.to_excel(writer, sheet_name="filter_summary", index=False)

        format_excel_sheet(writer, "all_LCR_features", output)
        format_excel_sheet(writer, "filter_summary", summary)

        for method, method_table in output.groupby("method"):
            sheet_name = method[:31] # type: ignore
            method_table.to_excel(writer, sheet_name=sheet_name, index=False) # type: ignore
            format_excel_sheet(writer, sheet_name, method_table) # type: ignore

    print(f"Saved: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()