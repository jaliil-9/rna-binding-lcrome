#!/usr/bin/env python3
"""
Create expanded LCR feature tables (v2) for downstream clustering.

Feature groups (Phase 2.1.2):
1. Amino acid composition      -- frac_A .. frac_Y
2. Entropy                     -- shannon_entropy, norm_entropy (length-corrected)
3. Periodicity                 -- mmpr_norm, best_period, dipeptide_dominance
4. Physicochemical composition -- property-class fractions (same AA sets as
                                  pc_analysis_annotation_v2_2.py)
5. Flank/end position          -- relative to protein length (% rules)

Protein length: taken from the RBP annotation workbook (--lengths-file,
sheet "Combined", columns uniprot_accession + uniprot_length). Proteins
missing from the mapping fall back to max(lcr_end) across all methods;
the source used is recorded in the protein_length_source column.

Inputs:
    lcr_methods_combined.xlsx        (sheet: all_results)
    rbp_rna_classification.xlsx      (sheet: all_combined)

Outputs:
    lcr_sequence_features_v2.xlsx
"""

import argparse
import math
from collections import Counter
from pathlib import Path

import pandas as pd

LCR_FILE = r"rbp_lcrs\lcr_methods_combined.xlsx"
RNA_FILE = r"rbp_superclasses\rbp_rna_classification.xlsx"
LENGTHS_FILE = r"datasets\combined_rbp_pfam38_rbpdb_modomics_uniprot.xlsx"
LENGTHS_SHEET = "Combined"
OUTPUT_FILE = Path(r"lcr_analyses\lcr_sequence_features_v2.xlsx")
MIN_LCR_LENGTH = 10

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# Property-class residue sets (aligned with aa-physicochemical-properties-v2.csv
# conventions used in pc_analysis_annotation_v2_2.py)
POSITIVE = set("KR")          # His excluded (largely uncharged at pH 7.4)
NEGATIVE = set("DE")
POLAR = set("STNQ")           # uncharged polar
STRONG_HYDRO = set("VILMFWY") # v2.2 aggregation rule set
AROMATIC = set("FWY")
DISORDER = set("PQESRKGDNA")  # disorder-promoting residues


def get_uniprot_accession(protein_id: str) -> str:
    """Extract leading UniProt accession from a composite protein identifier."""
    return str(protein_id).split("|")[0]


# ---------------------------------------------------------------------------
# Group 2: entropy
# ---------------------------------------------------------------------------

def shannon_entropy(sequence: str) -> float:
    """Return Shannon entropy in bits."""
    counts = Counter(sequence)
    length = len(sequence)
    return -sum(
        (count / length) * math.log2(count / length)
        for count in counts.values()
    )


def norm_entropy(sequence: str) -> float:
    """Entropy normalized by the maximum achievable at this length."""
    length = len(sequence)
    if length < 2:
        return 0.0
    max_entropy = math.log2(min(length, 20))
    return shannon_entropy(sequence) / max_entropy


# ---------------------------------------------------------------------------
# Group 1: amino acid composition
# ---------------------------------------------------------------------------

def dominant_residue_features(sequence: str) -> tuple[str, float]:
    """Return the dominant amino acid and its sequence fraction."""
    dominant_aa, dominant_count = Counter(sequence).most_common(1)[0]
    return dominant_aa, dominant_count / len(sequence)


def aa_fractions(sequence: str) -> dict:
    """Return the 20 amino acid fractions."""
    counts = Counter(sequence)
    length = len(sequence)
    return {f"frac_{aa}": counts.get(aa, 0) / length for aa in AA_ALPHABET}


# ---------------------------------------------------------------------------
# Group 3: periodicity
# ---------------------------------------------------------------------------

def mmpr_with_period(sequence: str) -> tuple[int, int]:
    """
    Minimum mutations to a perfect repeat, and the period achieving it.

    Tests all repeat periods from 1 to floor(L / 2).
    """
    length = len(sequence)
    best_mutations = length - 1
    best_period = 1

    for period in range(1, length // 2 + 1):
        mutations = 0
        for phase in range(period):
            residues = sequence[phase::period]
            modal_count = Counter(residues).most_common(1)[0][1]
            mutations += len(residues) - modal_count
        if mutations < best_mutations:
            best_mutations = mutations
            best_period = period

    return best_mutations, best_period


def dipeptide_dominance(sequence: str) -> float:
    """Return the fraction of the most common overlapping dipeptide."""
    if len(sequence) < 2:
        return 0.0
    dipeptides = [sequence[i:i + 2] for i in range(len(sequence) - 1)]
    modal_count = Counter(dipeptides).most_common(1)[0][1]
    return modal_count / len(dipeptides)


# ---------------------------------------------------------------------------
# Group 4: physicochemical composition (fractions only, no threshold rules)
# ---------------------------------------------------------------------------

def class_fraction(sequence: str, residue_set: set) -> float:
    """Fraction of residues belonging to a property class."""
    return sum(1 for res in sequence if res in residue_set) / len(sequence)


def physicochemical_fractions(sequence: str) -> dict:
    """Property-class fractions and Das & Pappu charge summaries."""
    f_pos = class_fraction(sequence, POSITIVE)
    f_neg = class_fraction(sequence, NEGATIVE)
    return {
        "frac_positive": f_pos,
        "frac_negative": f_neg,
        "fcr": f_pos + f_neg,
        "ncpr": f_pos - f_neg,
        "frac_polar": class_fraction(sequence, POLAR),
        "frac_strong_hydro": class_fraction(sequence, STRONG_HYDRO),
        "frac_aromatic": class_fraction(sequence, AROMATIC),
        "frac_disorder": class_fraction(sequence, DISORDER),
    }


# ---------------------------------------------------------------------------
# Sequence-level driver
# ---------------------------------------------------------------------------

def calculate_features(sequence: str) -> pd.Series:
    """Calculate all sequence-level LCR features (groups 1-4)."""
    sequence = str(sequence).upper().strip()
    dominant_aa, dominant_fraction = dominant_residue_features(sequence)
    mmpr, best_period = mmpr_with_period(sequence)

    features = {
        "dominant_aa": dominant_aa,
        "dominant_aa_fraction": dominant_fraction,
        "shannon_entropy": shannon_entropy(sequence),
        "norm_entropy": norm_entropy(sequence),
        "aa_richness": len(set(sequence)),
        "dipeptide_dominance": dipeptide_dominance(sequence),
        "mmpr": mmpr,
        "mmpr_norm": mmpr / len(sequence),
        "best_period": best_period,
    }
    features.update(aa_fractions(sequence))
    features.update(physicochemical_fractions(sequence))
    return pd.Series(features)


# ---------------------------------------------------------------------------
# Group 5: flank/end position (table level, needs protein length)
# ---------------------------------------------------------------------------

def add_position_features(table: pd.DataFrame,
                          lengths_file: Path,
                          lengths_sheet: str) -> pd.DataFrame:
    """
    Add protein-length-relative position features.

    protein_length    -- uniprot_length from the annotation workbook where
                         available, else max(lcr_end) per protein (proxy)
    rel_start/rel_end -- LCR boundaries as fraction of protein length
    lcr_coverage      -- LCR length as fraction of protein length
    terminus_dist_pct -- distance to the nearest terminus, per length:
                         min(lcr_start - 1, protein_length - lcr_end) / pl
    """
    lengths = pd.read_excel(lengths_file, sheet_name=lengths_sheet)
    length_map = (
        lengths[["uniprot_accession", "uniprot_length"]]
        .dropna()
        .drop_duplicates(subset="uniprot_accession")
        .set_index("uniprot_accession")["uniprot_length"]
    )

    table["protein_length"] = table["uniprot_accession"].map(length_map)
    proxy = table.groupby("uniprot_accession")["lcr_end"].transform("max")

    table["protein_length_source"] = "uniprot"
    missing = table["protein_length"].isna()
    table.loc[missing, "protein_length"] = proxy[missing]
    table.loc[missing, "protein_length_source"] = "proxy_max_lcr_end"

    pl = table["protein_length"]
    table["rel_start"] = table["lcr_start"] / pl
    table["rel_end"] = table["lcr_end"] / pl
    table["lcr_coverage"] = table["lcr_length"] / pl
    table["terminus_dist_pct"] = (
        pd.concat([table["lcr_start"] - 1, pl - table["lcr_end"]], axis=1)
        .min(axis=1) / pl
    )

    n_proxy = int(missing.sum())
    if n_proxy:
        print(f"WARNING: {n_proxy} LCRs use the max-lcr-end length proxy "
              f"(accession missing from {lengths_file.name}).")
    return table


# ---------------------------------------------------------------------------
# Output formatting (unchanged from v1)
# ---------------------------------------------------------------------------

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
    parser = argparse.ArgumentParser()
    parser.add_argument("--lcr-file", type=Path, default=Path(LCR_FILE))
    parser.add_argument("--rna-file", type=Path, default=Path(RNA_FILE))
    parser.add_argument(
        "--lengths-file",
        type=Path,
        default=Path(LENGTHS_FILE),
        help="RBP annotation workbook with a sheet containing "
             "uniprot_accession and uniprot_length columns.",
    )
    parser.add_argument("--lengths-sheet", default=LENGTHS_SHEET)
    parser.add_argument("--output-file", type=Path, default=OUTPUT_FILE)
    args = parser.parse_args()

    lcrs = pd.read_excel(args.lcr_file, sheet_name="all_results")
    rna_classes = pd.read_excel(args.rna_file, sheet_name="all_combined")

    rna_classes = rna_classes[
        ["uniprot_accession", "rna_primary_class", "rna_secondary_class"]
    ].drop_duplicates(subset="uniprot_accession")

    lcrs["uniprot_accession"] = lcrs["protein_id"].map(get_uniprot_accession)

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
        summary["lcrs_after_length_filter"].fillna(0).astype(int)
    )
    summary["retained_percent"] = (
        100 * summary["lcrs_after_length_filter"]
        / summary["lcrs_before_length_filter"]
    ).round(2)

    filtered_lcrs = filtered_lcrs.merge(
        rna_classes, on="uniprot_accession", how="left"
    )

    features = filtered_lcrs["sequence"].apply(calculate_features)

    output = pd.concat(
        [
            filtered_lcrs[
                [
                    "method", "protein_id", "uniprot_accession",
                    "start", "end", "length", "sequence",
                    "rna_primary_class", "rna_secondary_class",
                ]
            ].reset_index(drop=True),
            features.reset_index(drop=True),
        ],
        axis=1,
    ).rename(
        columns={"start": "lcr_start", "end": "lcr_end", "length": "lcr_length"}
    )

    output = add_position_features(output, args.lengths_file,
                                   args.lengths_sheet)

    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(args.output_file, engine="xlsxwriter") as writer:
        output.to_excel(writer, sheet_name="all_LCR_features", index=False)
        summary.to_excel(writer, sheet_name="filter_summary", index=False)
        format_excel_sheet(writer, "all_LCR_features", output)
        format_excel_sheet(writer, "filter_summary", summary)

        for method, method_table in output.groupby("method"):
            sheet_name = method[:31]  # type: ignore
            method_table.to_excel(writer, sheet_name=sheet_name, index=False)  # type: ignore
            format_excel_sheet(writer, sheet_name, method_table)  # type: ignore

    print(f"Saved: {args.output_file}")


if __name__ == "__main__":
    main()
