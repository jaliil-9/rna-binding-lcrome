#!/usr/bin/env python3
"""Summarise LCR sequence features by method and position class.

Input: lcr_position_classes.csv from classify_lcr_positions.py.
Optional metadata: CSV/XLSX containing protein_id (or protein_accession) and
rna_primary_class. This is required only for RNA-superclass distributions and
enrichment if that column is not already in the input.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

AA = list("ACDEFGHIKLMNPQRSTVWY")
MOTIFS = {
    "RGG_RGGbox": r"RGG",
    "RG": r"RG",
    "SR_RS": r"SR|RS",
    "YR": r"YR",
    "FG": r"FG",
    "PXXP": r"P..P",
}


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_excel(path) if path.suffix.lower() in {".xlsx", ".xls"} else pd.read_csv(path)


def protein_accession(value: str) -> str:
    return str(value).split("|")[0]


def shannon_entropy(sequence: str) -> float:
    counts = pd.Series(list(sequence)).value_counts()
    p = counts / counts.sum()
    return float(-(p * np.log2(p)).sum())


def add_sequence_features(df: pd.DataFrame) -> pd.DataFrame:
    seq = df["sequence"].fillna("").str.upper()
    df = df.copy()
    df["shannon_entropy"] = seq.map(shannon_entropy)
    for aa in AA:
        df[f"frac_{aa}"] = seq.str.count(aa) / df["lcr_length"]
    for name, pattern in MOTIFS.items():
        df[f"motif_{name}"] = seq.str.contains(pattern, regex=True, na=False)
    return df


def merge_superclass(df: pd.DataFrame, metadata_file: Path | None) -> pd.DataFrame:
    if "rna_primary_class" in df.columns:
        return df
    if metadata_file is None:
        return df

    meta = pd.read_excel(metadata_file, sheet_name="all_combined").copy()
    required = {"uniprot_accession", "rna_primary_class"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"Metadata sheet all_combined is missing: {sorted(missing)}")

    meta["protein_accession"] = meta["uniprot_accession"].astype(str).str.split("|").str[0]
    meta = meta[["protein_accession", "rna_primary_class"]].drop_duplicates("protein_accession")
    return df.merge(meta, on="protein_accession", how="left")


def main(input_file: Path, output_dir: Path, metadata_file: Path | None) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(input_file)
    required = {"method", "primary_class", "sequence", "lcr_length", "protein_accession"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    df = merge_superclass(df, metadata_file)
    df = add_sequence_features(df)
    group = ["method", "primary_class"]

    length_entropy = (df.groupby(group).agg(
        n_lcrs=("sequence", "size"),
        median_length=("lcr_length", "median"), mean_length=("lcr_length", "mean"),
        median_entropy=("shannon_entropy", "median"), mean_entropy=("shannon_entropy", "mean")
    ).reset_index())
    aa_columns = [f"frac_{aa}" for aa in AA]
    composition = (df.groupby(group)[aa_columns].mean().reset_index()
                   .melt(id_vars=group, var_name="amino_acid", value_name="mean_fraction"))
    composition["amino_acid"] = composition["amino_acid"].str.replace("frac_", "", regex=False)
    motif_columns = [f"motif_{name}" for name in MOTIFS]
    motifs = (df.groupby(group)[motif_columns].mean().reset_index()
              .melt(id_vars=group, var_name="motif", value_name="fraction_lcrs_with_motif"))
    motifs["motif"] = motifs["motif"].str.replace("motif_", "", regex=False)

        # Write every result table to one Excel workbook.
    workbook = output_dir / "lcr_position_class_analysis.xlsx"

    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        length_entropy.to_excel(writer, sheet_name="length_entropy", index=False)
        composition.to_excel(writer, sheet_name="aa_composition", index=False)
        motifs.to_excel(writer, sheet_name="motif_prevalence", index=False)

        if "rna_primary_class" in df.columns:
            valid = df.dropna(subset=["rna_primary_class"]).copy()

            observed = (
                valid.groupby(["method", "primary_class", "rna_primary_class"])
                .size()
                .reset_index(name="observed_n")
            )

            class_total = (
                valid.groupby(group)
                .size()
                .rename("class_n")
                .reset_index()
            )

            background = (
                valid.groupby(["method", "rna_primary_class"])
                .size()
                .rename("method_superclass_n")
                .reset_index()
            )

            method_total = (
                valid.groupby("method")
                .size()
                .rename("method_n")
                .reset_index()
            )

            enrichment = (
                observed
                .merge(class_total, on=group)
                .merge(background, on=["method", "rna_primary_class"])
                .merge(method_total, on="method")
            )

            enrichment["observed_fraction"] = (
                enrichment["observed_n"] / enrichment["class_n"]
            )

            enrichment["background_fraction"] = (
                enrichment["method_superclass_n"] / enrichment["method_n"]
            )

            enrichment["fold_enrichment"] = (
                enrichment["observed_fraction"]
                / enrichment["background_fraction"]
            )

            enrichment.to_excel(
                writer,
                sheet_name="superclass_enrichment",
                index=False
            )

            print("RNA-superclass analysis written.")
        else:
            print("No RNA superclass column: skipped superclass enrichment.")

    print(f"Wrote {workbook}")
    print(f"Analysed {len(df):,} method-specific LCRs.")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, default=Path(r"lcr_analyses\domain_function\lcr_position_classes.csv"))
    p.add_argument("--output-dir", type=Path, default=Path(r"lcr_analyses\domain_function\class_feature_analysis"))
    p.add_argument("--metadata", type=Path, default=Path(r"rbp_superclasses\rbp_rna_classification.xlsx"))
    args = p.parse_args()
    main(args.input, args.output_dir, args.metadata)
