#!/usr/bin/env python3
"""
Summarise LCR sequence-composition features by calling method and RNA superclass.

Input:
    lcr_sequence_features.xlsx
    Expected sheet: all_LCR_features

Outputs:
    lcr_feature_summary.xlsx
    lcr_feature_summary_csv/<method>.csv

Each output row represents:
    one calling method × one RNA superclass × one feature
"""

from pathlib import Path
import re
import argparse
import pandas as pd


FEATURES = [
    "lcr_length",
    "dominant_aa_fraction",
    "mmpr",
    "shannon_entropy",
    "aa_richness",
    "dipeptide_dominance",
]

GROUP_COLUMNS = ["method", "rna_primary_class"]


def safe_filename(value: str) -> str:
    """Create a filesystem-safe filename from a calling-method name."""
    return re.sub(r"[^\w.-]+", "_", str(value)).strip("_")


def format_worksheet(writer, sheet_name: str, dataframe: pd.DataFrame) -> None:
    """Apply simple, readable Excel formatting."""
    worksheet = writer.sheets[sheet_name]
    worksheet.freeze_panes(1, 0)
    worksheet.autofilter(0, 0, len(dataframe), len(dataframe.columns) - 1)

    for i, column in enumerate(dataframe.columns):
        width = max(len(column), dataframe[column].astype(str).str.len().max()) + 2
        worksheet.set_column(i, i, min(width, 28))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarise LCR features by method and RNA superclass."
    )
    parser.add_argument(
        "input_file",
        help="Path to lcr_sequence_features.xlsx",
    )
    parser.add_argument(
        "--sheet",
        default="all_LCR_features",
        help="Input worksheet name (default: all_LCR_features)",
    )
    parser.add_argument(
        "--output-dir",
        default="lcr_composition_summary",
        help="Directory for XLSX and method-specific CSV files",
    )
    args = parser.parse_args()

    input_file = Path(args.input_file)
    output_dir = Path(args.output_dir)
    csv_dir = output_dir / "csv_by_method"
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_dir.mkdir(parents=True, exist_ok=True)

    data = pd.read_excel(input_file, sheet_name=args.sheet)

    required_columns = GROUP_COLUMNS + FEATURES
    missing_columns = sorted(set(required_columns) - set(data.columns))
    if missing_columns:
        raise ValueError(
            "Input file is missing required columns: "
            + ", ".join(missing_columns)
        )

    data = data.copy()
    data["rna_primary_class"] = data["rna_primary_class"].fillna("unknown")

    for feature in FEATURES:
        data[feature] = pd.to_numeric(data[feature], errors="coerce")

    long_data = data.melt(
        id_vars=GROUP_COLUMNS,
        value_vars=FEATURES,
        var_name="feature",
        value_name="value",
    ).dropna(subset=["value"])

    summary = (
        long_data.groupby(GROUP_COLUMNS + ["feature"], as_index=False)["value"]
        .agg(
            n="count",
            mean="mean",
            median="median",
            std="std",
            minimum="min",
            maximum="max",
        )
        .sort_values(GROUP_COLUMNS + ["feature"])
        .reset_index(drop=True)
    )

    numeric_columns = ["mean", "median", "std", "minimum", "maximum"]
    summary[numeric_columns] = summary[numeric_columns].round(4)

    xlsx_file = output_dir / "lcr_feature_summary.xlsx"
    with pd.ExcelWriter(xlsx_file, engine="xlsxwriter") as writer:
        summary.to_excel(writer, sheet_name="all_methods", index=False)
        format_worksheet(writer, "all_methods", summary)

        for method, method_summary in summary.groupby("method", sort=True):
            sheet_name = str(method)[:31]
            method_summary.to_excel(writer, sheet_name=sheet_name, index=False)
            format_worksheet(writer, sheet_name, method_summary)

    for method, method_summary in summary.groupby("method", sort=True):
        csv_file = csv_dir / f"{safe_filename(method)}_feature_summary.csv"
        method_summary.to_csv(csv_file, index=False)

    print(f"Saved Excel summary: {xlsx_file}")
    print(f"Saved method-specific CSV files: {csv_dir}")


if __name__ == "__main__":
    main()