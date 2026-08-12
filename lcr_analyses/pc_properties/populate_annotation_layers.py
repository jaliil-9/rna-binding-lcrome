#!/usr/bin/env python3
"""
Populate empty rna_target_superclass and domain_position_class columns in
lcr_annotations.xlsx from a source xlsx ("Sheet1").

Mapping:
    rna_target_superclass  <- rna_primary_class
    domain_position_class  <- primary_class

Join key: protein_id + method + start + end (LCR-level unique).

Usage:
    python populate_annotation_layers.py lcr_annotations.xlsx source.xlsx
    python populate_annotation_layers.py lcr_annotations.xlsx source.xlsx -o out.xlsx
"""

import sys
import argparse
import pandas as pd

KEY = ["protein_id", "method", "start", "end"]
MAPPING = {
    "rna_target_superclass": "rna_primary_class",
    "domain_position_class": "primary_class",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("annotations", help="lcr_annotations.xlsx to fill")
    ap.add_argument("source", help="source xlsx containing Sheet1")
    ap.add_argument("-o", "--output", default=None,
                    help="output file (default: <annotations>_filled.xlsx)")
    ap.add_argument("--sheet", default="Sheet1", help="sheet in source file")
    args = ap.parse_args()

    out = args.output or args.annotations.replace(".xlsx", "_filled.xlsx")

    ann = pd.read_excel(args.annotations)
    src = pd.read_excel(args.source, sheet_name=args.sheet)

    # Normalize join-key dtypes on both sides
    for df in (ann, src):
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
        df["protein_id"] = df["protein_id"].astype(str)
        df["method"] = df["method"].astype(str)

    keep = KEY + list(MAPPING.values())
    missing = [c for c in keep if c not in src.columns]
    if missing:
        sys.exit(f"Source file is missing columns: {missing}")

    src_small = src[keep].drop_duplicates(subset=KEY)
    n_dup = len(src[keep]) - len(src_small)
    if n_dup:
        print(f"Warning: dropped {n_dup} duplicate key rows in source (kept first)")

    merged = ann.merge(src_small, on=KEY, how="left", suffixes=("", "_src"))

    filled = {t: 0 for t in MAPPING}
    for target, source_col in MAPPING.items():
        # object dtype required: pandas >=3 refuses str into all-NaN float col
        if target not in merged.columns:
            merged[target] = pd.Series([None] * len(merged), dtype="object")
        else:
            merged[target] = merged[target].astype("object")

        empty = merged[target].isna() | (
            merged[target].astype(str).str.strip().isin(["", "nan", "None"])
        )
        has_value = merged[source_col].notna()
        mask = empty & has_value
        merged.loc[mask, target] = merged.loc[mask, source_col].astype("object")
        filled[target] = int(mask.sum())
        merged = merged.drop(columns=[source_col])

    merged.to_excel(out, index=False)

    print(f"Rows: {len(merged)}")
    for t, n in filled.items():
        still_empty = int(
            (merged[t].isna()
             | merged[t].astype(str).str.strip().isin(["", "nan", "None"])).sum()
        )
        print(f"  {t}: filled {n} rows, still empty {still_empty}")
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
