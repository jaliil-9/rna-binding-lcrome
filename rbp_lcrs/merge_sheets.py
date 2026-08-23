#!/usr/bin/env python3
"""Combine CSV files into a single XLSX workbook.

Sheet 1 ("Combined"): all rows from all input files appended together.
Remaining sheets: one per input CSV, named after the file (without extension).

Usage:
    python csv_to_xlsx.py seg.csv cast.csv pfam_lcr.csv fLPS.csv PlaToLoCo.csv dot2d.csv
    python csv_to_xlsx.py *.csv -o results.xlsx
"""

import argparse
import re
from pathlib import Path

import pandas as pd

INVALID_SHEET_CHARS = re.compile(r"[\[\]:*?/\\]")
MAX_SHEET_NAME = 31  # Excel limit


def sheet_name_from(path: Path, used: set) -> str:
    """Derive a valid, unique Excel sheet name from the CSV filename."""
    name = INVALID_SHEET_CHARS.sub("_", path.stem)[:MAX_SHEET_NAME] or "Sheet"
    base, i = name, 1
    while name in used or name.lower() == "combined":
        suffix = f"_{i}"
        name = base[: MAX_SHEET_NAME - len(suffix)] + suffix
        i += 1
    used.add(name)
    return name


def main():
    parser = argparse.ArgumentParser(
        description="Combine CSVs into one XLSX with a 'Combined' sheet plus one sheet per file."
    )
    parser.add_argument("csvs", nargs="+", type=Path, help="Input CSV files (e.g. your 6 method files)")
    parser.add_argument("-o", "--output", type=Path, default=Path("combined.xlsx"),
                        help="Output workbook path (default: combined.xlsx)")
    args = parser.parse_args()

    frames = {}
    used_names = set()
    for csv_path in args.csvs:
        df = pd.read_csv(csv_path)
        frames[sheet_name_from(csv_path, used_names)] = df

    combined = pd.concat(frames.values(), ignore_index=True)

    with pd.ExcelWriter(args.output, engine="openpyxl") as writer:
        combined.to_excel(writer, sheet_name="Combined", index=False)
        for name, df in frames.items():
            df.to_excel(writer, sheet_name=name, index=False)

    print(f"Wrote '{args.output}' with {len(frames) + 1} sheets "
          f"({len(combined)} total rows in 'Combined').")


if __name__ == "__main__":
    main()
