#!/usr/bin/env python3
"""Overlap clustering assignments with RNA-target classes.

Creates, per detector method:
  overlap_pam.xlsx          (one worksheet per k)
  overlap_hierarchical.xlsx (one worksheet per k)
  overlap_hdbscan.xlsx      (one worksheet: hdbscan)

Each worksheet contains four labelled tables: counts, cluster_composition,
class_distribution, and enrichment.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

MIN_N = 5


def discover_inputs(path: Path) -> list[Path]:
    if path.is_file():
        return [path]
    files = sorted(path.rglob("protein_clusters.csv"))
    if not files:
        raise FileNotFoundError(f"No protein_clusters.csv files under: {path}")
    return files


def safe_sheet_name(name: str) -> str:
    return re.sub(r"[\\/*?:\[\]]", "_", name)[:31]


def log2_haldane_or(a: int, b: int, c: int, d: int) -> float:
    return float(np.log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))))


def make_tables(df: pd.DataFrame, cluster_col: str, class_col: str) -> dict[str, pd.DataFrame]:
    data = df[[cluster_col, class_col]].dropna().copy()
    if data.empty:
        raise ValueError(f"No usable rows for {cluster_col}.")

    data[cluster_col] = data[cluster_col].astype(int)
    data[class_col] = data[class_col].astype(str)
    clusters = sorted(data[cluster_col].unique())
    classes = sorted(data[class_col].unique())

    counts = pd.crosstab(data[cluster_col], data[class_col]).reindex(
        index=clusters, columns=classes, fill_value=0
    )
    counts.index.name = "cluster"

    cluster_composition = counts.div(counts.sum(axis=1), axis=0)
    cluster_composition.index.name = "cluster"

    class_distribution = counts.div(counts.sum(axis=0), axis=1)
    class_distribution.index.name = "cluster"

    n_total = int(counts.to_numpy().sum())
    rows = []
    for cluster in clusters:
        cluster_size = int(counts.loc[cluster].sum())
        for rna_class in classes:
            a = int(counts.loc[cluster, rna_class]) # type: ignore
            b = cluster_size - a
            class_total = int(counts[rna_class].sum())
            c = class_total - a
            d = n_total - a - b - c
            expected = cluster_size * class_total / n_total
            try:
                p_value = float(fisher_exact([[a, b], [c, d]], alternative="two-sided")[1]) # type: ignore
            except ValueError:
                p_value = np.nan
            rows.append({
                "cluster": cluster,
                "rna_class": rna_class,
                "observed": a,
                "expected": expected,
                "cluster_size": cluster_size,
                "class_total": class_total,
                "cluster_composition": a / cluster_size if cluster_size else np.nan,
                "class_distribution": a / class_total if class_total else np.nan,
                "log2_obs_expected": np.log2(a / expected) if a > 0 and expected > 0 else np.nan,
                "log2_or": log2_haldane_or(a, b, c, d),
                "p_value": p_value,
                "descriptive_only": min(a, b, c, d) < MIN_N,
            })
    enrichment = pd.DataFrame(rows)
    enrichment["p_fdr"] = np.nan
    valid = enrichment["p_value"].notna()
    if valid.any():
        enrichment.loc[valid, "p_fdr"] = multipletests(
            enrichment.loc[valid, "p_value"], method="fdr_bh"
        )[1]
    enrichment = enrichment.sort_values(
        ["p_fdr", "cluster", "rna_class"], na_position="last"
    ).reset_index(drop=True)
    return {
        "counts": counts.reset_index(),
        "cluster_composition": cluster_composition.reset_index(),
        "class_distribution": class_distribution.reset_index(),
        "enrichment": enrichment,
    }


def write_sheet(writer: pd.ExcelWriter, sheet_name: str, tables: dict[str, pd.DataFrame]) -> None:
    row = 0
    for title in ["counts", "cluster_composition", "class_distribution", "enrichment"]:
        pd.DataFrame({title: []}).to_excel(writer, sheet_name=sheet_name, startrow=row, index=False)
        row += 1
        tables[title].to_excel(writer, sheet_name=sheet_name, startrow=row, index=False)
        row += len(tables[title]) + 3


def columns_for_algorithm(df: pd.DataFrame, algorithm: str, requested_k: set[int] | None) -> list[tuple[str, str]]:
    if algorithm == "hdbscan":
        return [("hdbscan", "cluster_hdbscan")] if "cluster_hdbscan" in df.columns else []
    prefix = f"cluster_{algorithm}_k"
    found = []
    for col in df.columns:
        match = re.fullmatch(re.escape(prefix) + r"(\d+)", col)
        if match:
            k = int(match.group(1))
            if requested_k is None or k in requested_k:
                found.append((f"k{k}", col))
    return sorted(found, key=lambda x: int(x[0][1:]))


def process_file(path: Path, output_root: Path, class_col: str, requested_k: set[int] | None) -> list[dict]:
    df = pd.read_csv(path)
    if class_col not in df.columns:
        raise KeyError(f"{path}: missing required RNA-class column '{class_col}'.")
    if "method" not in df.columns or df["method"].dropna().nunique() != 1:
        raise ValueError(f"{path}: expected exactly one non-null detector method in column 'method'.")

    method = str(df["method"].dropna().iloc[0])
    method_dir = output_root / method
    method_dir.mkdir(parents=True, exist_ok=True)
    manifest = []

    for algorithm in ["pam", "hierarchical", "hdbscan"]:
        runs = columns_for_algorithm(df, algorithm, requested_k)
        if not runs:
            print(f"{method}: no {algorithm} assignment columns; skipped.")
            continue

        workbook = method_dir / f"overlap_{algorithm}.xlsx"
        with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
            for sheet_name, cluster_col in runs:
                tables = make_tables(df, cluster_col, class_col)
                write_sheet(writer, safe_sheet_name(sheet_name), tables)
                manifest.append({
                    "method": method,
                    "algorithm": algorithm,
                    "sheet": sheet_name,
                    "cluster_column": cluster_col,
                    "n_proteins": int(df[[cluster_col, class_col]].dropna().shape[0]),
                    "n_clusters": int(df[cluster_col].nunique(dropna=True)),
                    "workbook": str(workbook),
                })
        print(f"Wrote {workbook}")
    return manifest


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default="lcr_analyses/clustering/carriers", help="A protein_clusters.csv file or a directory containing method subdirectories.")
    ap.add_argument("--outdir", default="lcr_analyses/clustering/carriers/rna_class_overlap")
    ap.add_argument("--class-col", default="rnaprimaryclass")
    ap.add_argument("--k-values", nargs="*", type=int, default=None, help="Optional PAM/hierarchical k values, e.g. --k-values 2 3 4 5 6")
    args = ap.parse_args()

    requested_k = set(args.k_values) if args.k_values else None
    output_root = Path(args.outdir)
    output_root.mkdir(parents=True, exist_ok=True)
    manifest = []
    for file in discover_inputs(Path(args.input)):
        manifest.extend(process_file(file, output_root, args.class_col, requested_k))
    pd.DataFrame(manifest).to_csv(output_root / "overlap_manifest.csv", index=False)
    print(f"Done. Outputs under: {output_root.resolve()}")


if __name__ == "__main__":
    main()
