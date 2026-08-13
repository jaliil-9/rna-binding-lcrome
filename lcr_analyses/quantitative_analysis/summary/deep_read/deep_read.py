#!/usr/bin/env python3
"""Phase 2.4.1 deep-read extraction layer.

Reads the existing cross-method master tables produced by
phase_2_4_1_consolidate.py and emits three focused extraction tables plus two
charge-landscape figures. No recomputation of the underlying analyses.

Inputs  (under <outdir>/summary/tables/):
  master_lcr_continuous.csv, master_protein_continuous.csv,
  master_lcr_categorical.csv, master_protein_categorical.csv

Outputs (under <outdir>/summary/deep_read/):
  tables/
    charge_landscape_summary.csv        - per class x method: median frac_positive,
                                          frac_negative, fcr, ncpr + KW p/FDR (LCR level)
    charge_landscape_median_<metric>.csv- class x method median pivots (4 files)
    significant_coverage_cells.csv      - coverage + property-fraction KW cells with
                                          FDR < 0.05 (protein level, plus LCR-level
                                          coverage_per_lcr for reference)
    concordance_check.csv               - per method x feature x class: LCR-level log2
                                          (sequence-count and residue-count bases) next
                                          to protein-level log2, with concordance flags
  figures/
    ncpr_median_heatmap.png
    frac_positive_median_heatmap.png
  README_deep_read.md

Run:
    python phase_2_4_1_deep_read.py --outdir lcr_analyses/results
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_theme(style="whitegrid", context="notebook")

CHARGE_METRICS = ["frac_positive", "frac_negative", "fcr", "ncpr"]
FDR_ALPHA = 0.05
STRONG_LOG2 = 1.0  # |log2(obs/exp)| >= 1 means at least two-fold enrichment/depletion


def load_tables(root: Path) -> dict[str, pd.DataFrame]:
    tables_dir = root / "summary" / "tables"
    out = {}
    for key, fname in {
        "lcr_cont": "master_lcr_continuous.csv",
        "protein_cont": "master_protein_continuous.csv",
        "lcr_cat": "master_lcr_categorical.csv",
        "protein_cat": "master_protein_categorical.csv",
    }.items():
        path = tables_dir / fname
        if path.exists():
            out[key] = pd.read_csv(path)
        else:
            print(f"WARNING: missing {path} — related outputs will be skipped.")
    return out


def feature_value_column(df: pd.DataFrame) -> pd.Series:
    """Coalesce the analysis-dependent feature column into one series."""
    value = pd.Series(np.nan, index=df.index, dtype=object)
    for col in ["signature", "position", "has_pfam", "n_lcr_category"]:
        if col in df.columns:
            value = value.fillna(df[col])
    return value


def charge_landscape(lcr_cont: pd.DataFrame, tables: Path, figures: Path) -> None:
    sub = lcr_cont.loc[lcr_cont["metric"].isin(CHARGE_METRICS)].copy()
    if sub.empty:
        print("No charge metrics found in master_lcr_continuous.csv — skipping charge landscape.")
        return
    cols = ["method", "rna_class", "metric", "n", "median", "q1", "q3", "iqr", "mean", "sd", "p_value", "p_fdr"]
    sub[cols].sort_values(["metric", "method", "rna_class"]).to_csv(tables / "charge_landscape_summary.csv", index=False)
    for metric in CHARGE_METRICS:
        mat = sub.loc[sub["metric"] == metric].pivot_table(index="rna_class", columns="method", values="median", aggfunc="first")
        mat.to_csv(tables / f"charge_landscape_median_{metric}.csv")
        if metric in {"ncpr", "frac_positive"} and not mat.empty:
            fig, ax = plt.subplots(figsize=(max(6, 1.2 * mat.shape[1] + 2), max(4, .5 * mat.shape[0] + 1.5)))
            sns.heatmap(mat, cmap="vlag" if metric == "ncpr" else "viridis", center=0 if metric == "ncpr" else None,
                        annot=True, fmt=".2f", linewidths=.4, linecolor="white", ax=ax, cbar_kws={"label": f"median {metric}"})
            ax.set(title=f"Median {metric} of LCRs: RNA class x method", xlabel="Method", ylabel="RNA-target superclass")
            fig.tight_layout(); fig.savefig(figures / f"{metric}_median_heatmap.png", dpi=220); plt.close(fig)


def significant_coverage(lcr_cont: pd.DataFrame, protein_cont: pd.DataFrame, tables: Path) -> None:
    frames = []
    if protein_cont is not None:
        p = protein_cont.loc[protein_cont["metric"].eq("coverage") | protein_cont["metric"].str.startswith("property_fraction__", na=False)].copy()
        p["level"] = "protein"
        frames.append(p)
    if lcr_cont is not None:
        l = lcr_cont.loc[lcr_cont["metric"] == "coverage_per_lcr"].copy()
        l["level"] = "lcr"
        frames.append(l)
    if not frames:
        print("No continuous masters found — skipping significant coverage cells.")
        return
    cols = ["level", "method", "metric", "rna_class", "n", "median", "q1", "q3", "iqr", "mean", "sd", "p_value", "p_fdr"]
    out = pd.concat(frames, ignore_index=True)
    out = out.loc[out["p_fdr"] < FDR_ALPHA, cols].sort_values(["metric", "method", "rna_class"])
    out.to_csv(tables / "significant_coverage_cells.csv", index=False)
    print(f"Significant coverage/property-fraction cells (FDR < {FDR_ALPHA}): {len(out)}")


def concordance_for(feature_type: str, lcr_cat: pd.DataFrame, protein_cat: pd.DataFrame) -> pd.DataFrame:
    """Side-by-side LCR-level (two counting bases) and protein-level log2 for one feature type."""
    lcr_analyses = {"signature": "signature_sequence_count_lcr", "position": "position_lcr"}[feature_type]
    lcr_residue = {"signature": "signature_residue_count_lcr", "position": None}[feature_type]
    protein_analysis = {"signature": "signature_presence", "position": "position_presence"}[feature_type]

    lcr_seq = lcr_cat.loc[lcr_cat["analysis"] == lcr_analyses].copy()
    lcr_seq["feature_value"] = feature_value_column(lcr_seq)
    lcr_seq = lcr_seq.rename(columns={"log2_obs_exp": "lcr_log2_sequence_count", "observed": "lcr_observed"})
    lcr_seq = lcr_seq[["method", "feature_value", "rna_class", "lcr_log2_sequence_count", "lcr_observed"]]

    if lcr_residue:
        lcr_res = lcr_cat.loc[lcr_cat["analysis"] == lcr_residue].copy()
        lcr_res["feature_value"] = feature_value_column(lcr_res)
        lcr_res = lcr_res.rename(columns={"log2_obs_exp": "lcr_log2_residue_count"})
        lcr_res = lcr_res[["method", "feature_value", "rna_class", "lcr_log2_residue_count"]]
    else:
        lcr_res = lcr_seq[["method", "feature_value", "rna_class"]].assign(lcr_log2_residue_count=np.nan)

    prot = protein_cat.loc[(protein_cat["analysis"] == protein_analysis) & (protein_cat["present"] == "present")].copy()
    prot = prot.rename(columns={"log2_obs_exp": "protein_log2", "observed": "protein_observed", "p_fdr": "protein_table_p_fdr"})
    prot = prot[["method", "feature_value", "rna_class", "protein_log2", "protein_observed", "protein_table_p_fdr"]]

    keys = ["method", "feature_value", "rna_class"]
    merged = lcr_seq.merge(lcr_res, on=keys, how="outer").merge(prot, on=keys, how="outer")
    merged.insert(0, "feature_type", feature_type)

    def sign(s: pd.Series) -> pd.Series:
        return np.sign(s)

    merged["lcr_direction"] = sign(merged["lcr_log2_sequence_count"])
    merged["protein_direction"] = sign(merged["protein_log2"])
    both = merged["lcr_direction"].notna() & merged["protein_direction"].notna()
    merged["concordant"] = both & (merged["lcr_direction"] == merged["protein_direction"])
    merged["discordant"] = both & (merged["lcr_direction"] != merged["protein_direction"]) & (merged["lcr_direction"] != 0) & (merged["protein_direction"] != 0)
    abs_min = merged[["lcr_log2_sequence_count", "protein_log2"]].abs().min(axis=1)
    merged["strong_cell"] = merged["concordant"] & (abs_min >= STRONG_LOG2)
    return merged


def concordance(lcr_cat: pd.DataFrame, protein_cat: pd.DataFrame, tables: Path) -> None:
    frames = [concordance_for(ft, lcr_cat, protein_cat) for ft in ["signature", "position"]]
    out = pd.concat(frames, ignore_index=True)
    out = out.sort_values(["feature_type", "method", "feature_value", "rna_class"])
    out.to_csv(tables / "concordance_check.csv", index=False)
    n_conc = int(out["concordant"].sum()); n_disc = int(out["discordant"].sum()); n_strong = int(out["strong_cell"].sum())
    print(f"Concordance cells: {n_conc} concordant, {n_disc} discordant, {n_strong} strong (>=2-fold at both levels).")
    strong = out.loc[out["strong_cell"]].sort_values("protein_log2", key=lambda s: s.abs(), ascending=False)
    if not strong.empty:
        print("\nTop strong concordant cells (|protein log2| ranked):")
        for _, r in strong.head(15).iterrows():
            print(f"  {r['method']:<18} {r['feature_value']:<35} {r['rna_class']:<20} lcr={r['lcr_log2_sequence_count']:+.2f} protein={r['protein_log2']:+.2f}")


def write_readme(outdir: Path) -> None:
    (outdir / "README_deep_read.md").write_text(
        "# Phase 2.4.1 — Deep-read extraction layer\n\n"
        "Reads existing master tables only; no recomputation.\n\n"
        "- `charge_landscape_summary.csv`: median/IQR/mean/SD + KW p/FDR for frac_positive, frac_negative, fcr, ncpr per class x method (LCR level). Tests whether charge differences are distribution-level, not threshold artifacts.\n"
        "- `charge_landscape_median_<metric>.csv`: class x method median pivots; ncpr and frac_positive also rendered as heatmaps.\n"
        "- `significant_coverage_cells.csv`: coverage and property_fraction__* KW cells with FDR < 0.05 (protein level; LCR-level coverage_per_lcr included for reference).\n"
        "- `concordance_check.csv`: per method x feature x class, LCR-level log2 (sequence-count and residue-count bases) beside protein-level log2.\n"
        "  - `concordant`: same direction at both levels.\n"
        "  - `discordant`: opposite non-zero directions (investigate before quoting either level).\n"
        "  - `strong_cell`: concordant AND >= 2-fold (|log2| >= 1) at both levels — the machine-checked validated-cells list.\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default=r"lcr_analyses\quantitative_analysis", help="Results root containing summary/tables")
    args = parser.parse_args()

    root = Path(args.outdir)
    tables = root / "summary" / "deep_read" / "tables"
    figures = root / "summary" / "deep_read" / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    data = load_tables(root)
    if "lcr_cont" in data:
        charge_landscape(data["lcr_cont"], tables, figures)
    if "lcr_cont" in data or "protein_cont" in data:
        significant_coverage(data.get("lcr_cont"), data.get("protein_cont"), tables)
    if "lcr_cat" in data and "protein_cat" in data:
        concordance(data["lcr_cat"], data["protein_cat"], tables)
    write_readme(root / "summary" / "deep_read")
    print(f"\nDeep-read outputs written to: {(root / 'summary' / 'deep_read').resolve()}")


if __name__ == "__main__":
    main()
