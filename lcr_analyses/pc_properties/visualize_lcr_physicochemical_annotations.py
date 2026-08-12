#!/usr/bin/env python3
"""Create QC and biological-context visualizations for LCR physicochemical annotations.

Example:
python visualize_lcr_physicochemical_annotations.py \
  --annotations lcr_physicochemical_annotations.xlsx \
  --features lcr_physicochemical_features.xlsx \
  --position lcr_position_classes.xlsx \
  --rna rbp_rna_classification.xlsx \
  --output lcr_physicochemical_figures

Expected sheets:
  annotations workbook: Sheet1
  features workbook: composition, distribution, co_occurrence
  RNA workbook: all_combined

v2.2 update: figures 05/06 show signature prevalence per RNA-target superclass
and per domain-position class, split into one panel per calling method.
Feature/distribution column names aligned with v2 outputs (frac_disorder,
compact/dispersed/insufficient labels).
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_theme(style="whitegrid", context="talk")

ID_COLS = ["protein_id", "method", "start", "end", "length"]
FEATURE_COLUMNS = [
    "frac_polar", "frac_hydrophobic", "frac_strong_hydro", "frac_aromatic",
    "frac_disorder", "frac_positive", "frac_negative", "fcr", "ncpr",
    "frac_GS", "cooc_positive_aromatic",
]
DISTRIBUTION_PROPERTIES = ["polar", "hydrophobic", "aromatic", "disorder", "charged"]


def clean_text(value: object) -> str:
    if pd.isna(value):  # type: ignore
        return ""
    return str(value).strip()


def normalized_key(series: pd.Series) -> pd.Series:
    return series.map(clean_text).str.upper()


def coordinate_key(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    return numeric.map(lambda x: "" if pd.isna(x) else str(int(x)))


def safe_filename(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", label).strip("_")


def require_columns(df: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def add_join_keys(df: pd.DataFrame, start_col: str, end_col: str) -> pd.DataFrame:
    out = df.copy()
    out["_protein_key"] = normalized_key(out["protein_id"])
    out["_method_key"] = normalized_key(out["method"])
    out["_start_key"] = coordinate_key(out[start_col])
    out["_end_key"] = coordinate_key(out[end_col])
    return out


def add_lcr_instance_key(df: pd.DataFrame, join_keys: list[str]) -> pd.DataFrame:
    """Distinguish repeated LCR observations sharing protein/method/coords."""
    out = df.copy()
    out["_lcr_instance"] = out.groupby(join_keys, dropna=False).cumcount()
    return out


def save_figure(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()


def annotation_color_map(order: list[str]) -> dict[str, tuple]:
    palette = sns.color_palette("tab20", max(len(order), 1))
    return {sig: palette[i % len(palette)] for i, sig in enumerate(order)}


def stacked_prevalence(df: pd.DataFrame, group_col: str, output: Path, title: str) -> None:
    plot_df = df.dropna(subset=[group_col, "primary_physicochemical_annotation"]).copy()
    plot_df = plot_df[(plot_df[group_col].astype(str).str.strip() != "") &
                      (plot_df["primary_physicochemical_annotation"].astype(str).str.strip() != "")]
    if plot_df.empty:
        print(f"Skipped {title}: no usable data.")
        return
    counts = pd.crosstab(plot_df[group_col], plot_df["primary_physicochemical_annotation"])
    counts = counts.loc[counts.sum(axis=1).sort_values(ascending=False).index]
    proportions = counts.div(counts.sum(axis=1), axis=0) * 100
    ax = proportions.plot(kind="bar", stacked=True, figsize=(max(9, len(proportions) * 1.15), 7), colormap="tab20")
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("LCRs with primary annotation (%)")
    ax.set_ylim(0, 100)
    ax.tick_params(axis="x", rotation=35)
    ax.legend(title="Primary annotation", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    for i, total in enumerate(counts.sum(axis=1)):
        ax.text(i, 102, f"n={total}", ha="center", va="bottom", fontsize=9)
    save_figure(output)


def stacked_prevalence_per_method(
    df: pd.DataFrame, group_col: str, output: Path, title: str
) -> None:
    """
    One panel per calling method; each panel shows stacked signature
    prevalence (%) across the levels of group_col. Colors are consistent
    across panels (fixed by the global signature-frequency order).
    """
    use = df.dropna(subset=[group_col, "primary_physicochemical_annotation", "method"]).copy()
    use = use[(use[group_col].astype(str).str.strip() != "") &
              (use["primary_physicochemical_annotation"].astype(str).str.strip() != "") &
              (use["method"].astype(str).str.strip() != "")]
    if use.empty:
        print(f"Skipped {title}: no usable data.")
        return

    sig_order = use["primary_physicochemical_annotation"].value_counts().index.tolist()
    colors = annotation_color_map(sig_order)

    methods = use["method"].value_counts().index.tolist()
    n_panels = len(methods)

    fig, axes = plt.subplots(
        1, n_panels,
        figsize=(max(6 * n_panels, 8), 7.5),
        sharey=True,
    )
    if n_panels == 1:
        axes = [axes]

    for ax, method in zip(axes, methods):
        sub = use[use["method"] == method]
        counts = pd.crosstab(sub[group_col], sub["primary_physicochemical_annotation"])
        counts = counts.reindex(columns=sig_order, fill_value=0)
        counts = counts.loc[counts.sum(axis=1).sort_values(ascending=False).index]
        proportions = counts.div(counts.sum(axis=1), axis=0).fillna(0) * 100

        bottom = np.zeros(len(proportions))
        x = np.arange(len(proportions))
        for sig in sig_order:
            vals = proportions[sig].values
            ax.bar(x, vals, bottom=bottom, color=colors[sig], width=0.8)
            bottom += vals

        ax.set_title(f"{method}\n(n={len(sub)})", fontsize=12)
        ax.set_xticks(x)
        ax.set_xticklabels(proportions.index.astype(str), rotation=35, ha="right", fontsize=10)
        ax.set_ylim(0, 100)
        ax.set_xlabel("")
        for i, total in enumerate(counts.sum(axis=1)):
            ax.text(i, 102, f"n={total}", ha="center", va="bottom", fontsize=8)

    axes[0].set_ylabel("LCRs with primary annotation (%)")

    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[sig]) for sig in sig_order]  # type: ignore
    fig.legend(
        handles, sig_order,
        title="Primary annotation",
        loc="lower center",
        ncol=min(4, len(sig_order)),
        bbox_to_anchor=(0.5, -0.06),
        frameon=False,
        fontsize=9,
    )
    fig.suptitle(title, y=1.02)
    save_figure(output)


def lcr_counts_by_method(df: pd.DataFrame, output: Path) -> None:
    counts = df["method"].fillna("Unknown").value_counts().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(max(7, len(counts) * 1.05), 5.5))
    bars = ax.bar(counts.index.astype(str), counts.values, color=sns.color_palette("deep", len(counts)))  # type: ignore
    ax.set_title("LCR observations by calling method")
    ax.set_xlabel("Calling method")
    ax.set_ylabel("Number of LCRs")
    ax.tick_params(axis="x", rotation=35)
    for bar, value in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, value, str(value), ha="center", va="bottom", fontsize=10)
    save_figure(output)


def feature_heatmap(df: pd.DataFrame, output: Path) -> None:
    cols = [c for c in FEATURE_COLUMNS if c in df.columns]
    use = df.dropna(subset=["primary_physicochemical_annotation"]).copy()
    use = use[use["primary_physicochemical_annotation"].astype(str).str.strip() != ""]
    if use.empty or not cols:
        print("Skipped feature-profile heatmap: no usable data.")
        return
    medians = use.groupby("primary_physicochemical_annotation", observed=True)[cols].median()
    medians = medians.loc[use["primary_physicochemical_annotation"].value_counts().index.intersection(medians.index)]
    scaled = medians.copy()
    for col in scaled.columns:
        sd = scaled[col].std(ddof=0)
        scaled[col] = 0 if pd.isna(sd) or sd == 0 else (scaled[col] - scaled[col].mean()) / sd
    plt.figure(figsize=(max(12, len(cols) * 0.75), max(4.5, len(scaled) * 0.6)))
    sns.heatmap(scaled, cmap="vlag", center=0, linewidths=0.5, linecolor="white", cbar_kws={"label": "Standardized median"})
    plt.title("Physicochemical feature profile by primary annotation")
    plt.xlabel("Feature")
    plt.ylabel("Primary annotation")
    save_figure(output)


def feature_dotplot(df: pd.DataFrame, output: Path) -> None:
    cols = [c for c in ["frac_polar", "frac_hydrophobic", "frac_aromatic", "frac_disorder", "fcr", "ncpr"] if c in df.columns]
    use = df.dropna(subset=["primary_physicochemical_annotation"]).copy()
    use = use[use["primary_physicochemical_annotation"].astype(str).str.strip() != ""]
    if use.empty or not cols:
        print("Skipped selected-feature dot plot: no usable data.")
        return
    long = use.melt(id_vars="primary_physicochemical_annotation", value_vars=cols, var_name="feature", value_name="value")
    summary = long.groupby(["primary_physicochemical_annotation", "feature"], observed=True)["value"].agg(
        median="median", q1=lambda x: x.quantile(0.25), q3=lambda x: x.quantile(0.75)
    ).reset_index()
    order = use["primary_physicochemical_annotation"].value_counts().index.tolist()
    fig, axes = plt.subplots(1, len(cols), figsize=(max(15, len(cols) * 2.5), max(5.5, len(order) * 0.5)), sharey=True)
    if len(cols) == 1:
        axes = [axes]
    for ax, feature in zip(axes, cols):
        subset = summary[summary["feature"] == feature].set_index("primary_physicochemical_annotation").reindex(order).reset_index()
        y = np.arange(len(subset))
        ax.errorbar(subset["median"], y, xerr=[subset["median"] - subset["q1"], subset["q3"] - subset["median"]], fmt="o", color="#2b6cb0", capsize=3)
        ax.set_title(feature.replace("_", "\n"))
        ax.set_xlabel("Median (IQR)")
        ax.set_yticks(y)
        ax.set_yticklabels(order if ax is axes[0] else [])
        ax.axvline(0, color="grey", lw=0.7, zorder=0)
    fig.suptitle("Selected composition and charge features by annotation", y=1.02)
    save_figure(output)


def distribution_label_plot(df: pd.DataFrame, output: Path) -> None:
    label_cols = [f"{prop}_label" for prop in DISTRIBUTION_PROPERTIES]
    available = [col for col in label_cols if col in df.columns]
    if not available:
        print("Skipped distribution-label plot: no distribution labels found.")
        return

    annotation_col = "primary_physicochemical_annotation"
    use = df.dropna(subset=[annotation_col]).copy()
    use = use[use[annotation_col].astype(str).str.strip() != ""]
    if use.empty:
        print("Skipped distribution-label plot: no usable annotation data.")
        return

    annotation_order = use[annotation_col].value_counts().index.tolist()

    # v2 labels are lowercase
    label_order = ["compact", "dispersed", "insufficient", "Unassigned"]
    palette = {
        "compact": "#d73027",
        "dispersed": "#4575b4",
        "insufficient": "#bdbdbd",
        "Unassigned": "#f0f0f0",
    }

    fig, axes = plt.subplots(
        1, len(available),
        figsize=(max(18, len(available) * 4.2), max(6, len(annotation_order) * 0.45)),
        sharey=True,
    )
    if len(available) == 1:
        axes = [axes]

    for ax, label_col in zip(axes, available):
        property_name = label_col.replace("_label", "").replace("_", " ").title()

        plot_df = use[[annotation_col, label_col]].copy()
        plot_df[label_col] = plot_df[label_col].fillna("Unassigned")

        counts = pd.crosstab(plot_df[annotation_col], plot_df[label_col]).reindex(annotation_order, fill_value=0)
        counts = counts.reindex(columns=label_order, fill_value=0)
        proportions = counts.div(counts.sum(axis=1), axis=0).fillna(0) * 100

        proportions.plot(
            kind="barh", stacked=True, ax=ax,
            color=[palette[label] for label in proportions.columns],
            width=0.8, legend=False,
        )
        ax.set_title(f"{property_name} distribution")
        ax.set_xlabel("LCRs (%)")
        ax.set_xlim(0, 100)
        ax.set_ylabel("")
        ax.grid(axis="x", alpha=0.3)

        if ax is axes[0]:
            ax.set_yticklabels(annotation_order)
        else:
            ax.set_yticklabels([])

    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[label]) for label in label_order]  # type: ignore
    fig.legend(
        handles, label_order,
        title="Distribution label",
        loc="lower center",
        ncol=len(label_order),
        bbox_to_anchor=(0.5, -0.08),
        frameon=False,
    )
    fig.suptitle("Property-distribution labels by primary physicochemical annotation", y=1.02)
    save_figure(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--annotations", default="lcr_analyses/pc_properties/lcr_annotations.xlsx", help="Physicochemical annotation workbook.")
    parser.add_argument("--features", default="lcr_analyses/pc_properties/lcr_features.xlsx", help="Physicochemical feature workbook.")
    parser.add_argument("--position", default=r"lcr_analyses\domain_function\lcr_position_classes.csv", help="Domain-position classification workbook or CSV.")
    parser.add_argument("--rna", default=r"rbp_superclasses\rbp_rna_classification.xlsx", help="RNA target superclass workbook.")
    parser.add_argument("--output", default="lcr_analyses/pc_properties/plots_v2", help="Output directory.")
    args = parser.parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    annotations = pd.read_excel(args.annotations, sheet_name="Sheet1")
    composition = pd.read_excel(args.features, sheet_name="composition")
    distribution = pd.read_excel(args.features, sheet_name="distribution")
    cooc = pd.read_excel(args.features, sheet_name="co_occurrence")
    position = pd.read_csv(args.position) if str(args.position).lower().endswith(".csv") else pd.read_excel(args.position)
    rna = pd.read_excel(args.rna, sheet_name="all_combined")

    require_columns(annotations, ID_COLS + ["sequence", "primary_physicochemical_annotation"], "annotations sheet")
    require_columns(composition, ID_COLS, "composition sheet")
    require_columns(distribution, ID_COLS, "distribution sheet")
    require_columns(cooc, ID_COLS, "co_occurrence sheet")
    require_columns(position, ["protein_id", "method", "lcr_start", "lcr_end", "protein_accession", "primary_class"], "position classification")
    require_columns(rna, ["uniprot_accession", "rna_primary_class"], "RNA classification")

    annotations = add_join_keys(annotations, "start", "end")
    composition = add_join_keys(composition, "start", "end")
    distribution = add_join_keys(distribution, "start", "end")
    cooc = add_join_keys(cooc, "start", "end")
    position = add_join_keys(position, "lcr_start", "lcr_end")

    join_keys = ["_protein_key", "_method_key", "_start_key", "_end_key"]

    # Preserve repeated LCR observations without creating many-to-many joins.
    annotations = add_lcr_instance_key(annotations, join_keys)
    composition = add_lcr_instance_key(composition, join_keys)
    distribution = add_lcr_instance_key(distribution, join_keys)
    cooc = add_lcr_instance_key(cooc, join_keys)

    feature_join_keys = join_keys + ["_lcr_instance"]

    position_keep = join_keys + [c for c in ["protein_accession", "primary_class", "class_labels", "max_lcr_domain_fraction", "nearest_domain_distance", "nearest_pfam_accession", "nearest_pfam_name", "flanking_left_pfam", "flanking_right_pfam"] if c in position.columns]
    position_keep = position[position_keep].drop_duplicates(join_keys, keep="first")

    feature_cols = [c for c in composition.columns if c.startswith("frac_") or c in ["fcr", "ncpr"]]
    comp_keep = feature_join_keys + feature_cols
    cooc_cols = [c for c in cooc.columns if c.startswith("cooc_")]
    cooc_keep = feature_join_keys + cooc_cols
    distribution_label_cols = [c for c in distribution.columns if c.endswith("_label")]
    distribution_keep = feature_join_keys + distribution_label_cols

    merged = annotations.merge(composition[comp_keep], on=feature_join_keys, how="left", validate="one_to_one")
    merged = merged.merge(distribution[distribution_keep], on=feature_join_keys, how="left", validate="one_to_one")
    merged = merged.merge(cooc[cooc_keep], on=feature_join_keys, how="left", validate="one_to_one")
    merged = merged.merge(position_keep, on=join_keys, how="left", validate="many_to_one")

    rna = rna.copy()
    rna["_accession_key"] = normalized_key(rna["uniprot_accession"])
    rna_keep = [c for c in ["_accession_key", "rna_primary_class", "rna_secondary_class", "confidence", "evidence_source", "evidence_matched"] if c in rna.columns]
    rna_keep = rna[rna_keep].drop_duplicates("_accession_key", keep="first")
    merged["_accession_key"] = normalized_key(merged["protein_accession"])
    merged = merged.merge(rna_keep, on="_accession_key", how="left", validate="many_to_one")

    merged["domain_position_class"] = merged["primary_class"].fillna("Unassigned")
    merged["rna_target_superclass"] = merged["rna_primary_class"].fillna("Unassigned")

    merged.drop(columns=[c for c in merged.columns if c.startswith("_")], errors="ignore").to_excel(outdir / "lcr_physicochemical_visualization_data.xlsx", index=False)

    stacked_prevalence(merged, "method", outdir / "01_annotations_by_calling_method", "Primary physicochemical annotations by calling method")
    stacked_prevalence(merged, "rna_target_superclass", outdir / "02_annotations_by_rna_target_superclass", "Primary physicochemical annotations by RNA-target superclass")
    stacked_prevalence(merged, "domain_position_class", outdir / "03_annotations_by_domain_position", "Primary physicochemical annotations by domain-position class")
    feature_heatmap(merged, outdir / "04_annotation_feature_profile_heatmap")

    # NEW: per-calling-method panels
    stacked_prevalence_per_method(
        merged, "rna_target_superclass",
        outdir / "05_annotations_by_rna_target_per_method",
        "Primary physicochemical annotations by RNA-target superclass, per calling method",
    )
    stacked_prevalence_per_method(
        merged, "domain_position_class",
        outdir / "06_annotations_by_domain_position_per_method",
        "Primary physicochemical annotations by domain-position class, per calling method",
    )

    report = [
        f"Total LCR observations: {len(merged)}",
        f"Unique proteins: {merged['protein_id'].nunique(dropna=True)}",
        f"With domain-position match: {merged['primary_class'].notna().sum()} ({merged['primary_class'].notna().mean():.1%})",
        f"With RNA superclass match: {merged['rna_primary_class'].notna().sum()} ({merged['rna_primary_class'].notna().mean():.1%})",
        "",
        "Primary annotation counts:",
        merged['primary_physicochemical_annotation'].fillna('Unassigned').value_counts().to_string(),
    ]
    (outdir / "run_summary.txt").write_text("\n".join(report), encoding="utf-8")
    print(f"Completed. Outputs written to: {outdir.resolve()}")


if __name__ == "__main__":
    main()
