#!/usr/bin/env python3
"""
Plot LCR composition summaries by calling method and RNA superclass.

Inputs
------
lcr_sequence_features.xlsx
Expected sheet: all_LCR_features

Outputs
-------
lcr_figures/
├── figure_1_lcr_counts.png
├── figure_2_median_heatmaps.png
├── figure_3_top_feature_boxplots.png
└── feature_separation_ranking.csv
"""

from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


FEATURES = [
    "lcr_length",
    "dominant_aa_fraction",
    "mmpr",
    "shannon_entropy",
    "aa_richness",
    "dipeptide_dominance",
]

FEATURE_LABELS = {
    "lcr_length": "LCR length",
    "dominant_aa_fraction": "Dominant AA fraction",
    "mmpr": "MMPR",
    "shannon_entropy": "Shannon entropy",
    "aa_richness": "AA richness",
    "dipeptide_dominance": "Dipeptide dominance",
}

MIN_CLASS_SIZE = 10
TOP_N_FEATURES = 3


def save_figure(fig, output_file):
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Create LCR composition figures by method and RNA superclass."
    )
    parser.add_argument("input_file", help="Path to lcr_sequence_features.xlsx")
    parser.add_argument(
        "--sheet",
        default="all_LCR_features",
        help="Excel sheet name (default: all_LCR_features)",
    )
    parser.add_argument(
        "--output-dir",
        default="lcr_figures",
        help="Output directory (default: lcr_figures)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)

    data = pd.read_excel(args.input_file, sheet_name=args.sheet)

    required = ["method", "rna_primary_class"] + FEATURES
    missing = sorted(set(required) - set(data.columns))
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    data = data.copy()
    data["rna_primary_class"] = data["rna_primary_class"].fillna("unknown")

    for feature in FEATURES:
        data[feature] = pd.to_numeric(data[feature], errors="coerce")

    # Retain classes with enough LCRs for stable distributions.
    class_counts = (
        data.groupby(["method", "rna_primary_class"])
        .size()
        .reset_index(name="n")
    )

    valid_groups = class_counts[class_counts["n"] >= MIN_CLASS_SIZE]
    data = data.merge(
        valid_groups[["method", "rna_primary_class"]],
        on=["method", "rna_primary_class"],
        how="inner",
    )

    class_order = (
        data["rna_primary_class"]
        .value_counts()
        .sort_values(ascending=False)
        .index
        .tolist()
    )

    method_order = sorted(data["method"].unique())

    # ------------------------------------------------------------------
    # Figure 1: LCR counts per RNA superclass and calling method
    # ------------------------------------------------------------------
    counts = (
        data.groupby(["method", "rna_primary_class"])
        .size()
        .reset_index(name="LCR count")
    )

    fig, ax = plt.subplots(figsize=(12, 6))
    sns.barplot(
        data=counts,
        x="rna_primary_class",
        y="LCR count",
        hue="method",
        order=class_order,
        hue_order=method_order,
        ax=ax,
    )
    ax.set_title("LCR counts by RNA superclass and caller")
    ax.set_xlabel("RNA superclass")
    ax.set_ylabel("Number of LCRs")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(title="Calling method", bbox_to_anchor=(1.02, 1), loc="upper left")
    save_figure(fig, output_dir / "figure_1_lcr_counts.png")

    # ------------------------------------------------------------------
    # Figure 2: Median-feature heatmap, one panel per calling method.
    # Values are z-scored within each feature across RNA superclasses.
    # ------------------------------------------------------------------
    medians = (
        data.groupby(["method", "rna_primary_class"])[FEATURES]
        .median()
        .reset_index()
    )

    n_methods = len(method_order)
    ncols = min(3, n_methods)
    nrows = int(np.ceil(n_methods / ncols))

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(6 * ncols, 5 * nrows),
        squeeze=False,
    )

    for ax, method in zip(axes.flat, method_order):
        matrix = (
            medians[medians["method"] == method]
            .set_index("rna_primary_class")[FEATURES]
            .reindex(class_order)
            .dropna(how="all")
        )

        z_matrix = matrix.copy()
        for feature in FEATURES:
            std = z_matrix[feature].std()
            if pd.notna(std) and std > 0:
                z_matrix[feature] = (
                    z_matrix[feature] - z_matrix[feature].mean()
                ) / std
            else:
                z_matrix[feature] = 0

        sns.heatmap(
            z_matrix,
            cmap="vlag",
            center=0,
            vmin=-2,
            vmax=2,
            linewidths=0.4,
            linecolor="white",
            cbar_kws={"label": "Median z-score"},
            ax=ax,
        )

        ax.set_title(str(method))
        ax.set_xlabel("Feature")
        ax.set_ylabel("RNA superclass")
        ax.set_xticklabels(
            [FEATURE_LABELS[x] for x in FEATURES],
            rotation=45,
            ha="right",
        )

    for ax in axes.flat[n_methods:]:
        ax.remove()

    fig.suptitle(
        "Median LCR composition by RNA superclass and caller",
        y=1.02,
        fontsize=15,
    )
    fig.tight_layout()
    save_figure(fig, output_dir / "figure_2_median_heatmaps.png")

    # ------------------------------------------------------------------
    # Select features with strongest class separation.
    # Metric: average between-class variance of class medians across methods.
    # ------------------------------------------------------------------
    separation_rows = []

    for feature in FEATURES:
        per_method = []

        for method in method_order:
            method_medians = (
                data[data["method"] == method]
                .groupby("rna_primary_class")[feature]
                .median()
                .dropna()
            )

            if len(method_medians) >= 2:
                per_method.append(method_medians.var())

        separation_rows.append(
            {
                "feature": feature,
                "feature_label": FEATURE_LABELS[feature],
                "median_variance_between_classes": np.mean(per_method),
            }
        )

    ranking = (
        pd.DataFrame(separation_rows)
        .sort_values("median_variance_between_classes", ascending=False)
        .reset_index(drop=True)
    )
    ranking.to_csv(output_dir / "feature_separation_ranking.csv", index=False)

    top_features = ranking["feature"].head(TOP_N_FEATURES).tolist()

    # ------------------------------------------------------------------
    # Figure 3: Boxplots for the three strongest separating features.
    # One row per feature, one column per calling method.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(
        nrows=len(top_features),
        ncols=len(method_order),
        figsize=(5.5 * len(method_order), 4.5 * len(top_features)),
        squeeze=False,
    )

    for row, feature in enumerate(top_features):
        for col, method in enumerate(method_order):
            ax = axes[row, col]

            subset = data[data["method"] == method].copy()

            sns.boxplot(
                data=subset,
                x="rna_primary_class",
                y=feature,
                order=class_order,
                color="#9ecae1",
                showfliers=False,
                ax=ax,
            )

            sns.stripplot(
                data=subset,
                x="rna_primary_class",
                y=feature,
                order=class_order,
                color="black",
                alpha=0.20,
                size=2,
                jitter=0.25,
                ax=ax,
            )

            ax.set_title(str(method))
            ax.set_xlabel("RNA superclass")
            ax.set_ylabel(FEATURE_LABELS[feature])
            ax.tick_params(axis="x", rotation=45)

    fig.suptitle(
        "Top LCR features separating RNA superclasses",
        y=1.01,
        fontsize=15,
    )
    fig.tight_layout()
    save_figure(fig, output_dir / "figure_3_top_feature_boxplots.png")

    print(f"Saved figures to: {output_dir}")
    print("Selected features:")
    for feature in top_features:
        print(f"  - {FEATURE_LABELS[feature]}")


if __name__ == "__main__":
    main()