#!/usr/bin/env python3
"""Create compact plots from lcr_position_class_analysis.xlsx."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import plotly.express as px

CLASS_ORDER = [
    "domain_intrinsic", "domain_edge", "domain_adjacent",
    "interdomain_linker", "distal_terminal", "unclassified_no_pfam"
]


def save(fig, path):
    fig.write_image(path, scale=2)


def main(input_file, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    length_entropy = pd.read_excel(input_file, sheet_name="length_entropy")
    motifs = pd.read_excel(input_file, sheet_name="motif_prevalence")
    enrichment = pd.read_excel(input_file, sheet_name="superclass_enrichment")

    # 1. Relative abundance of positional classes within each caller.
    classes = length_entropy.copy()
    classes["fraction"] = classes["n_lcrs"] / classes.groupby("method")["n_lcrs"].transform("sum")
    fig = px.bar(
        classes, x="method", y="fraction", color="primary_class",
        category_orders={"primary_class": CLASS_ORDER},
        labels={"method": "Method", "fraction": "Fraction of LCRs", "primary_class": "Position class"},
        title="Position-class distribution within each LCR caller"
    )
    fig.update_layout(barmode="stack")
    fig.update_yaxes(tickformat=".0%")
    save(fig, output_dir / "01_class_distribution.png")

    # 2–3. Heatmaps of robust sequence summaries.
    for column, title, file_name, color_title, fmt in [
        ("median_length", "Median LCR length by position class and caller", "02_median_length.png", "Residues", ".0f"),
        ("median_entropy", "Median LCR entropy by position class and caller", "03_median_entropy.png", "Bits", ".2f"),
    ]:
        matrix = (length_entropy.pivot(index="method", columns="primary_class", values=column)
                  .reindex(columns=CLASS_ORDER))
        fig = px.imshow(
            matrix, text_auto=fmt, aspect="auto", color_continuous_scale="Viridis", # type: ignore
            labels={"x": "Position class", "y": "Method", "color": color_title}, title=title
        )
        save(fig, output_dir / file_name)

    # 4. Motif plot: one figure per motif avoids an unreadable all-motif plot.
    for motif, data in motifs.groupby("motif", sort=False):
        fig = px.bar(
            data, x="primary_class", y="fraction_lcrs_with_motif", color="method", barmode="group",
            category_orders={"primary_class": CLASS_ORDER},
            labels={"primary_class": "Position class", "fraction_lcrs_with_motif": "LCRs with motif", "method": "Method"},
            title=f"{motif} motif prevalence by class and caller"
        )
        fig.update_yaxes(tickformat=".0%")
        safe_name = motif.replace("/", "_").replace(" ", "_") # type: ignore
        save(fig, output_dir / f"motif_{safe_name}.png")

    # 5. Keep statistically practical enrichment rows and plot one heatmap per method.
    ranked = enrichment[(enrichment["observed_n"] >= 10) & (enrichment["class_n"] >= 30)].copy()
    ranked["log2_fold_enrichment"] = np.log2(ranked["fold_enrichment"])
    ranked["abs_log2_fold"] = ranked["log2_fold_enrichment"].abs()
    ranked.to_csv(output_dir / "superclass_enrichment_filtered.csv", index=False)

    for method, data in ranked.groupby("method", sort=False):
        keep_classes = (data.groupby("primary_class")["abs_log2_fold"].max()
                        .sort_values(ascending=False).head(6).index)
        matrix = (data[data["primary_class"].isin(keep_classes)]
                  .pivot(index="rna_primary_class", columns="primary_class", values="log2_fold_enrichment")
                  .reindex(columns=[x for x in CLASS_ORDER if x in keep_classes]))
        fig = px.imshow(
            matrix, text_auto=".2f", aspect="auto", color_continuous_scale="RdBu_r", zmin=-2, zmax=2, # type: ignore
            labels={"x": "Position class", "y": "RNA superclass", "color": "log2 fold enrichment"},
            title=f"RNA-superclass enrichment: {method}"
        )
        save(fig, output_dir / f"enrichment_{method}.png")

    print(f"Wrote plots and filtered enrichment table to: {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path(r"lcr_analyses\domain_function\class_feature_analysis\lcr_position_class_analysis.xlsx"))
    parser.add_argument("--output-dir", type=Path, default=Path(r"lcr_analyses\domain_function\class_feature_analysis\plots"))
    args = parser.parse_args()
    main(args.input, args.output_dir)
