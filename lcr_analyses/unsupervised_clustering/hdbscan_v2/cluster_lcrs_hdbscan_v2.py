#!/usr/bin/env python3
"""
Cluster LCRs independently within each RNA superclass using HDBSCAN.

Phase 2.1.2 -- expanded feature set (minimal-change version of
cluster_lcrs_hdbscan.py: only the feature list is updated to the v2
feature table produced by create_lcr_features_v2.py).

Example:
python cluster_lcrs_hdbscan_v2.py \
--input lcr_sequence_features_v2.xlsx \
--method fLPS_strict \
--min-cluster-size 10
"""

import argparse
from pathlib import Path

import hdbscan
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

FEATURES = [
    # entropy / periodicity (normalized forms)
    "norm_entropy",
    "mmpr_norm",
    "best_period",
    "dipeptide_dominance",
    "dominant_aa_fraction",
    # physicochemical composition
    "fcr",
    "ncpr",
    "frac_polar",
    "frac_strong_hydro",
    "frac_aromatic",
    "frac_disorder",
    # flank/end position
    "rel_start",
    "rel_end",
    "lcr_coverage",
    "terminus_dist_pct",
    # size
    "log_lcr_length",
]

MIN_INPUT_SIZE = 30


def read_feature_table(path: Path) -> pd.DataFrame:
    """Read the Excel feature workbook or a CSV feature table."""
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)

    return pd.read_excel(path, sheet_name="all_LCR_features")


def make_plot(data: pd.DataFrame, output_path: Path, title: str) -> None:
    """Create a PCA scatter plot coloured by RNA superclass."""
    fig, ax = plt.subplots(figsize=(8, 6), dpi=180)

    for superclass, subset in data.groupby("rna_primary_class", dropna=False):
        label = "unclassified" if pd.isna(superclass) else str(superclass)

        ax.scatter(
            subset["pca_1"],
            subset["pca_2"],
            s=18,
            alpha=0.7,
            label=label,
            edgecolors="none",
        )

    ax.set_title(title)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(
        title="RNA superclass",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Feature-table XLSX or CSV file.",
    )

    parser.add_argument(
        "--method",
        required=True,
        help="LCR detection method to analyse, e.g. CAST.",
    )

    parser.add_argument(
        "--min-cluster-size",
        default=5,
        type=int,
        help="Minimum number of LCRs required for an HDBSCAN cluster.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(r"lcr_analyses\unsupervised_clustering\hdbscan_v2\output"),
        help="Directory for CSV and plot outputs.",
    )

    args = parser.parse_args()

    if args.min_cluster_size < 2:
        raise ValueError("--min-cluster-size must be at least 2.")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    plots_dir = args.output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    table = read_feature_table(args.input)

    method_data = table[table["method"] == args.method].copy()

    if method_data.empty:
        raise ValueError(f"Method not found: {args.method}")

    method_data = method_data.dropna(
        subset=[
            "rna_primary_class",
            "lcr_length",
            "norm_entropy",
            "mmpr_norm",
            "best_period",
            "dipeptide_dominance",
            "dominant_aa_fraction",
            "fcr",
            "ncpr",
            "frac_polar",
            "frac_strong_hydro",
            "frac_aromatic",
            "frac_disorder",
            "rel_start",
            "rel_end",
            "lcr_coverage",
            "terminus_dist_pct",
        ]
    )

    method_data = method_data[method_data["lcr_length"] > 0].copy()
    method_data["log_lcr_length"] = np.log(method_data["lcr_length"])

    result_rows = []
    summary_rows = []

    for superclass, subset in method_data.groupby(
        "rna_primary_class",
        sort=True,
    ):
        n_input = len(subset)

        if n_input < MIN_INPUT_SIZE:
            summary_rows.append(
                {
                    "method": args.method,
                    "rna_primary_class": superclass,
                    "n_input": n_input,
                    "status": "skipped_n_lt_30",
                    "min_cluster_size": args.min_cluster_size,
                    "n_clusters": np.nan,
                    "noise_fraction": np.nan,
                }
            )
            continue

        X = subset[FEATURES]
        X_scaled = StandardScaler().fit_transform(X)

        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=args.min_cluster_size,
            min_samples=2,
            metric="euclidean",
            cluster_selection_method="eom",
        )

        labels = clusterer.fit_predict(X_scaled)

        clustered = subset[
            [
                "protein_id",
                "uniprot_accession",
                "lcr_start",
                "lcr_end",
                "sequence",
                "rna_primary_class",
            ]
        ].copy()

        clustered.insert(0, "method", args.method)
        clustered["cluster_id"] = labels

        result_rows.append(clustered)

        n_clusters = len(set(labels)) - int(-1 in labels)

        summary_rows.append(
            {
                "method": args.method,
                "rna_primary_class": superclass,
                "n_input": n_input,
                "status": "clustered",
                "min_cluster_size": args.min_cluster_size,
                "n_clusters": n_clusters,
                "noise_fraction": round(float(np.mean(labels == -1)), 4),
            }
        )

    results = (
        pd.concat(result_rows, ignore_index=True)
        if result_rows
        else pd.DataFrame()
    )

    summary = pd.DataFrame(summary_rows)

    results.to_csv(
        args.output_dir / f"{args.method}_hdbscan_clusters.csv",
        index=False,
    )

    summary.to_csv(
        args.output_dir / f"{args.method}_hdbscan_summary.csv",
        index=False,
    )

    eligible_classes = summary.loc[
        summary["status"] == "clustered",
        "rna_primary_class",
    ]

    plot_data = method_data[
        method_data["rna_primary_class"].isin(eligible_classes)
    ].copy()

    plot_data = method_data[
        method_data["rna_primary_class"].isin(eligible_classes)
    ].copy()

    # --- Exclude unknown class ---
    EXCLUDED_FROM_PLOTS = {"unknown"}
    plot_data = plot_data[
        ~plot_data["rna_primary_class"]
        .fillna("")
        .str.lower()
        .isin(EXCLUDED_FROM_PLOTS)
    ]
    # ----------------

    if len(plot_data) >= 2:
        X_plot = StandardScaler().fit_transform(plot_data[FEATURES])

        coordinates = PCA(
            n_components=2,
            random_state=42,
        ).fit_transform(X_plot)

        plot_data["pca_1"] = coordinates[:, 0]
        plot_data["pca_2"] = coordinates[:, 1]

        make_plot(
            plot_data,
            plots_dir / f"{args.method}_all_superclasses_pca.png",
            f"{args.method}: LCR feature space by RNA superclass",
        )

        for superclass, subset in plot_data.groupby(
            "rna_primary_class",
            sort=True,
        ):
            safe_name = (
                str(superclass)
                .replace("/", "_")
                .replace(" ", "_")
            )

            make_plot(
                subset,
                plots_dir / f"{args.method}_{safe_name}_pca.png",
                f"{args.method}: {superclass} LCR feature space",
            )


if __name__ == "__main__":
    main()
