#!/usr/bin/env python3
"""Plot HDBSCAN LCR results in a 2D UMAP feature space.

Example:
python plot_hdbscan_umap.py \
    --clusters AlcoR_hdbscan_clusters.csv \
    --features lcr_sequence_features.xlsx \
    --output-dir umap_AlcoR
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap.umap_ as umap
from sklearn.preprocessing import StandardScaler

FEATURES = [
    "log_lcr_length",
    "dominant_aa_fraction",
    "mmpr",
    "shannon_entropy",
    "aa_richness",
    "dipeptide_dominance",
]
KEYS = [
    "method",
    "protein_id",
    "uniprot_accession",
    "lcr_start",
    "lcr_end",
    "sequence",
    "rna_primary_class",
]


def read_features(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    return pd.read_excel(path, sheet_name="all_LCR_features")


def plot_clusters(data: pd.DataFrame, path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 6), dpi=180)
    noise = data[data["cluster_id"] == -1]
    clustered = data[data["cluster_id"] != -1]

    if not noise.empty:
        ax.scatter(noise["umap_1"], noise["umap_2"], s=12, c="lightgrey",
                   alpha=0.55, label="noise (-1)", edgecolors="none")

    if not clustered.empty:
        labels = clustered["cluster_id"].astype(int)
        scatter = ax.scatter(clustered["umap_1"], clustered["umap_2"], s=16,
                             c=labels, cmap="tab20", alpha=0.8, edgecolors="none")
        fig.colorbar(scatter, ax=ax, label="HDBSCAN cluster ID")

    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    if not noise.empty:
        ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clusters", required=True, type=Path,
                        help="Output CSV from cluster_lcrs_hdbscan.py")
    parser.add_argument("--features", required=True, type=Path,
                        help="Original XLSX or CSV feature table")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--min-dist", type=float, default=0.1)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = args.output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    clusters = pd.read_csv(args.clusters)
    features = read_features(args.features)

    required_cluster_columns = [*KEYS, "cluster_id"]
    missing = set(required_cluster_columns).difference(clusters.columns)
    if missing:
        raise ValueError(f"Missing columns in cluster CSV: {sorted(missing)}")

    required_feature_columns = [*KEYS, "lcr_length", *FEATURES[1:]]
    missing = set(required_feature_columns).difference(features.columns)
    if missing:
        raise ValueError(f"Missing columns in feature table: {sorted(missing)}")

    feature_columns = [*KEYS, "lcr_length", *FEATURES[1:]]

    clusters_for_merge = (
        clusters
        .drop_duplicates(subset=KEYS, keep="first")
        .copy()
    )

    features_for_merge = (
        features[feature_columns]
        .drop_duplicates(subset=KEYS, keep="first")
        .copy()
    )

    data = clusters_for_merge.merge(
        features_for_merge,
        on=KEYS,
        how="left",
        validate="one_to_one",
    )
    data = data.dropna(subset=["lcr_length", *FEATURES[1:]]).copy()
    data = data[data["lcr_length"] > 0].copy()
    data["log_lcr_length"] = np.log(data["lcr_length"])

    if len(data) < 3:
        raise ValueError("Fewer than 3 LCRs remain after merging feature values.")

    X = StandardScaler().fit_transform(data[FEATURES])
    embedding = umap.UMAP(n_components=2, n_neighbors=args.n_neighbors,
                          min_dist=args.min_dist, metric="euclidean",
                          random_state=42).fit_transform(X)
    data["umap_1"] = embedding[:, 0] # type: ignore
    data["umap_2"] = embedding[:, 1] # type: ignore

    method = data["method"].iloc[0]
    plot_clusters(data, plots_dir / f"{method}_all_hdbscan_umap.png",
                  f"{method}: HDBSCAN clusters in UMAP feature space")

    for superclass, subset in data.groupby("rna_primary_class", sort=True):
        safe_name = str(superclass).replace("/", "_").replace(" ", "_")
        plot_clusters(subset, plots_dir / f"{method}_{safe_name}_hdbscan_umap.png",
                      f"{method}: {superclass} HDBSCAN clusters in UMAP")

    data.to_csv(args.output_dir / f"{method}_hdbscan_umap_coordinates.csv", index=False)


if __name__ == "__main__":
    main()
