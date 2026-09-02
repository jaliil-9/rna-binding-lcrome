#!/usr/bin/env python3
"""Plot carrier PAM k=5 clusters in a 2D classical-MDS embedding.

The script reconstructs the carrier feature matrix, recomputes the same
three-block Gower distance used for PAM, derives a deterministic 2D classical
metric-MDS / PCoA embedding, and highlights recalculated PAM medoids.
"""

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
METHODS = ["AlcoR", "CAST", "FLPS", "LCRFinder", "SEG", "SEG_intermediate"]
CONFIG = {
    "method": METHODS[5],
    "cluster_file": "lcr_analyses/clustering/carriers/" + METHODS[5] + "/protein_clusters.csv",
    "annotations_file": "lcr_analyses/pc_properties/lcr_annotations.xlsx",
    "features_file": "lcr_analyses/pc_properties/lcr_features.xlsx",
    "output_dir": "lcr_analyses/clustering/visualization/scatter_plots/carriers/" + METHODS[5],
    "cluster_column": "cluster_pam_k5",
    "dpi": 250,
    "point_size": 20,
    "point_alpha": 0.55,
}

# Annotation and composition tables may have different pipeline labels for the
# same segmentation. Keys/values are normalized to lowercase.
COMPOSITION_METHOD_MAP = {
    "flps": "flps_strict",
    "seg": "seg_strict",
    "seg_intermediate": "seg_intermediate",
}

POSITION_ORDER = [
    "domain_intrinsic",
    "domain_edge",
    "domain_adjacent",
    "interdomain_linker",
    "distal_terminal",
]

CARRIER_COMPOSITION_METRICS = [
    "frac_polar",
    "frac_hydrophobic",
    "frac_strong_hydro",
    "frac_aromatic",
    "frac_disorder",
    "frac_positive",
    "frac_negative",
    "fcr",
    "ncpr",
    "frac_gs",
]

CLUSTER_COLOURS = [
    "#4c78a8",  # blue
    "#f58518",  # orange
    "#54a24b",  # green
    "#e45756",  # red
    "#b279a2",  # purple
    "#72b7b2",  # teal
    "#ff9da6",  # pink
    "#9d755d",  # brown
    "#bab0ac",  # grey
    "#edc949",  # yellow
]


# -----------------------------------------------------------------------------
# Input and standardization
# -----------------------------------------------------------------------------
def read_table(path, sheet_name=0):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input not found: {path}")
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name=sheet_name)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"Unsupported file type: {path}")


def normalise_text(value):
    if pd.isna(value):
        return np.nan
    value = str(value).strip().lower()
    value = re.sub(r"[ /-]+", "_", value)
    return re.sub(r"_+", "_", value)


def standardise_columns(df, name):
    df = df.copy()
    df.columns = [str(c).replace("\ufeff", "").strip().lower() for c in df.columns]
    aliases = {
        "uniprot_accession": "protein_id",
        "proteinid": "protein_id",
        "uniprot_length": "protein_length",
        "sourcemethod": "source_method",
    }
    df = df.rename(columns={c: aliases[c] for c in df.columns if c in aliases})
    if "protein_id" not in df.columns:
        raise ValueError(f"{name}: no protein ID column. Found: {df.columns.tolist()}")

    # Annotation IDs may contain UniProt plus additional provenance fields.
    df["protein_id"] = (
        df["protein_id"].astype(str).str.strip().str.split("|", n=1).str[0].str.upper()
    )
    for col in ("source_method", "domain_position_class", "primary_physicochemical_annotation"):
        if col in df.columns:
            df[col] = df[col].map(normalise_text)
    return df


def require_columns(df, columns, name):
    missing = [x for x in columns if x not in df.columns]
    if missing:
        raise ValueError(f"{name}: missing columns {missing}")


# -----------------------------------------------------------------------------
# Carrier feature reconstruction
# -----------------------------------------------------------------------------
def weighted_protein_properties(lcr, metrics):
    result = pd.DataFrame(index=lcr["protein_id"].drop_duplicates())
    for metric in metrics:
        valid = lcr[["protein_id", "length", metric]].dropna().copy()
        numerator = (
            valid[metric].astype(float).mul(valid["length"].astype(float))
            .groupby(valid["protein_id"], observed=True).sum()
        )
        denominator = valid["length"].astype(float).groupby(valid["protein_id"], observed=True).sum()
        result[metric] = numerator.div(denominator)
    return result


def build_carrier_features(clusters, annotations, composition, method):
    method_key = normalise_text(method)
    composition_method = COMPOSITION_METHOD_MAP.get(method_key, method_key) # type: ignore
    protein_ids = pd.Index(clusters["protein_id"].unique(), name="protein_id")

    annotations = annotations.loc[
        (annotations["source_method"] == method_key)
        & annotations["protein_id"].isin(protein_ids)
    ].copy()
    composition = composition.loc[
        (composition["source_method"] == composition_method)
        & composition["protein_id"].isin(protein_ids)
    ].copy()

    # Harmonize only the provenance label used in the coordinate join.
    composition["source_method"] = method_key

    print(f"Cluster proteins: {len(protein_ids):,}")
    print(f"Matched annotation proteins: {annotations['protein_id'].nunique():,}")
    print(f"Matched annotation LCRs: {len(annotations):,}")
    print(f"Matched composition proteins: {composition['protein_id'].nunique():,}")
    print(f"Matched composition LCRs: {len(composition):,}")

    key = ["protein_id", "source_method", "start", "end"]
    require_columns(
        annotations,
        key + ["length", "domain_position_class", "primary_physicochemical_annotation"],
        "annotations",
    )
    require_columns(composition, key + CARRIER_COMPOSITION_METRICS, "composition")
    if annotations.duplicated(key).any() or composition.duplicated(key).any():
        raise ValueError("Duplicate LCR coordinate keys found in annotations or composition.")

    lcr = annotations.merge(
        composition[key + CARRIER_COMPOSITION_METRICS],
        on=key,
        how="left",
        validate="one_to_one",
        indicator="_composition_merge",
    )
    missing = lcr.loc[lcr["_composition_merge"] != "both"]
    if not missing.empty:
        print(missing[key + ["length"]].to_string(index=False))
        raise ValueError(f"{len(missing)} LCRs lack matching composition values.")
    lcr = lcr.drop(columns="_composition_merge")
    lcr["length"] = pd.to_numeric(lcr["length"], errors="raise")

    base = clusters[["protein_id", "protein_length"]].drop_duplicates().set_index("protein_id")
    if base.index.duplicated().any() or (base["protein_length"] <= 0).any():
        raise ValueError("Protein IDs must be unique and protein lengths must be positive.")

    architecture = lcr.groupby("protein_id", observed=True).agg(
        n_lcr=("length", "size"),
        lcr_residues=("length", "sum"),
        mean_lcr_length=("length", "mean"),
    )
    protein = base.join(architecture)
    if protein[["n_lcr", "lcr_residues", "mean_lcr_length"]].isna().any().any():
        raise ValueError("Carrier cluster file includes a protein without matched LCR annotations.")
    protein["n_lcr"] = protein["n_lcr"].astype(int)
    protein["coverage"] = protein["lcr_residues"] / protein["protein_length"]

    observed_positions = [x for x in POSITION_ORDER if x in set(lcr["domain_position_class"].dropna())]
    observed_signatures = sorted(set(lcr["primary_physicochemical_annotation"].dropna()))
    for position in observed_positions:
        present = lcr.loc[lcr["domain_position_class"] == position, "protein_id"].unique()
        protein[f"pos__{position}"] = protein.index.isin(present).astype(int)
    for signature in observed_signatures:
        present = lcr.loc[lcr["primary_physicochemical_annotation"] == signature, "protein_id"].unique()
        protein[f"sig__{signature}"] = protein.index.isin(present).astype(int)

    properties = weighted_protein_properties(lcr, CARRIER_COMPOSITION_METRICS)
    protein = protein.join(properties)
    if protein[CARRIER_COMPOSITION_METRICS].isna().any().any():
        raise ValueError("At least one carrier protein has an undefined weighted composition property.")

    return clusters.merge(protein.reset_index(), on=["protein_id", "protein_length"], how="inner", validate="one_to_one")


# -----------------------------------------------------------------------------
# Gower distance, medoids, and classical MDS
# -----------------------------------------------------------------------------
def gower_mixed_distance(df, binary_cols, categorical_cols, continuous_cols):
    """Exact equal-weight three-block mixed distance used for PAM."""
    components = []
    if binary_cols:
        components.append(pairwise_distances(df[binary_cols], metric="hamming"))
    if categorical_cols:
        Xc = df[categorical_cols].astype(str).to_numpy()
        components.append(np.mean(Xc[:, None, :] != Xc[None, :, :], axis=2, dtype=float))
    if continuous_cols:
        X = df[continuous_cols].astype(float).to_numpy()
        lo, hi = np.nanmin(X, axis=0), np.nanmax(X, axis=0)
        span = hi - lo
        keep = span > 0
        if keep.any():
            X = (X[:, keep] - lo[keep]) / span[keep]
            components.append(pairwise_distances(X, metric="manhattan") / X.shape[1])
    if not components:
        raise ValueError("No variable clustering features remain.")
    distance = np.mean(components, axis=0)
    np.fill_diagonal(distance, 0.0)
    return distance


def feature_blocks(features):
    binary = [c for c in features.columns if c.startswith(("pos__", "sig__"))]
    continuous = ["n_lcr", "coverage", "mean_lcr_length", *CARRIER_COMPOSITION_METRICS]
    binary = [c for c in binary if features[c].nunique(dropna=False) > 1]
    continuous = [c for c in continuous if features[c].nunique(dropna=False) > 1]
    return binary, [], continuous


def pam_medoids(distance, cluster_labels):
    medoid_indices = []
    labels = np.asarray(cluster_labels)
    for cluster in sorted(pd.unique(labels), key=lambda x: int(x) if str(x).isdigit() else str(x)):
        indices = np.flatnonzero(labels == cluster)
        within = distance[np.ix_(indices, indices)]
        medoid_indices.append(indices[np.argmin(within.sum(axis=1))])
    return np.asarray(medoid_indices, dtype=int)


def classical_mds(distance, n_components=2):
    """Deterministic classical metric MDS, also called PCoA."""
    n = distance.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ np.square(distance) @ J
    eigenvalues, eigenvectors = np.linalg.eigh(B)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    positive = np.clip(eigenvalues[:n_components], 0, None)
    coordinates = eigenvectors[:, :n_components] * np.sqrt(positive)
    positive_total = eigenvalues[eigenvalues > 0].sum()
    explained = positive / positive_total if positive_total > 0 else np.full(n_components, np.nan)
    return coordinates, eigenvalues, explained


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def cluster_sort_key(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return str(value)


def plot_embedding(df, medoid_indices, explained, config):
    cluster_values = sorted(df["cluster_id"].unique(), key=cluster_sort_key)
    fig, ax = plt.subplots(figsize=(8.2, 7.2), constrained_layout=True)

    for i, cluster in enumerate(cluster_values):
        subset = df.loc[df["cluster_id"] == cluster]
        colour = CLUSTER_COLOURS[i % len(CLUSTER_COLOURS)]
        ax.scatter(
            subset["mds_1"], subset["mds_2"],
            s=config["point_size"], c=colour, alpha=config["point_alpha"],
            linewidths=0, label=f"Cluster {cluster} (n={len(subset)})", rasterized=True,
        )

    medoids = df.iloc[medoid_indices]
    ax.scatter(
        medoids["mds_1"], medoids["mds_2"],
        marker="*", s=240, c="white", edgecolors="black", linewidths=1.0,
        zorder=5, label="PAM medoid",
    )
    for _, row in medoids.iterrows():
        ax.annotate(
            f"C{row['cluster_id']}", (row["mds_1"], row["mds_2"]),
            xytext=(6, 6), textcoords="offset points", fontsize=9, fontweight="bold",
            bbox={"boxstyle": "round,pad=0.16", "fc": "white", "ec": "none", "alpha": 0.75},
            zorder=6,
        )

    ax.axhline(0, color="#dddddd", lw=0.7, zorder=0)
    ax.axvline(0, color="#dddddd", lw=0.7, zorder=0)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel(f"Classical MDS 1 ({100 * explained[0]:.1f}% positive-eigenvalue variance)")
    ax.set_ylabel(f"Classical MDS 2 ({100 * explained[1]:.1f}% positive-eigenvalue variance)")
    ax.set_title(
        f"{config['method']} | carriers | PAM k=5\n"
        "2D classical-MDS embedding of the PAM Gower distance",
        loc="left", fontsize=12, fontweight="bold",
    )
    ax.legend(title="PAM assignment", bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False,
              fontsize=8, title_fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    return fig


def main():
    config = CONFIG.copy()
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    clusters = standardise_columns(read_table(config["cluster_file"]), "cluster file")
    annotations = standardise_columns(read_table(config["annotations_file"]), "annotations")
    composition = standardise_columns(read_table(config["features_file"], sheet_name="composition"), "composition") # type: ignore

    cluster_col = config["cluster_column"].lower()
    require_columns(clusters, ["protein_id", "protein_length", cluster_col], "cluster file")
    clusters["protein_length"] = pd.to_numeric(clusters["protein_length"], errors="raise")
    clusters = clusters[["protein_id", "protein_length", cluster_col]].rename(columns={cluster_col: "cluster_id"})
    if clusters["protein_id"].duplicated().any() or clusters["cluster_id"].isna().any():
        raise ValueError("Cluster input must have one non-missing PAM assignment per protein.")

    features = build_carrier_features(clusters, annotations, composition, config["method"])
    binary, categorical, continuous = feature_blocks(features)
    print(f"\nFeature blocks: binary={len(binary)}, categorical={len(categorical)}, continuous={len(continuous)}")
    print(f"Continuous features: {continuous}")

    distance = gower_mixed_distance(features, binary, categorical, continuous)
    coordinates, eigenvalues, explained = classical_mds(distance)
    medoid_indices = pam_medoids(distance, features["cluster_id"])

    result = features[["protein_id", "protein_length", "cluster_id"]].copy()
    result["mds_1"] = coordinates[:, 0]
    result["mds_2"] = coordinates[:, 1]
    result["is_medoid"] = False
    result.loc[result.index[medoid_indices], "is_medoid"] = True

    prefix = f"{config['method']}_carriers_{config['cluster_column']}_gower_mds"
    result.to_csv(output_dir / f"{prefix}_coordinates.csv", index=False)
    pd.DataFrame({"eigenvalue": eigenvalues}).to_csv(output_dir / f"{prefix}_eigenvalues.csv", index=False)

    fig = plot_embedding(result, medoid_indices, explained, config)
    fig.savefig(output_dir / f"{prefix}.png", dpi=config["dpi"], bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote: {output_dir / f'{prefix}.png'}")
    print(f"Wrote: {output_dir / f'{prefix}_coordinates.csv'}")
    print(f"Wrote: {output_dir / f'{prefix}_eigenvalues.csv'}")


if __name__ == "__main__":
    main()
