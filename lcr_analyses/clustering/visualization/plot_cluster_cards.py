#!/usr/bin/env python3
"""Create minimal PAM cluster cards for LCR clustering analyses.

Edit CONFIG paths and run once for Analysis A (global, k=4) and once for
Analysis B (carriers, k=5). The script reconstructs the feature matrix,
recomputes exact PAM medoids under the supplied equal-block Gower function,
and writes one PNG card and one auditable CSV per cluster.
"""

from pathlib import Path
import re
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from sklearn.metrics import pairwise_distances


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
CONFIG = {
    "analysis": "global",  # "global" -> PAM k=4; "carriers" -> PAM k=5
    "method": "SEG_intermediate",         # source_method
    "cluster_file": "lcr_analyses/clustering/global/SEG_intermediate/protein_clusters.csv",
    "annotations_file": "lcr_analyses/pc_properties/lcr_annotations.xlsx",
    "features_file": "lcr_analyses/pc_properties/lcr_features.xlsx",
    "pfam_file": "rbp_lcrs/pfam_lcr_overlap.xlsx",
    "pfam_sheet": "Pfam_hits",
    "output_dir": "lcr_analyses/clustering/visualization/cluster_cards/global/SEG_intermediate",
    "cluster_id_column": "cluster_pam_k4", 
    "coverage_mode": "sum_lcr_lengths",
    "dpi": 220,
    "max_domain_labels": 8,
}

COMPOSITION_METHOD_MAP = {
    "flps": "flps_strict",
    "seg": "seg_strict",
}

POSITION_ORDER = [
    "domain_intrinsic",
    "domain_edge",
    "domain_adjacent",
    "interdomain_linker",
    "distal_terminal",
]

DISPLAY_POSITION = {
    "domain_intrinsic": "domain-intrinsic",
    "domain_edge": "domain-edge",
    "domain_adjacent": "domain-adjacent",
    "interdomain_linker": "interdomain-linker",
    "distal_terminal": "distal-terminal",
}

SIGNATURE_COLOURS = {
    "rg_rgg_repeat_region": "#8e63ce",
    "sr_rs_repeat_region": "#e06c9f",
    "gs_rich_neutral_region": "#44a6a1",
    "basic_enriched_region": "#d88927",
    "arginine_rich_rna_contact_region": "#c85454",
    "acidic_region": "#4a90c2",
    "mixed_charge_polyampholyte": "#897862",
    "polar_linker": "#65a86e",
    "hydrophobic_region": "#8a6d3b",
    "aromatic_sticker_region": "#a879cc",
    "aromatic_aggregation_prone_region": "#87595b",
    "cation_pi_rich_neighborhood": "#b965a4",
    "disorder_rich_spacer_region": "#74a9cf",
    "serine_rich_region": "#f08b61",
    "glutamine_rich_region": "#d5aa3d",
    "glycine_rich_region": "#5b9bd5",
    "proline_rich_disordered_region": "#a89f4d",
    "unmapped_physicochemical_properties": "#8c8c8c",
}
FALLBACK_COLOURS = list(plt.get_cmap("tab20").colors) # type: ignore


# -----------------------------------------------------------------------------
# Input helpers and validation
# -----------------------------------------------------------------------------
def read_table(path, sheet_name=0):
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Input not found: {path}")

    if path.suffix.lower() in {".xlsx", ".xls"}:
        print(f"\nReading: {path}")
        print(f"Requested sheet: {sheet_name}")

        workbook = pd.ExcelFile(path)
        print(f"Available sheets: {workbook.sheet_names}")

        df = pd.read_excel(path, sheet_name=sheet_name)

        print("Parsed columns:")
        print([repr(col) for col in df.columns.tolist()])
        print("First 3 rows:")
        print(df.head(3).to_string())

        return df

    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path)
        print(f"\nReading CSV: {path}")
        print([repr(col) for col in df.columns.tolist()])
        return df

    raise ValueError(f"Unsupported input format: {path}")


def normalise_text(value):
    if pd.isna(value):
        return np.nan
    value = str(value).strip().lower()
    value = re.sub(r"[ /-]+", "_", value)
    value = re.sub(r"_+", "_", value)
    return value


def standardise_columns(df, source_name):
    df = df.copy()

    def clean_column_name(col):
        col = str(col)
        col = col.replace("\ufeff", "")
        col = col.replace("\xa0", " ")
        col = col.replace("\n", " ")
        col = col.strip().lower()
        col = re.sub(r"\s+", "_", col)
        col = re.sub(r"_+", "_", col)
        return col

    df.columns = [clean_column_name(col) for col in df.columns]

    aliases = {
        "uniprot_accession": "protein_id",
        "uniprot_id": "protein_id",
        "accession": "protein_id",
        "proteinid": "protein_id",
        "protein_id": "protein_id",

        "sourcemethod": "source_method",
        "source_method": "source_method",
        "method": "source_method",

        "uniprot_length": "protein_length",
        "proteinlength": "protein_length",
        "protein_length": "protein_length",
        "length_protein": "protein_length",
    }

    df = df.rename(
        columns={
            col: aliases[col]
            for col in df.columns
            if col in aliases
        }
    )

    if "protein_id" not in df.columns:
        raise ValueError(
            f"{source_name}: no protein ID column found.\n"
            f"Parsed columns: {df.columns.tolist()}"
        )

    # Convert composite annotation IDs to canonical UniProt accessions.
    # Example:
    # Q9Y2T7|ensembl_protein_id=...  -->  Q9Y2T7
    df["protein_id"] = (
        df["protein_id"]
        .astype(str)
        .str.strip()
        .str.split("|", n=1)
        .str[0]
        .str.upper()
    )

    for col in (
        "source_method",
        "domain_position_class",
        "primary_physicochemical_annotation",
    ):
        if col in df.columns:
            df[col] = df[col].map(normalise_text)

    return df


def require_columns(df, columns, source_name):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"{source_name}: missing required columns: {missing}")


def selected_cluster_column(analysis, configured):
    if configured:
        return configured
    return "cluster_pam_k4" if analysis == "global" else "cluster_pam_k5"


# -----------------------------------------------------------------------------
# Feature reconstruction
# -----------------------------------------------------------------------------
def weighted_average(group, metric):
    valid = group[[metric, "length"]].dropna()
    if valid.empty or valid["length"].sum() <= 0:
        return np.nan
    return np.average(valid[metric].astype(float), weights=valid["length"].astype(float))


def top_two_presence(cluster_features, columns, prefix, labels):
    records = []
    n = len(cluster_features)
    for col in columns:
        count = int(cluster_features[col].sum())
        if count:
            key = col.removeprefix(prefix)
            records.append((labels.get(key, key.replace("_", " ")), count, count / n, key))
    records.sort(key=lambda x: (-x[1], x[0]))
    while len(records) < 2:
        records.append(("none detected", 0, 0.0, "none"))
    return records[:2]


def build_protein_features(clusters, annotations, composition, analysis, method):
    method_key = normalise_text(method)
    composition_method_key = COMPOSITION_METHOD_MAP.get(method_key, method_key) # type: ignore
    annotations = annotations.loc[annotations["source_method"] == method_key].copy()
    composition = composition.loc[
    composition["source_method"] == composition_method_key].copy()
    composition["source_method"] = method_key

    protein_ids = pd.Index(clusters["protein_id"].drop_duplicates(), name="protein_id")
    annotations = annotations.loc[annotations["protein_id"].isin(protein_ids)].copy()
    composition = composition.loc[composition["protein_id"].isin(protein_ids)].copy()

    print("\n--- Protein-ID matching check ---")
    print(f"Cluster proteins: {len(protein_ids):,}")
    print(f"Annotation LCR rows after method + ID filtering: {len(annotations):,}")
    print(f"Annotation proteins after method + ID filtering: {annotations['protein_id'].nunique():,}")
    print(f"Composition rows after method + ID filtering: {len(composition):,}")
    print(f"Composition proteins after method + ID filtering: {composition['protein_id'].nunique():,}")

    if annotations.empty:
        cluster_example = list(protein_ids[:5])
        annotation_example = annotations["protein_id"].head(5).tolist()

        raise ValueError(
            "No annotation rows remain after matching annotation protein IDs "
            "to cluster protein IDs. This indicates an ID-format mismatch.\n"
            f"Example cluster IDs: {cluster_example}\n"
            f"Annotation IDs after filtering: {annotation_example}"
        )

    key = ["protein_id", "source_method", "start", "end"]
    require_columns(annotations, key + ["domain_position_class", "primary_physicochemical_annotation"], "annotations")
    require_columns(composition, key + ["ncpr", "fcr"], "composition sheet")
    if annotations.duplicated(key).any():
        raise ValueError("annotations has duplicate LCR coordinate keys.")
    if composition.duplicated(key).any():
        raise ValueError("composition sheet has duplicate LCR coordinate keys.")

    lcr = annotations.merge(
        composition[key + ["ncpr", "fcr"]],
        on=key,
        how="left",
        validate="one_to_one",
        indicator="_composition_merge",
    )

    missing_composition = lcr.loc[
        lcr["_composition_merge"] == "left_only"
    ].copy()

    if not missing_composition.empty:
        print("\n--- LCRs missing composition after join ---")
        print(
            missing_composition[
                [
                    "protein_id",
                    "source_method",
                    "start",
                    "end",
                    "length",
                    "primary_physicochemical_annotation",
                ]
            ].to_string(index=False)
        )

    lcr = lcr.drop(columns="_composition_merge")
    if analysis == "carriers" and lcr[["ncpr", "fcr"]].isna().any().any():
        n_missing = int(lcr[["ncpr", "fcr"]].isna().any(axis=1).sum())
        raise ValueError(f"{n_missing} LCRs lack ncpr/fcr after annotation-composition join.")

    base = clusters[["protein_id", "protein_length"]].drop_duplicates().set_index("protein_id")
    if base["protein_length"].isna().any() or (base["protein_length"] <= 0).any():
        raise ValueError("protein_length must be present and positive for every clustered protein.")

    lcr["length"] = pd.to_numeric(lcr["length"], errors="raise")
    architecture = lcr.groupby("protein_id", observed=True).agg(
        n_lcr=("length", "size"),
        lcr_residues=("length", "sum"),
        mean_lcr_length=("length", "mean"),
    )
    protein = base.join(architecture).fillna({"n_lcr": 0, "lcr_residues": 0})
    protein["n_lcr"] = protein["n_lcr"].astype(int)
    protein["coverage"] = protein["lcr_residues"] / protein["protein_length"]
    protein["multiplicity_cat"] = np.select(
        [protein["n_lcr"].eq(0), protein["n_lcr"].eq(1), protein["n_lcr"].eq(2)],
        ["0", "1", "2"], default="3+"
    )

    observed_positions = [x for x in POSITION_ORDER if x in set(annotations["domain_position_class"].dropna())]
    observed_signatures = sorted(set(annotations["primary_physicochemical_annotation"].dropna()))
    for position in observed_positions:
        ids = annotations.loc[annotations["domain_position_class"] == position, "protein_id"].unique()
        protein[f"pos__{position}"] = protein.index.isin(ids).astype(int)
    for signature in observed_signatures:
        ids = annotations.loc[
            annotations["primary_physicochemical_annotation"] == signature, "protein_id"
        ].unique()
        protein[f"sig__{signature}"] = protein.index.isin(ids).astype(int)

    if analysis == "carriers":
        carriers = protein.index[protein["n_lcr"] > 0]
        protein = protein.loc[carriers].copy()

        for metric in ("ncpr", "fcr"):
            valid = lcr[["protein_id", "length", metric]].dropna().copy()

            numerator = (
                (valid[metric].astype(float) * valid["length"].astype(float))
                .groupby(valid["protein_id"], observed=True)
                .sum()
            )

            denominator = (
                valid["length"]
                .astype(float)
                .groupby(valid["protein_id"], observed=True)
                .sum()
            )

            values = numerator / denominator

            protein[metric] = values.reindex(protein.index).to_numpy()

        if protein[["ncpr", "fcr"]].isna().any().any():
            raise ValueError("Carrier protein has undefined weighted ncpr or fcr.")

    protein = protein.reset_index()
    clusters = clusters.merge(protein, on=["protein_id", "protein_length"], how="inner", validate="one_to_one")
    if len(clusters) != len(protein):
        raise ValueError("Cluster membership and reconstructed protein features do not match.")
    return clusters, lcr


# -----------------------------------------------------------------------------
# Exact supplied Gower implementation and medoids
# -----------------------------------------------------------------------------
def gower_mixed_distance(df, binary_cols, categorical_cols, continuous_cols):
    """Equal-weight three-block mixed distance in [0, 1]."""
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


def infer_feature_blocks(features, analysis):
    binary = [c for c in features.columns if c.startswith(("pos__", "sig__"))]
    categorical = ["multiplicity_cat"] if analysis == "global" else []
    continuous = ["coverage"] if analysis == "global" else ["n_lcr", "coverage", "mean_lcr_length", "ncpr", "fcr"]
    binary = [c for c in binary if features[c].nunique(dropna=False) > 1]
    categorical = [c for c in categorical if features[c].nunique(dropna=False) > 1]
    continuous = [c for c in continuous if features[c].nunique(dropna=False) > 1]
    return binary, categorical, continuous


def calculate_medoids(features, cluster_col, binary, categorical, continuous):
    distance = gower_mixed_distance(features, binary, categorical, continuous)
    medoids = []
    for cluster_id, row_idx in features.groupby(cluster_col, sort=True).groups.items():
        indices = np.asarray(list(row_idx))
        within = distance[np.ix_(indices, indices)]
        totals = within.sum(axis=1)
        best_local_index = np.argmin(totals)
        chosen = indices[best_local_index]

        medoid_total_distance = float(totals[best_local_index])

        if len(indices) > 1:
            medoid_mean_distance = medoid_total_distance / (len(indices) - 1)
        else:
            medoid_mean_distance = 0.0

        medoids.append({
            "cluster_id": cluster_id,
            "medoid_row": chosen,
            "medoid_total_within_distance": medoid_total_distance,
            "medoid_mean_within_distance": medoid_mean_distance,
        })
    medoids = pd.DataFrame(medoids)
    features = features.copy()
    features["is_medoid"] = False
    features.loc[medoids["medoid_row"], "is_medoid"] = True
    return features, medoids


# -----------------------------------------------------------------------------
# Drawing
# -----------------------------------------------------------------------------
def signature_colour(signature, fallback_index):
    return SIGNATURE_COLOURS.get(signature, FALLBACK_COLOURS[fallback_index % len(FALLBACK_COLOURS)])


def format_fraction(x):
    return f"{x:.3f}"


def plot_profile(ax, summary, analysis):
    rows = [
        (
            "Coverage (median)",
            summary["coverage_median"],
            "continuous_01",
            None,
        ),
        (
            summary["top_position_1"],
            summary["top_position_1_fraction"],
            "fraction",
            "top_position_1",
        ),
        (
            summary["top_position_2"],
            summary["top_position_2_fraction"],
            "fraction",
            "top_position_2",
        ),
        (
            summary["top_signature_1"],
            summary["top_signature_1_fraction"],
            "fraction",
            "top_signature_1",
        ),
        (
            summary["top_signature_2"],
            summary["top_signature_2_fraction"],
            "fraction",
            "top_signature_2",
        ),
    ]

    if analysis == "carriers":
        rows.extend([
            (
                "NCPR (median)",
                summary["ncpr_median"],
                "signed",
                None,
            ),
            (
                "FCR (median)",
                summary["fcr_median"],
                "continuous_01",
                None,
            ),
        ])

    ax.set_xlim(-1.05, 1.22)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.axis("off")

    for i, (label, value, kind, summary_prefix) in enumerate(rows[::-1]):
        y = i

        ax.text(
            -1.03,
            y,
            label,
            ha="left",
            va="center",
            fontsize=8,
        )

        if kind in {"continuous_01", "fraction"}:
            shown = max(0.0, min(1.0, float(value)))

            ax.barh(
                y,
                1.0,
                left=0,
                color="#eeeeee",
                height=0.48,
            )

            ax.barh(
                y,
                shown,
                left=0,
                color="#4c78a8",
                height=0.48,
            )

            if kind == "fraction" and summary_prefix is not None:
                count = summary[f"{summary_prefix}_count"]
                text = (
                    f"{100 * value:.1f}% "
                    f"({count}/{summary['cluster_n']})"
                )
            else:
                text = format_fraction(value)

            ax.text(
                1.03,
                y,
                text,
                ha="left",
                va="center",
                fontsize=8,
            )

        else:
            # NCPR: bar length = absolute magnitude; colour = charge direction.
            magnitude = min(abs(float(value)), 1.0)

            if value > 0:
                colour = "#c65b5b"   # positive net charge: muted red
            elif value < 0:
                colour = "#e6c84f"   # negative net charge: muted yellow
            else:
                colour = "#bdbdbd"   # charge-balanced: grey

            ax.barh(
                y,
                1.0,
                left=0,
                color="#eeeeee",
                height=0.48,
            )

            ax.barh(
                y,
                magnitude,
                left=0,
                color=colour,
                height=0.48,
            )

            ax.text(
                1.03,
                y,
                f"{value:+.3f}",
                ha="left",
                va="center",
                fontsize=8,
            )

def plot_architecture(ax, medoid, lcr, domains, max_domain_labels):
    protein_id = medoid["protein_id"]
    length = float(medoid["protein_length"])
    m_lcr = lcr.loc[lcr["protein_id"] == protein_id].sort_values(["start", "end"])
    m_domains = domains.loc[domains["protein_id"] == protein_id].sort_values(["domain_start", "domain_end"])

    ax.set_xlim(0, length)
    ax.set_ylim(-1.05, 1.35)
    ax.set_yticks([])
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.spines["bottom"].set_position(("data", -0.73))
    ax.set_xlabel(f"Residue coordinate (protein length: {int(length)} aa)", fontsize=8)
    ax.tick_params(axis="x", labelsize=7, length=3)
    # Main protein backbone
    ax.hlines(
        0,
        1,
        length,
        color="#666666",
        linewidth=2.0,
        zorder=1,
    )

    # Domain annotation track: connects all Pfam blocks visually
    domain_track_y = 0.60

    ax.hlines(
        domain_track_y,
        1,
        length,
        color="#b0b0b0",
        linewidth=1.0,
        zorder=1,
    )

    ax.text(1, -0.45, "1", ha="left", va="top", fontsize=7)
    ax.text(length, -0.45, str(int(length)), ha="right", va="top", fontsize=7)

    for i, (_, row) in enumerate(m_domains.head(max_domain_labels).iterrows()):
        start, end = float(row["domain_start"]), float(row["domain_end"])
        ax.add_patch(Rectangle((start, 0.45), end - start + 1, 0.30,
                               facecolor="#bdbdbd", edgecolor="#555555", lw=0.5, zorder=2))
        label = str(row["pfam_name"])
        ax.text((start + end) / 2, 0.90 + (i % 2) * 0.22, f"{label}\n{int(start)}–{int(end)}",
                ha="center", va="bottom", fontsize=6, clip_on=True)
    if len(m_domains) > max_domain_labels:
        ax.text(length, 1.25, f"+{len(m_domains) - max_domain_labels} domains", ha="right", fontsize=6)

    # All LCRs are drawn directly on the same protein-backbone line.
    lcr_y = -0.10
    lcr_height = 0.20

    for i, (_, row) in enumerate(m_lcr.iterrows()):
        start = float(row["start"])
        end = float(row["end"])

        signature = row["primary_physicochemical_annotation"]
        colour = signature_colour(signature, i)

        ax.add_patch(Rectangle(
            (start, lcr_y),
            end - start + 1,
            lcr_height,
            facecolor=colour,
            edgecolor="#333333",
            lw=0.45,
            zorder=3,
        ))

        label = (
            f"{signature.replace('_', ' ')}\n"
            f"{int(start)}–{int(end)}"
        )

        ax.text(
            (start + end) / 2,
            lcr_y - 0.08,
            label,
            ha="center",
            va="top",
            fontsize=5.5,
            clip_on=True,
        )


def make_card(summary, medoid, lcr, domains, config):
    fig = plt.figure(figsize=(10.2, 5.8), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[1.35, 1.0])
    profile_ax = fig.add_subplot(grid[0])
    architecture_ax = fig.add_subplot(grid[1])

    title = (
        f"{config['method']} | {config['analysis']} | {summary['cluster_column']} | "
        f"Cluster {summary['cluster_id']} | n = {summary['cluster_n']}"
    )
    subtitle = (
        f"PAM medoid: {medoid['protein_id']}  |  "
        f"within-cluster mean distance: {summary['medoid_mean_within_distance']:.3f}"
    )
    fig.suptitle(title, x=0.02, ha="left", fontsize=12, fontweight="bold")
    profile_ax.set_title(subtitle, loc="left", fontsize=8.5, pad=8)
    plot_profile(profile_ax, summary, config["analysis"])
    plot_architecture(architecture_ax, medoid, lcr, domains, config["max_domain_labels"])
    return fig


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    config = CONFIG.copy()
    analysis = config["analysis"].strip().lower()
    if analysis not in {"global", "carriers"}:
        raise ValueError("CONFIG['analysis'] must be 'global' or 'carriers'.")
    cluster_col = selected_cluster_column(analysis, config["cluster_id_column"])
    outdir = Path(config["output_dir"])
    outdir.mkdir(parents=True, exist_ok=True)

    clusters = standardise_columns(read_table(config["cluster_file"]), "cluster file")
    annotations = standardise_columns(read_table(config["annotations_file"]), "annotations")
    composition = standardise_columns(read_table(config["features_file"], sheet_name="composition"), "composition sheet") # type: ignore
    domains = standardise_columns(read_table(config["pfam_file"], sheet_name=config["pfam_sheet"]), "Pfam table")

    require_columns(clusters, ["protein_id", "protein_length", cluster_col], "cluster file")
    require_columns(domains, ["protein_id", "pfam_name", "domain_start", "domain_end"], "Pfam table")
    clusters["protein_length"] = pd.to_numeric(clusters["protein_length"], errors="raise")
    clusters = clusters[["protein_id", "protein_length", cluster_col]].drop_duplicates()
    clusters = clusters.rename(columns={cluster_col: "cluster_id"})
    if clusters["cluster_id"].isna().any():
        raise ValueError(f"{cluster_col} contains missing cluster assignments.")
    if clusters["protein_id"].duplicated().any():
        raise ValueError("The selected cluster file has duplicate protein IDs.")

    features, lcr = build_protein_features(clusters, annotations, composition, analysis, config["method"])
    binary, categorical, continuous = infer_feature_blocks(features, analysis)
    print("\n--- Reconstructed feature matrix ---")
    print(f"Proteins entering medoid calculation: {len(features):,}")
    print(f"Binary features ({len(binary)}): {binary}")
    print(f"Categorical features ({len(categorical)}): {categorical}")
    print(f"Continuous features ({len(continuous)}): {continuous}")

    if features.empty:
        raise ValueError(
            "No proteins remain after feature reconstruction. "
            "For a carrier analysis, this almost certainly means no LCR "
            "annotations matched the protein_clusters accessions."
        )
    features, medoids = calculate_medoids(features, "cluster_id", binary, categorical, continuous)

    position_columns = [c for c in features.columns if c.startswith("pos__")]
    signature_columns = [c for c in features.columns if c.startswith("sig__")]
    summaries = []
    for cluster_id, group in features.groupby("cluster_id", sort=True):
        medoid = group.loc[group["is_medoid"]].iloc[0]
        top_pos = top_two_presence(group, position_columns, "pos__", DISPLAY_POSITION)
        top_sig = top_two_presence(group, signature_columns, "sig__", {})
        distance_row = medoids.loc[medoids["cluster_id"] == cluster_id].iloc[0]
        summary = {
            "method": config["method"],
            "analysis": analysis,
            "cluster_column": cluster_col,
            "cluster_id": cluster_id,
            "cluster_n": len(group),
            "medoid_protein_id": medoid["protein_id"],
            "medoid_protein_length": int(medoid["protein_length"]),
            "medoid_total_within_distance": distance_row["medoid_total_within_distance"],
            "medoid_mean_within_distance": distance_row["medoid_mean_within_distance"],
            "coverage_median": group["coverage"].median(),
            "ncpr_median": group["ncpr"].median() if analysis == "carriers" else np.nan,
            "fcr_median": group["fcr"].median() if analysis == "carriers" else np.nan,
        }
        for rank, (label, count, fraction, _) in enumerate(top_pos, 1):
            summary[f"top_position_{rank}"] = label
            summary[f"top_position_{rank}_count"] = count
            summary[f"top_position_{rank}_fraction"] = fraction
        for rank, (label, count, fraction, _) in enumerate(top_sig, 1):
            summary[f"top_signature_{rank}"] = label
            summary[f"top_signature_{rank}_count"] = count
            summary[f"top_signature_{rank}_fraction"] = fraction
        summaries.append(summary)

        fig = make_card(summary, medoid, lcr, domains, config)
        output = outdir / f"cluster_card_{config['method']}_{analysis}_{cluster_col}_cluster_{cluster_id}.png"
        fig.savefig(output, dpi=config["dpi"], bbox_inches="tight")
        plt.close(fig)

    summary_df = pd.DataFrame(summaries).sort_values("cluster_id")
    summary_df.to_csv(outdir / f"cluster_card_summary_{config['method']}_{analysis}_{cluster_col}.csv", index=False)
    features.to_csv(outdir / f"reconstructed_features_{config['method']}_{analysis}_{cluster_col}.csv", index=False)
    print(f"Wrote {len(summary_df)} cards and summary tables to: {outdir}")
    print("Feature blocks used for medoid recalculation:")
    print(f"  binary ({len(binary)}): {binary}")
    print(f"  categorical ({len(categorical)}): {categorical}")
    print(f"  continuous ({len(continuous)}): {continuous}")


if __name__ == "__main__":
    main()
