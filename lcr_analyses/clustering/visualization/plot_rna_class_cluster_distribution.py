#!/usr/bin/env python3
"""Plot RNA-class-normalized cluster distributions: P(cluster | RNA class).

One vertical 100% stacked bar is drawn per RNA-target class. Each class is
normalized independently, so large mRNA/unknown classes do not dominate the
visualization. A CSV summary records class totals, represented clusters, and
the dominant cluster for each RNA class.
"""

from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import pandas as pd


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
METHOD = "AlcoR"  # source_method
CONFIG = {
    "analysis": "global",  # "global" uses PAM k=4; "carriers" uses PAM k=5
    "method": METHOD,  # source_method
    "cluster_file": "lcr_analyses/clustering/global/" + METHOD + "/protein_clusters.csv",
    "output_dir": "lcr_analyses/clustering/visualization/cluster_distribution/global/" + METHOD,
    "cluster_column": None,  # None -> cluster_pam_k4 (global), cluster_pam_k5 (carriers)
    "class_column": "rnaprimaryclass",
    "representation_threshold": 0.05,  # P(cluster | class) >= 5%
    "dpi": 220,
}

# Keep this biological order across all method-specific plots.
RNA_CLASS_ORDER = [
    "mRNA",
    "tRNA",
    "pre-rRNA",
    "snRNA",
    "snoRNA",
    "ncRNA",
    "ribosomal protein",
    "diverse",
    "unknown",
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
# Helpers
# -----------------------------------------------------------------------------
def normalise_class(value):
    """Map superficial formatting variants to the fixed RNA-class labels."""
    if pd.isna(value):
        return "unknown"

    text = str(value).strip().lower()
    text = re.sub(r"[ _]+", " ", text)

    aliases = {
        "mrna": "mRNA",
        "trna": "tRNA",
        "pre rrna": "pre-rRNA",
        "pre-rrna": "pre-rRNA",
        "prerrna": "pre-rRNA",
        "snrna": "snRNA",
        "snorna": "snoRNA",
        "ncrna": "ncRNA",
        "ribosomal protein": "ribosomal protein",
        "ribosomal proteins": "ribosomal protein",
        "diverse": "diverse",
        "unknown": "unknown",
        "unclassified": "unknown",
        "nan": "unknown",
    }
    return aliases.get(text, str(value).strip())


def cluster_sort_key(value):
    """Sort numeric cluster labels naturally, then fall back to text."""
    try:
        return (0, int(value))
    except (ValueError, TypeError):
        return (1, str(value))


def selected_cluster_column(analysis, configured):
    if configured:
        return configured
    if analysis == "global":
        return "cluster_pam_k4"
    if analysis == "carriers":
        return "cluster_pam_k5"
    raise ValueError("analysis must be 'global' or 'carriers'.")


def class_display_label(rna_class, n):
    return f"{rna_class}\n(n={n:,})"


# -----------------------------------------------------------------------------
# Analysis
# -----------------------------------------------------------------------------
def build_tables(df, class_col, cluster_col, threshold):
    df = df[[class_col, cluster_col]].dropna(subset=[cluster_col]).copy()
    df["rna_class"] = df[class_col].map(normalise_class)
    df["cluster_id"] = df[cluster_col]

    observed_classes = set(df["rna_class"])
    class_order = [x for x in RNA_CLASS_ORDER if x in observed_classes]
    class_order += sorted(observed_classes - set(class_order))

    clusters = sorted(df["cluster_id"].unique(), key=cluster_sort_key)

    counts = pd.crosstab(df["rna_class"], df["cluster_id"])
    counts = counts.reindex(index=class_order, columns=clusters, fill_value=0)

    class_totals = counts.sum(axis=1)
    fractions = counts.div(class_totals, axis=0)

    summary_rows = []
    for rna_class in class_order:
        n_class = int(class_totals.loc[rna_class])
        row_counts = counts.loc[rna_class]
        row_fractions = fractions.loc[rna_class]
        dominant_cluster = row_counts.idxmax()
        dominant_n = int(row_counts.max()) # type: ignore
        dominant_fraction = float(row_fractions.loc[dominant_cluster]) # type: ignore
        represented = row_fractions[row_fractions >= threshold]

        summary_rows.append({
            "rna_class": rna_class,
            "n_class": n_class,
            "n_clusters_total": len(clusters),
            "representation_threshold": threshold,
            "n_clusters_represented": int(len(represented)),
            "clusters_represented": "; ".join(
                f"{cluster} ({100 * fraction:.1f}%)"
                for cluster, fraction in represented.items()
            ),
            "dominant_cluster": dominant_cluster,
            "dominant_cluster_n": dominant_n,
            "dominant_cluster_fraction": dominant_fraction,
        })

    summary = pd.DataFrame(summary_rows)
    return counts, fractions, summary


def plot_distribution(fractions, summary, config, cluster_col):
    classes = fractions.index.tolist()
    clusters = fractions.columns.tolist()
    class_n = summary.set_index("rna_class")["n_class"]

    fig_width = max(8.0, 0.95 * len(classes) + 2.5)
    fig, ax = plt.subplots(figsize=(fig_width, 7.0), constrained_layout=True)

    bottoms = pd.Series(0.0, index=classes)
    for i, cluster in enumerate(clusters):
        values = fractions[cluster]
        colour = CLUSTER_COLOURS[i % len(CLUSTER_COLOURS)]

        bars = ax.bar(
            classes,
            values,
            bottom=bottoms,
            color=colour,
            edgecolor="white",
            linewidth=0.6,
            width=0.76,
            label=f"Cluster {cluster}",
        )

        for bar, value, bottom in zip(bars, values, bottoms):
            if value >= 0.08:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bottom + value / 2,
                    f"{100 * value:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="white" if value >= 0.18 else "#222222",
                )

        bottoms = bottoms + values

    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_ylabel("Proteins assigned to cluster within RNA class", fontsize=10)
    ax.set_xlabel("RNA-target superclass", fontsize=10)
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(
        [class_display_label(rna_class, int(class_n.loc[rna_class])) for rna_class in classes],
        fontsize=8,
    )
    ax.set_title(
        f"{config['method']} | {config['analysis']} | {cluster_col}\n"
        r"RNA-class distribution across clusters: $P(\mathrm{cluster}\mid\mathrm{RNA\ class})$",
        fontsize=12,
        fontweight="bold",
        loc="left",
    )
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(
        title="PAM cluster",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=False,
        fontsize=8,
        title_fontsize=8,
    )
    return fig


def main():
    config = CONFIG.copy()
    analysis = config["analysis"].strip().lower()
    cluster_col = selected_cluster_column(analysis, config["cluster_column"])
    class_col = config["class_column"]

    path = Path(config["cluster_file"])
    if not path.exists():
        raise FileNotFoundError(f"Cluster file not found: {path}")

    df = pd.read_csv(path)
    df.columns = [str(col).strip().lower() for col in df.columns]
    class_col = class_col.lower()
    cluster_col = cluster_col.lower()

    required = [class_col, cluster_col]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(
            f"Missing columns: {missing}\n"
            f"Available columns: {df.columns.tolist()}"
        )

    counts, fractions, summary = build_tables(
        df=df,
        class_col=class_col,
        cluster_col=cluster_col,
        threshold=float(config["representation_threshold"]),
    )

    outdir = Path(config["output_dir"])
    outdir.mkdir(parents=True, exist_ok=True)
    stem = f"{config['method']}_{analysis}_{cluster_col}_rna_class_distribution"

    fig = plot_distribution(fractions, summary, config, cluster_col)
    fig.savefig(outdir / f"{stem}.png", dpi=config["dpi"], bbox_inches="tight")
    plt.close(fig)

    counts.to_csv(outdir / f"{stem}_counts.csv")
    fractions.to_csv(outdir / f"{stem}_fractions.csv")
    summary.to_csv(outdir / f"{stem}_summary.csv", index=False)

    print(f"Wrote: {outdir / f'{stem}.png'}")
    print(f"Wrote: {outdir / f'{stem}_counts.csv'}")
    print(f"Wrote: {outdir / f'{stem}_fractions.csv'}")
    print(f"Wrote: {outdir / f'{stem}_summary.csv'}")
    print("\nRNA-class summary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
