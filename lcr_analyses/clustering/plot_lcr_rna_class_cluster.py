#!/usr/bin/env python3
"""Clean PNG visualizations for LCR-cluster × RNA-class overlap results.

Reads:
  <cluster-root>/<analysis>/<method>/protein_clusters.csv
  <cluster-root>/<analysis>/rna_class_overlap/<method>/overlap_pam.xlsx
  <cluster-root>/<analysis>/rna_class_overlap/<method>/overlap_hierarchical.xlsx

Writes PNG figures only. Heatmaps contain no point, star, or text overlays:
cell colour is log2(Haldane-corrected odds ratio); grey means low-information.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.stats.multitest import multipletests

CLASS_ORDER = [
    "mRNA", "pre-rRNA", "ribosomal protein", "tRNA", "snRNA", "snoRNA",
    "ncRNA", "diverse", "unknown",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cluster-root", type=Path, default=Path("lcr_analyses/clustering"))
    p.add_argument("--outdir", type=Path, default=Path("lcr_analyses/clustering/visualization_clean"))
    p.add_argument("--analysis", choices=["global", "carriers", "both"], default="both")
    p.add_argument("--methods", nargs="*", help="Optional exact method names")
    p.add_argument("--global-k", nargs="+", type=int, default=[4])
    p.add_argument("--carrier-k", nargs="+", type=int, default=[4, 5, 6])
    p.add_argument("--min-observed", type=int, default=5)
    p.add_argument("--vmax", type=float, default=3.0)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--significant-only", action="store_true", help="Grey-out q >= alpha cells")
    p.add_argument("--include-carrier-hierarchical", action="store_true")
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def find_methods(analysis_dir: Path, requested: list[str] | None) -> list[str]:
    methods = sorted(p.parent.name for p in analysis_dir.rglob("protein_clusters.csv"))
    if not methods:
        raise FileNotFoundError(f"No protein_clusters.csv files under {analysis_dir}")
    if requested is None:
        return methods
    missing = sorted(set(requested) - set(methods))
    if missing:
        raise FileNotFoundError(f"Requested methods not found: {', '.join(missing)}")
    return requested


def read_enrichment_table(workbook: Path, sheet: str) -> pd.DataFrame:
    raw = pd.read_excel(workbook, sheet_name=sheet, header=None)
    title_rows = raw.index[raw.iloc[:, 0].astype(str).str.strip().eq("enrichment")].tolist()
    if not title_rows:
        raise ValueError(f"No enrichment table in {workbook} [{sheet}]")
    header_row = title_rows[0] + 1
    columns = raw.iloc[header_row].astype(str).tolist()
    table = raw.iloc[header_row + 1:].copy()
    table.columns = columns
    return table.dropna(how="all")


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def load_algorithm(
    workbook: Path, method: str, analysis: str, algorithm: str, requested_k: list[int]
) -> pd.DataFrame:
    if not workbook.exists():
        return pd.DataFrame()
    sheets = pd.ExcelFile(workbook).sheet_names
    wanted = [f"k{k}" for k in requested_k]
    frames = []
    for sheet in wanted:
        if sheet not in sheets:
            continue
        tab = read_enrichment_table(workbook, sheet)
        tab["analysis"] = analysis
        tab["method"] = method
        tab["algorithm"] = algorithm
        tab["run"] = sheet
        frames.append(tab)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    for col in ["cluster", "observed", "expected", "cluster_size", "log2_or", "p_value", "p_fdr"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out["descriptive_only"] = as_bool(out["descriptive_only"])
    out["rna_class"] = out["rna_class"].astype(str)
    return out


def add_family_fdr(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["q_family"] = np.nan
    for _, indices in out.groupby(["analysis", "algorithm"], dropna=False).groups.items():
        pvals = out.loc[indices, "p_value"]
        valid = pvals[pvals.notna()].index
        if len(valid):
            out.loc[valid, "q_family"] = multipletests(out.loc[valid, "p_value"], method="fdr_bh")[1]
    return out


def ordered_classes(df: pd.DataFrame) -> list[str]:
    present = set(df["rna_class"].dropna())
    return [x for x in CLASS_ORDER if x in present] + sorted(present - set(CLASS_ORDER))


def row_labels(panel: pd.DataFrame, clusters: list[int]) -> list[str]:
    sizes = panel.groupby("cluster")["cluster_size"].first().to_dict()
    return [f"C{cluster}  (n={int(sizes[cluster])})" for cluster in clusters]


def draw_heatmap(
    ax: plt.Axes, panel: pd.DataFrame, classes: list[str], min_observed: int,
    significant_only: bool, alpha: float, vmax: float,
) -> None:
    clusters = sorted(panel["cluster"].dropna().astype(int).unique())
    effect = panel.pivot(index="cluster", columns="rna_class", values="log2_or").reindex(clusters, columns=classes)
    observed = panel.pivot(index="cluster", columns="rna_class", values="observed").reindex(clusters, columns=classes)
    low_info = panel.pivot(index="cluster", columns="rna_class", values="descriptive_only").reindex(clusters, columns=classes).fillna(True)
    qvals = panel.pivot(index="cluster", columns="rna_class", values="q_family").reindex(clusters, columns=classes)
    mask = low_info | observed.lt(min_observed) | effect.isna()
    if significant_only:
        mask = mask | qvals.ge(alpha) | qvals.isna()
    sns.heatmap(
        effect, ax=ax, mask=mask, cmap="vlag", center=0, vmin=-vmax, vmax=vmax,
        cbar=False, square=False, linewidths=0.35, linecolor="white",
        xticklabels=classes, yticklabels=row_labels(panel, clusters),
    )
    ax.set_facecolor("#E5E5E5")
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    ax.tick_params(axis="y", labelsize=7, length=0)
    ax.set_xlabel("")
    ax.set_ylabel("")


def save_atlas(
    df: pd.DataFrame, analysis: str, algorithm: str, k_values: list[int], outdir: Path,
    min_observed: int, significant_only: bool, alpha: float, vmax: float, dpi: int,
) -> None:
    subset = df[(df["analysis"] == analysis) & (df["algorithm"] == algorithm) & (df["run"].isin([f"k{k}" for k in k_values]))]
    if subset.empty:
        return
    methods = sorted(subset["method"].unique())
    classes = ordered_classes(subset)
    fig, axes = plt.subplots(
        len(methods), len(k_values), squeeze=False,
        figsize=(3.25 * len(k_values) + 1.15, max(3.0, 1.7 * len(methods) + 1.2)),
        layout="constrained",
    )
    for r, method in enumerate(methods):
        for c, k in enumerate(k_values):
            ax = axes[r, c]
            panel = subset[(subset["method"] == method) & (subset["run"] == f"k{k}")]
            if panel.empty:
                ax.axis("off")
                continue
            draw_heatmap(ax, panel, classes, min_observed, significant_only, alpha, vmax)
            if r == 0:
                ax.set_title(f"{algorithm.capitalize()}  k={k}", fontsize=10, pad=8)
            if c == 0:
                ax.text(-0.58, 0.5, method, transform=ax.transAxes, rotation=90,
                        va="center", ha="center", fontsize=9, fontweight="bold")
    sm = plt.cm.ScalarMappable(cmap="vlag", norm=plt.Normalize(-vmax, vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), shrink=0.72, pad=0.02)
    cbar.set_label("log2 odds ratio", fontsize=8)
    note = f"Grey: observed n < {min_observed} or low-information contingency table."
    if significant_only:
        note += f" Also grey: BH-FDR q ≥ {alpha}."
    fig.text(0.01, 0.008, note, fontsize=8)
    title = f"{analysis.capitalize()} RNA-class overlap: {algorithm}"
    fig.suptitle(title, fontsize=13)
    outdir.mkdir(parents=True, exist_ok=True)
    suffix = "significant" if significant_only else "all_effects"
    filename = f"{analysis}_{algorithm}_atlas_{'_'.join(map(str, k_values))}_{suffix}.png"
    fig.savefig(outdir / filename, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def save_hdbscan_summary(analysis_dir: Path, methods: list[str], analysis: str, outdir: Path, dpi: int) -> None:
    rows = []
    for method in methods:
        path = analysis_dir / method / "protein_clusters.csv"
        data = pd.read_csv(path)
        if "cluster_hdbscan" not in data.columns:
            continue
        vc = data["cluster_hdbscan"].value_counts()
        for label, count in vc.items():
            rows.append({"method": method, "group": "noise" if int(label) == -1 else f"cluster {int(label)}", "fraction": count / len(data)})
    if not rows:
        return
    tab = pd.DataFrame(rows)
    group_order = ["noise"] + sorted([x for x in tab.group.unique() if x != "noise"], key=lambda x: int(x.split()[1]))
    wide = tab.pivot(index="method", columns="group", values="fraction").fillna(0).reindex(columns=group_order)
    fig, ax = plt.subplots(figsize=(8, max(3, 0.55 * len(wide) + 1.4)), layout="constrained")
    palette = ["#9A9A9A"] + list(sns.color_palette("Set2", n_colors=max(1, len(group_order) - 1)))
    left = np.zeros(len(wide))
    for group, color in zip(group_order, palette):
        ax.barh(wide.index, wide[group], left=left, color=color, edgecolor="white", linewidth=0.5, label=group)
        left += wide[group].to_numpy()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of proteins")
    ax.set_title(f"{analysis.capitalize()} HDBSCAN assignment diagnostic")
    ax.invert_yaxis()
    ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{analysis}_hdbscan_assignment_diagnostic.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def save_pam_transitions(analysis_dir: Path, methods: list[str], k_values: list[int], outdir: Path, dpi: int) -> None:
    pairs = list(zip(k_values[:-1], k_values[1:]))
    if not pairs:
        return
    fig, axes = plt.subplots(len(methods), len(pairs), squeeze=False, figsize=(4.1 * len(pairs), max(3, 2.35 * len(methods))), layout="constrained")
    for r, method in enumerate(methods):
        data = pd.read_csv(analysis_dir / method / "protein_clusters.csv")
        for c, (ka, kb) in enumerate(pairs):
            ax = axes[r, c]
            a, b = f"cluster_pam_k{ka}", f"cluster_pam_k{kb}"
            if a not in data.columns or b not in data.columns:
                ax.axis("off")
                continue
            tab = pd.crosstab(data[a], data[b])
            sns.heatmap(tab, cmap="Blues", annot=True, fmt="d", cbar=False, linewidths=0.35, linecolor="white", ax=ax)
            ax.set_xlabel(f"PAM k={kb}")
            ax.set_ylabel(f"PAM k={ka}")
            if r == 0:
                ax.set_title(f"k={ka} → k={kb}")
            if c == 0:
                ax.text(-0.62, 0.5, method, transform=ax.transAxes, rotation=90,
                        va="center", ha="center", fontsize=9, fontweight="bold")
    fig.suptitle("Carrier PAM membership transitions", fontsize=13)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / "carriers_pam_k_transitions.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    sns.set_theme(style="white", context="paper")
    analyses = ["global", "carriers"] if args.analysis == "both" else [args.analysis]
    tables = []
    analysis_info = []
    for analysis in analyses:
        analysis_dir = args.cluster_root / analysis
        methods = find_methods(analysis_dir, args.methods)
        k_values = args.global_k if analysis == "global" else args.carrier_k
        analysis_info.append((analysis, analysis_dir, methods, k_values))
        for method in methods:
            overlap_dir = analysis_dir / "rna_class_overlap" / method
            for algorithm in ("pam", "hierarchical"):
                workbook = overlap_dir / f"overlap_{algorithm}.xlsx"
                table = load_algorithm(workbook, method, analysis, algorithm, k_values)
                if not table.empty:
                    tables.append(table)
    if not tables:
        raise FileNotFoundError("No overlap PAM/hierarchical workbooks found in the expected directories.")
    enrichment = add_family_fdr(pd.concat(tables, ignore_index=True))
    args.outdir.mkdir(parents=True, exist_ok=True)
    enrichment.to_csv(args.outdir / "visualization_enrichment_table.csv", index=False)
    for analysis, analysis_dir, methods, k_values in analysis_info:
        target = args.outdir / analysis
        save_atlas(enrichment, analysis, "pam", k_values, target, args.min_observed,
                   args.significant_only, args.alpha, args.vmax, args.dpi)
        if analysis == "global" or args.include_carrier_hierarchical:
            save_atlas(enrichment, analysis, "hierarchical", k_values, target, args.min_observed,
                       args.significant_only, args.alpha, args.vmax, args.dpi)
        save_hdbscan_summary(analysis_dir, methods, analysis, target, args.dpi)
        if analysis == "carriers":
            save_pam_transitions(analysis_dir, methods, k_values, target, args.dpi)
    print(f"Done. PNG figures written under: {args.outdir.resolve()}")


if __name__ == "__main__":
    main()
