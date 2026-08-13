#!/usr/bin/env python3
"""Phase 2.4.1 consolidation layer.

Walks the per-method outputs produced by phase_2_4_1_lcr_analysis_by_method.py
and builds a single cross-method summary under <outdir>/summary/:

  summary/tables/
    master_lcr_continuous.csv          - all methods' LCR-level continuous summaries stacked
    master_protein_continuous.csv      - all methods' protein-level continuous summaries stacked
    master_lcr_categorical.csv         - all methods' LCR-level categorical enrichments stacked
    master_protein_categorical.csv     - all methods' protein-level categorical enrichments stacked
    master_protein_core_metrics.csv    - per-protein n_lcr / coverage, all methods
    master_qc_overview.csv             - QC overview, all methods
    headline_numbers.csv               - one row per method: key figures at a glance
    median_coverage_by_class_method.csv     - RNA class x method (median coverage)
    median_length_by_class_method.csv       - RNA class x method (median LCR length)
    signature_enrichment_key.csv       - log2(obs/exp) for key signatures, class x method
    top_signature_per_class.csv        - strongest enriched signature per RNA class per method

  summary/figures/
    coverage_heatmap.png               - median coverage, RNA class x method
    length_heatmap.png                 - median LCR length, RNA class x method
    signature_enrichment_heatmap.png   - log2 enrichment of key signatures, class x method facets
    n_lcr_by_method.png                - mean LCRs per protein per method
    coverage_by_method.png             - median coverage per method

  RESULTS_GUIDE.md                     - which output answers which analysis question

Run after the per-method analysis:
    python phase_2_4_1_consolidate.py --outdir lcr_analyses/results
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_theme(style="whitegrid", context="notebook")

# Signatures highlighted in the cross-method enrichment heatmap (edit freely).
KEY_SIGNATURES = [
    "basicenrichedregion", "argininerichrnacontactregion", "acidicregion",
    "mixedchargepolyampholyte", "disorderrichspacerregion",
    "rgrggrepeatregion", "srrsrepeatregion", "gsrichneutralregion",
    "aromaticstickerregion", "prolinerichdisorderedregion",
]


def read_method_table(method_dir: Path, name: str) -> pd.DataFrame | None:
    path = method_dir / "tables" / name
    if not path.exists():
        return None
    return pd.read_csv(path)


def collect(results_root: Path) -> dict[str, pd.DataFrame]:
    methods_root = results_root
    stores: dict[str, list[pd.DataFrame]] = {
        "lcr_continuous": [], "protein_continuous": [],
        "lcr_categorical": [], "protein_categorical": [],
        "protein_core": [], "qc": [],
    }
    for method_dir in sorted(p for p in methods_root.iterdir() if p.is_dir()):
        mapping = {
            "lcr_continuous": "lcr_level_continuous_summary.csv",
            "protein_continuous": "protein_level_continuous_summary.csv",
            "lcr_categorical": "lcr_level_categorical_enrichment.csv",
            "protein_categorical": "protein_level_categorical_enrichment.csv",
            "protein_core": "protein_level_core_metrics.csv",
            "qc": "qc_overview.csv",
        }
        for key, fname in mapping.items():
            df = read_method_table(method_dir, fname)
            if df is not None:
                if "method" not in df.columns:
                    df["method"] = method_dir.name
                stores[key].append(df)
    return {k: pd.concat(v, ignore_index=True) for k, v in stores.items() if v}


def pivot_metric(df: pd.DataFrame, metric: str, value: str) -> pd.DataFrame:
    sub = df.loc[df["metric"] == metric]
    return sub.pivot_table(index="rna_class", columns="method", values=value, aggfunc="first")


def plot_heatmap(mat: pd.DataFrame, title: str, outfile: Path, cmap: str = "viridis", center: float | None = None, fmt: str = ".2f") -> None:
    if mat.empty:
        return
    fig, ax = plt.subplots(figsize=(max(6, 1.2 * mat.shape[1] + 2), max(4, .5 * mat.shape[0] + 1.5)))
    kwargs = dict(cmap=cmap, annot=True, fmt=fmt, linewidths=.4, linecolor="white", ax=ax)
    if center is not None:
        vmax = max(.5, float(np.nanmax(np.abs(mat.to_numpy()))))
        kwargs.update(center=center, vmin=-vmax, vmax=vmax)
    sns.heatmap(mat, **kwargs)
    ax.set(title=title, xlabel="Method", ylabel="RNA-target superclass")
    fig.tight_layout(); fig.savefig(outfile, dpi=220); plt.close(fig)


def build_headline(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    qc = data.get("qc")
    core = data.get("protein_core")
    cont = data.get("lcr_continuous")
    cat = data.get("lcr_categorical")
    methods = qc["method"].tolist() if qc is not None else sorted(core["method"].unique())
    for method in methods:
        row: dict[str, object] = {"method": method}
        if qc is not None:
            q = qc.loc[qc["method"] == method].iloc[0]
            row.update(n_lcr=int(q["n_lcr"]), n_proteins=int(q["n_proteins"]),
                       length_match=f"{q['fraction_matched_uniprot_length']:.1%}")
        if core is not None:
            c = core.loc[core["method"] == method]
            row.update(median_coverage=c["coverage"].median(),
                       pct_multi_lcr=f"{(c['n_lcr'] >= 2).mean():.1%}",
                       mean_n_lcr=c["n_lcr"].mean())
        if cont is not None:
            l = cont.loc[(cont["method"] == method) & (cont["metric"] == "length")]
            row["median_lcr_length"] = l["median"].median() if not l.empty else np.nan
        if cat is not None:
            s = cat.loc[(cat["method"] == method) & (cat["analysis"] == "signature_sequence_count_lcr")]
            if not s.empty:
                top = s.loc[s.groupby("rna_class")["log2_obs_exp"].idxmax()]
                row["n_class_signature_peaks"] = len(top)
        rows.append(row)
    return pd.DataFrame(rows)


def build_top_signatures(cat: pd.DataFrame) -> pd.DataFrame:
    s = cat.loc[cat["analysis"] == "signature_sequence_count_lcr"].copy()
    if s.empty:
        return s
    idx = s.groupby(["method", "rna_class"])["log2_obs_exp"].idxmax()
    top = s.loc[idx, ["method", "rna_class", "signature", "log2_obs_exp", "p_fdr"]]
    return top.sort_values(["method", "rna_class"]).reset_index(drop=True)


def write_results_guide(outdir: Path) -> None:
    (outdir / "RESULTS_GUIDE.md").write_text(
        "# Phase 2.4.1 — Results Guide\n\n"
        "## Reading order\n"
        "1. `summary/tables/headline_numbers.csv` — per-method key figures; check methods detect comparable amounts before interpreting biology.\n"
        "2. `summary/figures/coverage_heatmap.png` and `length_heatmap.png` — do methods grossly agree on scale?\n"
        "3. `summary/figures/signature_enrichment_heatmap.png` — cross-method view of physicochemical enrichment by RNA class.\n"
        "4. Per-method folders (`methods/<method>/`) — drill-down archive; open only when a summary signal needs detail.\n\n"
        "## Which output answers which question\n"
        "| Analysis question | Summary layer | Per-method source |\n"
        "|---|---|---|\n"
        "| 1. Domain position x RNA class | `master_lcr_categorical.csv` (analysis=`position_lcr`), `master_protein_categorical.csv` (`position_presence`) | `lcr_position_stacked.png` |\n"
        "| 2. Physicochemical signature x RNA class | `signature_enrichment_key.csv`, `top_signature_per_class.csv`, `signature_enrichment_heatmap.png` | `lcr_signature_*_stacked.png` |\n"
        "| 3. Sequence-level properties | `master_lcr_continuous.csv`, `length_heatmap.png` | `lcr_length_histogram.png`, `lcr_fraction_histograms.png` |\n"
        "| 4. LCR multiplicity per protein | `headline_numbers.csv` (`pct_multi_lcr`, `mean_n_lcr`), `n_lcr_by_method.png` | `protein_multiplicity_stacked.png`, `protein_n_lcr_histogram.png` |\n"
        "| 5. Coverage and property fractions | `coverage_heatmap.png`, `master_protein_continuous.csv` | `protein_coverage_histogram.png`, `protein_signature_fraction_histograms.png` |\n\n"
        "## Conventions\n"
        "- Categorical enrichment: log2(observed/expected); 0 = expected, +1 = 2x enriched, -1 = 2x depleted.\n"
        "- Continuous: median, IQR, mean, SD, Kruskal-Wallis p (BH-FDR within level and method).\n"
        "- Robustness rule: a signal is robust if it holds at protein level; LCR level is supporting evidence.\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default=r"lcr_analyses\quantitative_analysis", help="Results root containing methods/")
    args = parser.parse_args()

    results_root = Path(args.outdir)
    summary_tables = results_root / "summary" / "tables"
    summary_figures = results_root / "summary" / "figures"
    summary_tables.mkdir(parents=True, exist_ok=True)
    summary_figures.mkdir(parents=True, exist_ok=True)

    data = collect(results_root)
    if not data:
        raise SystemExit(f"No per-method tables found under {results_root / 'methods'} — run the per-method analysis first.")

    # 1. Master stacked tables
    name_map = {
        "lcr_continuous": "master_lcr_continuous.csv",
        "protein_continuous": "master_protein_continuous.csv",
        "lcr_categorical": "master_lcr_categorical.csv",
        "protein_categorical": "master_protein_categorical.csv",
        "protein_core": "master_protein_core_metrics.csv",
        "qc": "master_qc_overview.csv",
    }
    for key, fname in name_map.items():
        if key in data:
            data[key].to_csv(summary_tables / fname, index=False)

    # 2. Headline numbers
    headline = build_headline(data)
    headline.to_csv(summary_tables / "headline_numbers.csv", index=False)

    # 3. Class x method pivots + heatmaps
    if "protein_continuous" in data:
        cov = pivot_metric(data["protein_continuous"], "coverage", "median")
        cov.to_csv(summary_tables / "median_coverage_by_class_method.csv")
        plot_heatmap(cov, "Median LCR coverage per protein", summary_figures / "coverage_heatmap.png")
    if "lcr_continuous" in data:
        length = pivot_metric(data["lcr_continuous"], "length", "median")
        length.to_csv(summary_tables / "median_length_by_class_method.csv")
        plot_heatmap(length, "Median LCR length (aa)", summary_figures / "length_heatmap.png", fmt=".0f")

    # 4. Signature enrichment, key signatures, class x method
    if "lcr_categorical" in data:
        cat = data["lcr_categorical"]
        sig = cat.loc[(cat["analysis"] == "signature_sequence_count_lcr") & (cat["signature"].isin(KEY_SIGNATURES))].copy()
        if not sig.empty:
            sig["class_method"] = sig["rna_class"] + " | " + sig["method"]
            mat = sig.pivot_table(index="signature", columns="class_method", values="log2_obs_exp", aggfunc="first")
            mat.to_csv(summary_tables / "signature_enrichment_key.csv")
            fig, ax = plt.subplots(figsize=(max(10, .55 * mat.shape[1] + 2), max(4, .5 * mat.shape[0] + 1.5)))
            vmax = max(.5, float(np.nanmax(np.abs(mat.to_numpy()))))
            sns.heatmap(mat, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=True, fmt=".1f", linewidths=.4, linecolor="white", ax=ax,
                        cbar_kws={"label": "log2(observed / expected)"})
            ax.set(title="Physicochemical-signature enrichment by RNA class and method", xlabel="RNA class | method", ylabel="Primary signature")
            ax.tick_params(axis="x", rotation=60)
            fig.tight_layout(); fig.savefig(summary_figures / "signature_enrichment_heatmap.png", dpi=220); plt.close(fig)
        top = build_top_signatures(cat)
        if not top.empty:
            top.to_csv(summary_tables / "top_signature_per_class.csv", index=False)

    # 5. Simple per-method comparison bars
    if not headline.empty:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        if "mean_n_lcr" in headline:
            axes[0].bar(headline["method"], headline["mean_n_lcr"], color="#87b7d9")
            axes[0].set(title="Mean LCRs per protein", ylabel="Mean n LCRs"); axes[0].tick_params(axis="x", rotation=40)
        if "median_coverage" in headline:
            axes[1].bar(headline["method"], headline["median_coverage"], color="#d9a787")
            axes[1].set(title="Median LCR coverage per protein", ylabel="Median coverage"); axes[1].tick_params(axis="x", rotation=40)
        fig.tight_layout(); fig.savefig(summary_figures / "n_lcr_and_coverage_by_method.png", dpi=220); plt.close(fig)

    # 6. Results guide
    write_results_guide(results_root)

    print(f"Summary tables:   {summary_tables.resolve()}")
    print(f"Summary figures:  {summary_figures.resolve()}")
    print(f"Results guide:    {(results_root / 'RESULTS_GUIDE.md').resolve()}")


if __name__ == "__main__":
    main()
