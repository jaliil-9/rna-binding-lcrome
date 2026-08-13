#!/usr/bin/env python3
"""Phase 2.4.1: quantitative and diversity analysis of RBP LCRs.

Each detection method is analyzed and written independently:
  <outdir>/methods/<method>/tables/
  <outdir>/methods/<method>/figures/

Plot policy: histograms only (no boxplots). The seven physicochemical fraction
histograms are grouped into one multi-panel figure per method
(lcr_fraction_histograms.png); protein-level signature fractions into
protein_signature_fraction_histograms.png.

Data-integration fixes (validated by debug_lcr_join.py):
  1. Accessions are extracted from pipe-delimited composite IDs
     (text before the first '|'), e.g.
     'Q9Y2T7|ensembl_protein_id=...' -> 'Q9Y2T7'.
  2. Method labels are canonicalized on both sides:
     'fLPS_strict' -> 'FLPS', 'SEG_strict' -> 'SEG'.
  3. ALL sheets of the measurements workbook (composition, distribution,
     co_occurrence) are read and merged on the interval keys, with metric
     columns renamed back to their canonical names (frac_A, frac_polar, ...).

No method, length, or overlap filters are applied. RNA-target superclass comes
from the annotations table; the master Combined sheet supplies uniprot_length.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import chi2_contingency, kruskal
from statsmodels.stats.multitest import multipletests

sns.set_theme(style="whitegrid", context="notebook")

AA = list("ACDEFGHIKLMNPQRSTVWY")
PROPERTY_METRICS = [
    "frac_polar", "frac_hydrophobic", "frac_strong_hydro", "frac_aromatic",
    "frac_disorder", "frac_positive", "frac_negative", "fcr", "ncpr", "frac_GS",
]
FRACTION_PANEL_METRICS = [
    "frac_GS", "frac_aromatic", "frac_disorder", "frac_hydrophobic",
    "frac_negative", "frac_polar", "frac_positive",
]
SINGLE_HISTOGRAM_METRICS = ["length", "coverage_per_lcr", "fcr", "ncpr"]
POSITION_NO_PFAM = "unclassifiednopfam"

# Canonical method labels: normalized (lowercase, alnum-only) -> canonical.
METHOD_ALIASES = {
    "cast": "CAST",
    "seg": "SEG",
    "segstrict": "SEG",
    "segintermediate": "SEG_intermediate",
    "flps": "FLPS",
    "flpsstrict": "FLPS",
    "lcrfinder": "LCRFinder",
    "alcor": "AlcoR",
}

# Canonical metric names: normalized (lowercase, alnum-only) -> canonical.
METRIC_ALIASES = {
    "fracpolar": "frac_polar", "frachydrophobic": "frac_hydrophobic",
    "fracstronghydro": "frac_strong_hydro", "fracaromatic": "frac_aromatic",
    "fracdisorder": "frac_disorder", "fracpositive": "frac_positive",
    "fracnegative": "frac_negative", "fcr": "fcr", "ncpr": "ncpr",
    "fracgs": "frac_GS",
    "nrgmotifs": "n_RG_motifs", "rgrepeatcluster": "RG_repeat_cluster", "rgrepeatregion": "RG_repeat_region",
    "nsrmotifs": "n_SR_motifs", "srrepeatcluster": "SR_repeat_cluster", "srrepeatregion": "SR_repeat_region",
    "ngsmotifs": "n_GS_motifs", "gsrepeatcluster": "GS_repeat_cluster", "gsrepeatregion": "GS_repeat_region",
    "coocpositivearomatic": "cooc_positive_aromatic",
}

DISTRIBUTION_SUBMETRICS = {
    "nruns": "n_runs", "meanrunlength": "mean_run_length",
    "maxrunlength": "max_run_length", "meangap": "mean_gap",
    "cvgap": "cv_gap", "label": "label",
}

INTERVAL_KEYS = ["uniprot_accession", "method", "start", "end"]
RAW_KEY_COLUMNS = {"proteinid", "sourcemethod", "start", "end", "length"}


def normalized_name(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).strip().lower())


def safe_name(x: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(x).strip()) or "unnamed_method"


def canonical_method(x: str) -> str:
    """Map any method spelling to the canonical label used in the analysis."""
    return METHOD_ALIASES.get(normalized_name(x), str(x).strip())


def canonical_metric_name(normalized: str) -> str:
    """Map a normalized measurements column name back to its canonical form."""
    if re.fullmatch(r"frac[a-z]", normalized):
        return f"frac_{normalized[-1].upper()}"
    if normalized in METRIC_ALIASES:
        return METRIC_ALIASES[normalized]
    m = re.fullmatch(r"(polar|hydrophobic|aromatic|disorder|charged)(nruns|meanrunlength|maxrunlength|meangap|cvgap|label)", normalized)
    if m:
        return f"{m.group(1)}_{DISTRIBUTION_SUBMETRICS[m.group(2)]}"
    return normalized


def extract_accession(x: str) -> str:
    """UniProt accession from a pipe-delimited composite protein ID.

    'Q9Y2T7|ensembl_protein_id=ENSP...|entry=...' -> 'Q9Y2T7'
    Also tolerates the legacy concatenated format and bare accessions.
    """
    text = str(x).strip()
    if "|" in text:
        return text.split("|", 1)[0].strip()
    m = re.match(r"^([A-Za-z0-9]+?)(?:ensemblproteinid|$)", text)
    return m.group(1) if m else text


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [normalized_name(c) for c in out.columns]
    return out


def locate_column(df: pd.DataFrame, choices: Iterable[str], required: bool = True) -> str | None:
    available = {normalized_name(c): c for c in df.columns}
    for choice in choices:
        if normalized_name(choice) in available:
            return available[normalized_name(choice)]
    if required:
        raise KeyError(f"Required column missing. Expected one of {list(choices)}; found {list(df.columns)}")
    return None


def standardize_annotations(path: str) -> pd.DataFrame:
    raw = normalize_columns(pd.read_excel(path))
    mapping = {
        "proteinid": ["proteinid", "protein_id"],
        "method": ["sourcemethod", "source_method", "method"],
        "start": ["start"], "end": ["end"], "length": ["length"],
        "rna_class": ["rnatargetsuperclass", "rna_target_superclass"],
        "position": ["domainpositionclass", "domain_position_class"],
        "signature": ["primaryphysicochemicalannotation", "primary_physicochemical_annotation"],
    }
    out = pd.DataFrame({new: raw[locate_column(raw, old)] for new, old in mapping.items()})
    out["uniprot_accession"] = out["proteinid"].map(extract_accession)
    for col in ["method", "rna_class", "position", "signature"]:
        out[col] = out[col].fillna("missing").astype(str).str.strip()
    out["method"] = out["method"].map(canonical_method)
    for col in ["start", "end", "length"]:
        out[col] = pd.to_numeric(out[col], errors="coerce").astype("Int64")
    return out.dropna(subset=["start", "end", "length"]).copy()


def standardize_measurements(path: str) -> pd.DataFrame:
    """Read ALL sheets and merge them on the interval keys into one wide table.

    Metric columns are renamed to canonical form (frac_A, frac_polar, ...);
    the raw key columns of each sheet are used only for joining and dropped.
    """
    sheets = pd.read_excel(path, sheet_name=None)
    frames = []
    for sheet_name, raw in sheets.items():
        raw = normalize_columns(raw)
        base = pd.DataFrame({
            "measurement_protein_id": raw[locate_column(raw, ["protein_id", "proteinid", "uniprot_accession", "uniprot"])],
            "method": raw[locate_column(raw, ["source_method", "sourcemethod", "method"])],
            "start": raw[locate_column(raw, ["start"])],
            "end": raw[locate_column(raw, ["end"])],
        })
        metrics = raw.drop(columns=[c for c in RAW_KEY_COLUMNS if c in raw.columns])
        metrics.columns = [canonical_metric_name(c) for c in metrics.columns]
        metrics = metrics.loc[:, ~metrics.columns.duplicated()]
        frame = pd.concat([base, metrics], axis=1)
        frame["_sheet"] = sheet_name
        frames.append(frame)
    merged = frames[0]
    for frame in frames[1:]:
        merged = merged.merge(frame, on=["measurement_protein_id", "method", "start", "end"], how="outer", suffixes=("", "_dup"))
        merged = merged.loc[:, ~merged.columns.str.endswith("_dup")]
    merged["measurement_protein_id"] = merged["measurement_protein_id"].astype(str).str.strip()
    merged["uniprot_accession"] = merged["measurement_protein_id"].map(extract_accession)
    merged["method"] = merged["method"].astype(str).str.strip().map(canonical_method)
    for col in ["start", "end"]:
        merged[col] = pd.to_numeric(merged[col], errors="coerce").astype("Int64")
    for col in merged.columns:
        if col not in {"measurement_protein_id", "uniprot_accession", "method", "_sheet"}:
            try:
                merged[col] = pd.to_numeric(merged[col], errors="raise")
            except (ValueError, TypeError):
                pass  # genuinely textual column (e.g. distribution labels) — keep as-is
    return merged


def standardize_master(path: str) -> pd.DataFrame:
    raw = normalize_columns(pd.read_excel(path, sheet_name="Combined"))
    acc = locate_column(raw, ["uniprot_accession", "uniprotaccession"])
    length = locate_column(raw, ["uniprot_length", "uniprotlength"])
    out = raw[[acc, length]].copy()
    out.columns = ["uniprot_accession", "uniprot_length"]
    out["uniprot_accession"] = out["uniprot_accession"].astype(str).str.strip()
    out["uniprot_length"] = pd.to_numeric(out["uniprot_length"], errors="coerce")
    return out.drop_duplicates("uniprot_accession")


def read_and_join(args: argparse.Namespace) -> pd.DataFrame:
    ann = standardize_annotations(args.annotations)
    meas = standardize_measurements(args.measurements)
    master = standardize_master(args.master)

    metric_cols = [c for c in meas.columns if c not in {"measurement_protein_id", "_sheet"} | set(INTERVAL_KEYS)]
    dup = meas.duplicated(subset=INTERVAL_KEYS).sum()
    if dup:
        print(f"WARNING: {dup} duplicated measurement rows on interval keys; keeping first occurrence.")
        meas = meas.drop_duplicates(subset=INTERVAL_KEYS)

    joined = ann.merge(meas[INTERVAL_KEYS + metric_cols], on=INTERVAL_KEYS, how="left", validate="one_to_one")
    joined = joined.merge(master, on="uniprot_accession", how="left", validate="many_to_one")
    joined["coverage_per_lcr"] = joined["length"] / joined["uniprot_length"]

    # Runtime join-rate QC so silent mismatches cannot recur.
    probe = "frac_polar" if "frac_polar" in joined.columns else metric_cols[0]
    print("\nJoin QC (fraction of annotation rows with measurements / with UniProt length):")
    for method, d in joined.groupby("method"):
        print(f"  {method:<18} measurements: {d[probe].notna().mean():.1%}   uniprot_length: {d['uniprot_length'].notna().mean():.1%}")
    overall = joined["uniprot_length"].notna().mean()
    if overall < 0.95:
        print(f"WARNING: only {overall:.1%} of rows matched a UniProt length — check accession extraction.")
    return joined


def bh_adjust(df: pd.DataFrame, pcol: str = "p_value") -> pd.DataFrame:
    out = df.copy()
    out["p_fdr"] = np.nan
    ok = out[pcol].notna()
    if ok.any():
        out.loc[ok, "p_fdr"] = multipletests(out.loc[ok, pcol], method="fdr_bh")[1]
    return out


def safe_chi_square(table: pd.DataFrame) -> float:
    table = table.loc[table.sum(axis=1) > 0, table.sum(axis=0) > 0]
    if table.shape[0] < 2 or table.shape[1] < 2:
        return np.nan
    try:
        return float(chi2_contingency(table.astype(float), correction=False).pvalue)  # type: ignore
    except ValueError:
        return np.nan


def categorical_enrichment(df: pd.DataFrame, feature: str, weight: str | None = None) -> tuple[pd.DataFrame, float]:
    cols = [feature, "rna_class"] + ([weight] if weight else [])
    data = df[cols].dropna().copy()
    if weight:
        table = pd.pivot_table(data, index=feature, columns="rna_class", values=weight, aggfunc="sum", fill_value=0)
    else:
        table = pd.crosstab(data[feature], data["rna_class"])
    # Force float dtype: nullable Int64 weights yield object-dtype arrays that
    # break NumPy's in-place divide/log2 below.
    observed = table.to_numpy(dtype=float)
    total = observed.sum()
    expected = np.outer(table.sum(axis=1).to_numpy(dtype=float), table.sum(axis=0).to_numpy(dtype=float)) / total
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.divide(observed, expected, out=np.full(expected.shape, np.nan), where=expected > 0)
        fold = np.log2(ratio, out=np.full(expected.shape, np.nan), where=ratio > 0)
    long = table.rename_axis(index=feature, columns="rna_class").stack().rename("observed").reset_index()  # type: ignore
    long["observed"] = observed.ravel()
    long["expected"] = expected.ravel()
    long["log2_obs_exp"] = fold.ravel()
    long["proportion_within_rna_class"] = long["observed"] / long.groupby("rna_class")["observed"].transform("sum")
    return long, safe_chi_square(table)


def continuous_summary(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    x = df[["rna_class", metric]].dropna().copy()
    summary = x.groupby("rna_class")[metric].agg(
        n="size", median="median", mean="mean", sd="std",
        q1=lambda s: s.quantile(.25), q3=lambda s: s.quantile(.75),
    ).reset_index()
    summary["iqr"] = summary["q3"] - summary["q1"]
    groups = [g[metric].to_numpy() for _, g in x.groupby("rna_class") if len(g) >= 2]
    try:
        p = float(kruskal(*groups).pvalue) if len(groups) >= 2 else np.nan
    except ValueError:
        p = np.nan
    summary["p_value"] = p
    return summary


def plot_stacked(long: pd.DataFrame, feature: str, title: str, outfile: Path) -> None:
    mat = long.pivot(index="rna_class", columns=feature, values="proportion_within_rna_class").fillna(0)
    if mat.empty:
        return
    ax = mat.plot(kind="bar", stacked=True, figsize=(max(7, .9 * mat.shape[0]), 5), colormap="tab20")
    ax.set(title=title, xlabel="RNA-target superclass", ylabel="Proportion", ylim=(0, 1))
    ax.legend(title=feature, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout(); plt.savefig(outfile, dpi=220); plt.close()


def plot_histogram(df: pd.DataFrame, metric: str, title: str, outfile: Path) -> None:
    """Single density histogram per RNA class. No boxplots are produced anywhere."""
    x = df[["rna_class", metric]].dropna()
    if x.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.histplot(data=x, x=metric, hue="rna_class", element="step", stat="density", common_norm=False, bins=35, ax=ax)
    ax.set(title=title, xlabel=metric, ylabel="Density")
    fig.tight_layout(); fig.savefig(outfile, dpi=220); plt.close(fig)


def plot_fraction_panels(df: pd.DataFrame, metrics: list[str], label_getter, title: str, outfile: Path) -> None:
    """One multi-panel figure with a separate histogram per metric."""
    available = [m for m in metrics if m in df.columns and df[m].notna().any()]
    if not available:
        return
    ncols = 3
    nrows = int(np.ceil(len(available) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 4 * nrows), squeeze=False)
    axes_flat = axes.ravel()
    for ax, metric in zip(axes_flat, available):
        x = df[["rna_class", metric]].dropna()
        sns.histplot(data=x, x=metric, hue="rna_class", element="step", stat="density", common_norm=False, bins=30, legend=False, ax=ax)
        ax.set(title=label_getter(metric, x), xlabel=metric, ylabel="Density")
    for ax in axes_flat[len(available):]:
        ax.axis("off")
    handles, labels = axes_flat[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, title="RNA-target superclass", loc="lower center", ncol=min(5, len(labels)), bbox_to_anchor=(.5, .005))
    fig.suptitle(title, y=.99, fontsize=14)
    fig.tight_layout(rect=(0, .06, 1, .96))
    fig.savefig(outfile, dpi=220)
    plt.close(fig)


def write_qc(d: pd.DataFrame, tables: Path, method: str) -> None:
    qc = pd.DataFrame({
        "method": [method], "n_lcr": [len(d)], "n_proteins": [d["uniprot_accession"].nunique()],
        "fraction_matched_uniprot_length": [d["uniprot_length"].notna().mean()],
    })
    qc.to_csv(tables / "qc_overview.csv", index=False)
    d.groupby("rna_class")["uniprot_accession"].nunique().rename("n_proteins").reset_index().to_csv(tables / "qc_proteins_per_rna_class.csv", index=False)


def analyze_lcr_categorical(d: pd.DataFrame, method: str, tables: Path, figures: Path) -> None:
    all_results = []
    positional = d.loc[d["position"] != POSITION_NO_PFAM]
    long, p = categorical_enrichment(positional, "position")
    long["analysis"] = "position_lcr"; long["p_value"] = p; all_results.append(long)
    plot_stacked(long, "position", f"Domain-position proportions (LCR): {method}", figures / "lcr_position_stacked.png")

    pfam = d.assign(has_pfam=np.where(d["position"] == POSITION_NO_PFAM, "no_pfam", "pfam_classified"))
    long, p = categorical_enrichment(pfam, "has_pfam")
    long["analysis"] = "no_pfam_lcr"; long["p_value"] = p; all_results.append(long)

    for basis, weight in [("sequence_count", None), ("residue_count", "length")]:
        long, p = categorical_enrichment(d, "signature", weight)
        long["analysis"] = f"signature_{basis}_lcr"; long["p_value"] = p; all_results.append(long)
        plot_stacked(long, "signature", f"Primary-signature proportions ({basis}): {method}", figures / f"lcr_signature_{basis}_stacked.png")

    result = pd.concat(all_results, ignore_index=True)
    result["method"] = method
    bh_adjust(result).to_csv(tables / "lcr_level_categorical_enrichment.csv", index=False)


def analyze_lcr_continuous(d: pd.DataFrame, method: str, tables: Path, figures: Path) -> None:
    metrics = ["length", "coverage_per_lcr"] + [m for m in PROPERTY_METRICS if m in d.columns]
    summaries = []
    for metric in metrics:
        summary = continuous_summary(d, metric).assign(method=method, level="lcr", metric=metric)
        summaries.append(summary)
    result = pd.concat(summaries, ignore_index=True)
    bh_adjust(result).to_csv(tables / "lcr_level_continuous_summary.csv", index=False)

    for metric in SINGLE_HISTOGRAM_METRICS:
        if metric in d.columns:
            plot_histogram(d, metric, f"{metric} distribution (LCR): {method}", figures / f"lcr_{metric}_histogram.png")
    plot_fraction_panels(
        d, FRACTION_PANEL_METRICS, lambda m, x: m,
        f"LCR physicochemical fraction distributions: {method}",
        figures / "lcr_fraction_histograms.png",
    )


def protein_presence_table(d: pd.DataFrame, feature_col: str, feature_value: str) -> tuple[pd.DataFrame, float]:
    base = d[["uniprot_accession", "rna_class"]].drop_duplicates()
    present = d.loc[d[feature_col] == feature_value, "uniprot_accession"].drop_duplicates()
    base["present"] = np.where(base["uniprot_accession"].isin(present), "present", "absent")
    return categorical_enrichment(base, "present")


def analyze_protein(d: pd.DataFrame, method: str, tables: Path, figures: Path) -> None:
    per_protein = d.groupby(["uniprot_accession", "rna_class", "uniprot_length"], dropna=False).agg(n_lcr=("length", "size"), lcr_residues=("length", "sum")).reset_index()
    per_protein["coverage"] = per_protein["lcr_residues"] / per_protein["uniprot_length"]
    per_protein["n_lcr_category"] = pd.cut(per_protein["n_lcr"], [0, 1, 2, np.inf], labels=["1", "2", "3+"]).astype(str)
    per_protein["method"] = method
    per_protein.to_csv(tables / "protein_level_core_metrics.csv", index=False)

    categorical = []
    long, p = categorical_enrichment(per_protein, "n_lcr_category")
    long["analysis"] = "multiplicity_protein"; long["p_value"] = p; categorical.append(long)
    plot_stacked(long, "n_lcr_category", f"LCR multiplicity: {method}", figures / "protein_multiplicity_stacked.png")
    fig, ax = plt.subplots(figsize=(7, 4.5)); sns.histplot(per_protein, x="n_lcr", discrete=True, ax=ax); ax.set(title=f"LCR count per protein: {method}", xlabel="Number of LCRs")
    fig.tight_layout(); fig.savefig(figures / "protein_n_lcr_histogram.png", dpi=220); plt.close(fig)

    for feature_col, values, label in [
        ("position", [x for x in d["position"].dropna().unique() if x != POSITION_NO_PFAM], "position_presence"),
        ("signature", d["signature"].dropna().unique(), "signature_presence"),
    ]:
        for value in values:
            long, p = protein_presence_table(d, feature_col, value)
            long = long.loc[long["present"] == "present"].copy()
            long["analysis"] = label; long["feature_value"] = value; long["p_value"] = p; categorical.append(long)
    cat = pd.concat(categorical, ignore_index=True); cat["method"] = method
    bh_adjust(cat).to_csv(tables / "protein_level_categorical_enrichment.csv", index=False)

    plot_histogram(per_protein, "coverage", f"LCR coverage distribution (protein): {method}", figures / "protein_coverage_histogram.png")
    fig, ax = plt.subplots(figsize=(7, 5)); sns.scatterplot(data=per_protein, x="n_lcr", y="coverage", hue="rna_class", alpha=.7, ax=ax)
    ax.set(title=f"Coverage vs LCR count: {method}")
    fig.tight_layout(); fig.savefig(figures / "protein_coverage_vs_n_lcr.png", dpi=220); plt.close(fig)

    continuous = [continuous_summary(per_protein, "coverage").assign(method=method, metric="coverage", level="protein")]
    top_signatures = d["signature"].value_counts().head(8).index.tolist()
    fraction_input = d.assign(signature_for_fraction=np.where(d["signature"].isin(top_signatures), d["signature"], "other"))
    fractions = fraction_input.groupby(["uniprot_accession", "rna_class", "uniprot_length", "signature_for_fraction"], dropna=False)["length"].sum().rename("signature_residues").reset_index()
    fractions["property_fraction"] = fractions["signature_residues"] / fractions["uniprot_length"]
    fractions["method"] = method
    fractions.to_csv(tables / "protein_level_signature_fractions.csv", index=False)
    for signature, x in fractions.groupby("signature_for_fraction"):
        continuous.append(continuous_summary(x, "property_fraction").assign(method=method, metric=f"property_fraction__{signature}", level="protein"))
    pivoted = fractions.pivot_table(index=["uniprot_accession", "rna_class"], columns="signature_for_fraction", values="property_fraction", aggfunc="first").reset_index()
    signature_metrics = [c for c in pivoted.columns if c not in {"uniprot_accession", "rna_class"}]
    plot_fraction_panels(
        pivoted, signature_metrics, lambda m, x: f"property_fraction: {m}",
        f"Protein-level signature fraction distributions: {method}",
        figures / "protein_signature_fraction_histograms.png",
    )
    result = pd.concat(continuous, ignore_index=True)
    bh_adjust(result).to_csv(tables / "protein_level_continuous_summary.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--annotations", default=r"lcr_analyses\pc_properties\lcr_annotations.xlsx")
    parser.add_argument("--measurements", default=r"lcr_analyses\pc_properties\lcr_features.xlsx")
    parser.add_argument("--master", default=r"datasets\combined_rbp_pfam38_rbpdb_modomics_uniprot.xlsx")
    parser.add_argument("--outdir", default=r"lcr_analyses\results")
    args = parser.parse_args()

    out = Path(args.outdir)
    methods_root = out / "methods"
    methods_root.mkdir(parents=True, exist_ok=True)
    joined = read_and_join(args)

    manifest = []
    for method, d in joined.groupby("method", sort=True):
        method_dir = methods_root / safe_name(method)  # type: ignore
        tables = method_dir / "tables"
        figures = method_dir / "figures"
        tables.mkdir(parents=True, exist_ok=True)
        figures.mkdir(parents=True, exist_ok=True)
        for old_png in figures.glob("*.png"):
            old_png.unlink()
        d = d.copy()
        d.to_csv(tables / "joined_input_used_for_analysis.csv", index=False)
        write_qc(d, tables, method)  # type: ignore
        analyze_lcr_categorical(d, method, tables, figures)  # type: ignore
        analyze_lcr_continuous(d, method, tables, figures)  # type: ignore
        analyze_protein(d, method, tables, figures)  # type: ignore
        manifest.append({"method": method, "folder": str(method_dir), "n_lcr": len(d), "n_proteins": d["uniprot_accession"].nunique()})
        print(f"Completed {method}: {method_dir}")

    pd.DataFrame(manifest).to_csv(out / "method_manifest.csv", index=False)
    print(f"All method-specific results written under: {methods_root.resolve()}")


if __name__ == "__main__":
    main()
