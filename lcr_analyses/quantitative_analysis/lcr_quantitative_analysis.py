#!/usr/bin/env python3
"""
Quantitative and diversity analysis of RBP LCRs.

Each detection method is analyzed and written independently:
  /methods/<method>/tables/
  /methods/<method>/figures/

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
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests

sns.set_theme(style="whitegrid", context="notebook")

# ---------------- configuration ----------------

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
MIN_N = 5  # n-floor: smaller groups are descriptive_only

# canonical column orders: identity -> counts -> rates -> effect -> stats -> flag
CAT_COLS = ["analysis", "feature_value", "rna_class",
            "n_pos_group", "n_tot_group", "n_pos_ref", "n_tot_ref",
            "rate_group", "rate_ref", "log2_or", "p_value", "p_fdr",
            "descriptive_only"]
CONT_COLS = ["rna_class", "metric", "n", "median", "mean", "sd",
             "cliffs_delta_vs_rest", "p_value", "p_fdr", "descriptive_only"]

# Canonical method labels.
METHOD_ALIASES = {
    "cast": "CAST",
    "seg": "SEG",
    "segintermediate": "SEG_intermediate",
    "flps": "FLPS",
    "lcrfinder": "LCRFinder",
    "alcor": "AlcoR",
}

# Canonical metric names.
METRIC_ALIASES = {
    "fracpolar": "frac_polar", "frachydrophobic": "frac_hydrophobic",
    "fracstronghydro": "frac_strong_hydro", "fracaromatic": "frac_aromatic",
    "fracdisorder": "frac_disorder", "fracpositive": "frac_positive",
    "fracnegative": "frac_negative", "fcr": "fcr", "ncpr": "ncpr",
    "fracgs": "frac_GS",
}
DISTRIBUTION_SUBMETRICS = {
    "nruns": "n_runs", "meanrunlength": "mean_run_length",
    "maxrunlength": "max_run_length", "meangap": "mean_gap",
    "cvgap": "cv_gap", "label": "label",
}
INTERVAL_KEYS = ["uniprot_accession", "method", "start", "end"]
RAW_KEY_COLUMNS = {"proteinid", "sourcemethod", "start", "end", "length"}

def order_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Apply the canonical column order; unlisted columns go last."""
    preferred = CAT_COLS if "log2_or" in df.columns else CONT_COLS
    cols = [c for c in preferred if c in df.columns]
    return df[cols + [c for c in df.columns if c not in cols]]

# ---------------- loading / standardization ----------------

def normalized_name(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).strip().lower())

def safe_name(x: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(x).strip()) or "unnamed_method"

def canonical_method(x: str) -> str:
    return METHOD_ALIASES.get(normalized_name(x), str(x).strip())

def canonical_metric_name(normalized: str) -> str:
    if re.fullmatch(r"frac[a-z]", normalized):
        return f"frac_{normalized[-1].upper()}"
    if normalized in METRIC_ALIASES:
        return METRIC_ALIASES[normalized]
    m = re.fullmatch(
        r"(polar|hydrophobic|aromatic|disorder|charged)"
        r"(nruns|meanrunlength|maxrunlength|meangap|cvgap|label)", normalized)
    if m:
        return f"{m.group(1)}_{DISTRIBUTION_SUBMETRICS[m.group(2)]}"
    return normalized

def extract_accession(x: str) -> str:
    """UniProt accession

    'Q9Y2T7|ensembl_protein_id=ENSP...|entry=...' -> 'Q9Y2T7'
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

def locate_column(df: pd.DataFrame, choices: Iterable[str], required: bool = True):
    available = {normalized_name(c): c for c in df.columns}
    for choice in choices:
        if normalized_name(choice) in available:
            return available[normalized_name(choice)]
    if required:
        raise KeyError(f"Expected one of {list(choices)}; found {list(df.columns)}")
    return None

def standardize_annotations(path: str) -> pd.DataFrame:
    raw = normalize_columns(pd.read_excel(path))
    mapping = {
        "proteinid": ["proteinid", "protein_id"],
        "method": ["sourcemethod", "source_method", "method"],
        "start": ["start"], "end": ["end"], "length": ["length"],
        "rna_class": ["rnatargetsuperclass", "rna_target_superclass"],
        "position": ["domainpositionclass", "domain_position_class"],
        "signature": ["primaryphysicochemicalannotation",
                      "primary_physicochemical_annotation"],
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
    """Merge all feature sheets on the interval keys; canonical metric names."""
    sheets = pd.read_excel(path, sheet_name=None)
    frames = []
    for sheet_name, raw in sheets.items():
        raw = normalize_columns(raw)
        base = pd.DataFrame({
            "uniprot_accession": raw[locate_column(
                raw, ["protein_id", "proteinid", "uniprot_accession"])].map(extract_accession),
            "method": raw[locate_column(
                raw, ["source_method", "sourcemethod", "method"])].map(canonical_method),
            "start": pd.to_numeric(raw[locate_column(raw, ["start"])],
                                   errors="coerce").astype("Int64"),
            "end": pd.to_numeric(raw[locate_column(raw, ["end"])],
                                 errors="coerce").astype("Int64"),
        })
        metrics = raw.drop(columns=[c for c in RAW_KEY_COLUMNS if c in raw.columns])
        metrics.columns = [canonical_metric_name(c) for c in metrics.columns]
        metrics = metrics.loc[:, ~metrics.columns.duplicated()]
        frames.append(pd.concat([base, metrics], axis=1))
    merged = frames[0]
    for frame in frames[1:]:
        merged = merged.merge(frame, on=INTERVAL_KEYS, how="outer")
    for col in merged.columns:
        if col not in INTERVAL_KEYS:
            try:
                merged[col] = pd.to_numeric(merged[col], errors="raise")
            except (ValueError, TypeError):
                pass
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
    """Join annotations + measurements + protein lengths; runtime QC printed."""
    ann = standardize_annotations(args.annotations)
    meas = standardize_measurements(args.measurements)
    master = standardize_master(args.master)

    metric_cols = [c for c in meas.columns if c not in INTERVAL_KEYS]
    dup = meas.duplicated(subset=INTERVAL_KEYS).sum()
    if dup:
        print(f"WARNING: {dup} duplicated measurement rows; keeping first.")
        meas = meas.drop_duplicates(subset=INTERVAL_KEYS)

    joined = ann.merge(meas[INTERVAL_KEYS + metric_cols], on=INTERVAL_KEYS,
                       how="left", validate="one_to_one")
    joined = joined.merge(master, on="uniprot_accession", how="left",
                          validate="many_to_one")
    joined["coverage_per_lcr"] = joined["length"] / joined["uniprot_length"]

    probe = "frac_polar" if "frac_polar" in joined.columns else metric_cols[0]
    print("\nJoin QC (rows with measurements / with UniProt length):")
    for method, d in joined.groupby("method"):
        print(f"  {method:<18} {d[probe].notna().mean():.1%} / "
              f"{d['uniprot_length'].notna().mean():.1%}")
    return joined

# ---------------- statistics ----------------

def bh_adjust(df: pd.DataFrame) -> pd.DataFrame:
    """Benjamini-Hochberg FDR within the given table (one analysis block)."""
    out = df.copy()
    out["p_fdr"] = np.nan
    ok = out["p_value"].notna()
    if ok.any():
        out.loc[ok, "p_fdr"] = multipletests(out.loc[ok, "p_value"],
                                             method="fdr_bh")[1]
    return out

def bh_adjust_per_block(df: pd.DataFrame, block_col: str) -> pd.DataFrame:
    """FDR scope = one analysis block, even when blocks share a sheet."""
    return (df.groupby(block_col, group_keys=False)[df.columns]
            .apply(bh_adjust))

def fisher_cell(a_pos: int, a_tot: int, b_pos: int, b_tot: int) -> dict:
    """One 2x2 cell: group vs. rest. Fisher p + Haldane log2 odds ratio.
    log2_or is NaN when either arm has zero carriers (Haldane sign artifact),
    matching compare_binary() in the background analysis."""
    a_neg, b_neg = a_tot - a_pos, b_tot - b_pos
    try:
        p_value = float(fisher_exact([[a_pos, a_neg], [b_pos, b_neg]]).pvalue)
    except (ValueError, ZeroDivisionError):
        p_value = np.nan
    if a_pos == 0 or b_pos == 0:
        log2_or = np.nan
    else:
        log2_or = float(np.log2(((a_pos + 0.5) * (b_neg + 0.5)) /
                                ((a_neg + 0.5) * (b_pos + 0.5))))
    return {
        "n_pos_group": a_pos, "n_tot_group": a_tot,
        "n_pos_ref": b_pos, "n_tot_ref": b_tot,
        "rate_group": a_pos / a_tot if a_tot else np.nan,
        "rate_ref": b_pos / b_tot if b_tot else np.nan,
        "log2_or": log2_or, "p_value": p_value,
        "descriptive_only": bool(min(a_pos, a_neg, b_pos, b_neg) < MIN_N),
    }


def categorical_enrichment(data: pd.DataFrame, feature: str,
                           weight: str | None = None,
                           report_values: list[str] | None = None) -> pd.DataFrame:
    """One-vs-rest Fisher per class x feature-value cell.

    Units are rows of `data` (LCRs or proteins); `weight` switches the basis
    from sequence count to residue count (sum of LCR lengths).
    `report_values` restricts the reported feature values -- used to drop
    mirror rows of binary features (e.g. keep only "no_pfam" or "present").
    """
    cols = [feature, "rna_class"] + ([weight] if weight else [])
    d = data[cols].dropna()
    values = report_values if report_values is not None else sorted(
        d[feature].astype(str).unique())
    rows = []
    for value in values:
        for cls in sorted(d["rna_class"].unique()):
            def count(sub):  # noqa: B023 -- late binding intended per value/cls
                return float(sub[weight].sum()) if weight else float(len(sub))
            in_class = d[d["rna_class"] == cls]
            in_value = d[d[feature].astype(str) == value]
            a_pos = int(round(count(in_class[in_class[feature].astype(str) == value])))
            a_tot = int(round(count(in_class)))
            b_pos = int(round(count(in_value[in_value["rna_class"] != cls])))
            b_tot = int(round(count(d[d["rna_class"] != cls])))
            r = fisher_cell(a_pos, a_tot, b_pos, b_tot)
            r.update({"feature_value": value, "rna_class": cls})
            rows.append(r)
    return pd.DataFrame(rows)

def continuous_summary(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Per class: n, median, mean, SD + Cliff's delta vs. rest + KW omnibus p."""
    x = df[["rna_class", metric]].dropna()
    groups = [g[metric].to_numpy() for _, g in x.groupby("rna_class") if len(g) >= 2]
    try:
        p_kw = float(mannwhitneyu(*groups).pvalue) if len(groups) >= 2 else np.nan
    except ValueError:
        p_kw = np.nan
    rows = []
    for cls, g in x.groupby("rna_class"):
        vals = g[metric]
        rest = x.loc[x["rna_class"] != cls, metric]
        delta = np.nan
        if len(vals) >= 2 and len(rest) >= 2:
            u = mannwhitneyu(vals, rest, alternative="two-sided")
            delta = float(2 * u.statistic / (len(vals) * len(rest)) - 1)  # type: ignore
        rows.append({
            "rna_class": cls, "metric": metric, "n": len(vals),
            "median": vals.median(), "mean": vals.mean(), "sd": vals.std(),
            "cliffs_delta_vs_rest": delta, "p_value": p_kw,
            "descriptive_only": bool(len(vals) < MIN_N),
        })
    return pd.DataFrame(rows)

# ---------------- figures (plot policy: histograms + stacked bars) ----------------

def plot_stacked_from_counts(data: pd.DataFrame, feature: str, title: str,
                             outfile: Path, weight: str | None = None) -> None:
    """Stacked proportion bars per class; proportions computed on the fly."""
    cols = [feature, "rna_class"] + ([weight] if weight else [])
    d = data[cols].dropna()
    if weight:
        table = pd.pivot_table(d, index="rna_class", columns=feature,
                               values=weight, aggfunc="sum", fill_value=0)
    else:
        table = pd.crosstab(d["rna_class"], d[feature])
    if table.empty:
        return
    mat = table.div(table.sum(axis=1), axis=0)
    ax = mat.plot(kind="bar", stacked=True,
                  figsize=(max(7, .9 * mat.shape[0]), 5), colormap="tab20")
    ax.set(title=title, xlabel="RNA-target superclass", ylabel="Proportion",
           ylim=(0, 1))
    ax.legend(title=feature, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(outfile, dpi=220)
    plt.close()

def plot_histogram(df: pd.DataFrame, metric: str, title: str, outfile: Path) -> None:
    x = df[["rna_class", metric]].dropna()
    if x.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.histplot(data=x, x=metric, hue="rna_class", element="step",
                 stat="density", common_norm=False, bins=35, ax=ax)
    ax.set(title=title, xlabel=metric, ylabel="Density")
    fig.tight_layout()
    fig.savefig(outfile, dpi=220)
    plt.close(fig)

def plot_fraction_panels(df, metrics, label_getter, title, outfile) -> None:
    available = [m for m in metrics if m in df.columns and df[m].notna().any()]
    if not available:
        return
    ncols = 3
    nrows = int(np.ceil(len(available) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 4 * nrows),
                             squeeze=False)
    axes_flat = axes.ravel()
    for ax, metric in zip(axes_flat, available):
        x = df[["rna_class", metric]].dropna()
        sns.histplot(data=x, x=metric, hue="rna_class", element="step",
                     stat="density", common_norm=False, bins=30, legend=False,
                     ax=ax)
        ax.set(title=label_getter(metric), xlabel=metric, ylabel="Density")
    for ax in axes_flat[len(available):]:
        ax.axis("off")
    handles, labels = axes_flat[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, title="RNA-target superclass",
                   loc="lower center", ncol=min(5, len(labels)),
                   bbox_to_anchor=(.5, .005))
    fig.suptitle(title, y=.99, fontsize=14)
    fig.tight_layout(rect=(0, .06, 1, .96))
    fig.savefig(outfile, dpi=220)
    plt.close(fig)

# ---------------- per-method analysis ----------------

def analyze_lcr_level(d: pd.DataFrame, figures: Path) -> dict[str, pd.DataFrame]:
    """LCR level (descriptive): categorical enrichment + continuous summaries."""
    method = d.name
    cat_frames = []
    pos_norm = d["position"].map(normalized_name)
    positional = d.loc[pos_norm != POSITION_NO_PFAM]
    cat = categorical_enrichment(positional, "position")
    cat["analysis"] = "position_lcr"
    cat_frames.append(cat)
    plot_stacked_from_counts(positional, "position",
                             f"Domain-position proportions (LCR): {method}",
                             figures / "lcr_position_stacked.png")

    # binary feature: report the no_pfam arm only (pfam_classified is its mirror)
    pfam = d.assign(has_pfam=np.where(pos_norm == POSITION_NO_PFAM,
                                     "no_pfam", "pfam_classified"))
    cat = categorical_enrichment(pfam, "has_pfam", report_values=["no_pfam"])
    cat["analysis"] = "no_pfam_lcr"
    cat_frames.append(cat)

    for basis, weight in [("sequence_count", None), ("residue_count", "length")]:
        cat = categorical_enrichment(d, "signature", weight)
        cat["analysis"] = f"signature_{basis}_lcr"
        cat_frames.append(cat)
        plot_stacked_from_counts(d, "signature",
                                 f"Primary-signature proportions ({basis}): {method}",
                                 figures / f"lcr_signature_{basis}_stacked.png",
                                 weight)

    cont_frames = []
    metrics = ["length", "coverage_per_lcr"] + [m for m in PROPERTY_METRICS
                                                if m in d.columns]
    for metric in metrics:
        cont_frames.append(continuous_summary(d, metric))

    for metric in SINGLE_HISTOGRAM_METRICS:
        if metric in d.columns:
            plot_histogram(d, metric, f"{metric} distribution (LCR): {method}",
                           figures / f"lcr_{metric}_histogram.png")
    plot_fraction_panels(d, FRACTION_PANEL_METRICS, lambda m: m,
                         f"LCR physicochemical fraction distributions: {method}",
                         figures / "lcr_fraction_histograms.png")
    return {"categorical_lcr": pd.concat(cat_frames, ignore_index=True),
            "continuous_lcr": pd.concat(cont_frames, ignore_index=True)}

def analyze_protein_level(d: pd.DataFrame, figures: Path) -> dict[str, pd.DataFrame]:
    """Protein level (inferential): multiplicity, presence, coverage, fractions."""
    method = d.name
    per_protein = (d.groupby(["uniprot_accession", "rna_class", "uniprot_length"],
                             dropna=False)
                   .agg(n_lcr=("length", "size"), lcr_residues=("length", "sum"))
                   .reset_index())
    per_protein["coverage"] = (per_protein["lcr_residues"]
                               / per_protein["uniprot_length"])
    per_protein["n_lcr_category"] = pd.cut(
        per_protein["n_lcr"], [0, 1, 2, np.inf], labels=["1", "2", "3+"]).astype(str)

    cat_frames = []
    cat = categorical_enrichment(per_protein, "n_lcr_category")
    cat["analysis"] = "multiplicity_protein"
    cat_frames.append(cat)
    plot_stacked_from_counts(per_protein, "n_lcr_category",
                             f"LCR multiplicity: {method}",
                             figures / "protein_multiplicity_stacked.png")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    sns.histplot(per_protein, x="n_lcr", discrete=True, ax=ax)
    ax.set(title=f"LCR count per protein: {method}", xlabel="Number of LCRs")
    fig.tight_layout()
    fig.savefig(figures / "protein_n_lcr_histogram.png", dpi=220)
    plt.close(fig)

    # presence tables: carrier vs non-carrier per feature value (report present only)
    base = d[["uniprot_accession", "rna_class"]].drop_duplicates()
    for feature_col, values, label in [
        ("position", [x for x in d["position"].dropna().unique()
                      if normalized_name(x) != POSITION_NO_PFAM], "position_presence"),
        ("signature", sorted(d["signature"].dropna().unique()),
         "signature_presence"),
    ]:
        for value in values:
            carriers = set(d.loc[d[feature_col] == value, "uniprot_accession"]) # type: ignore
            tbl = base.assign(present=base["uniprot_accession"].isin(carriers)
                              .map({True: "present", False: "absent"}))
            cat = categorical_enrichment(tbl, "present",
                                         report_values=["present"])
            cat["analysis"] = label
            cat["feature_value"] = value
            cat_frames.append(cat)

    plot_histogram(per_protein, "coverage",
                   f"LCR coverage distribution (protein): {method}",
                   figures / "protein_coverage_histogram.png")
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.scatterplot(data=per_protein, x="n_lcr", y="coverage", hue="rna_class",
                    alpha=.7, ax=ax)
    ax.set(title=f"Coverage vs LCR count: {method}")
    fig.tight_layout()
    fig.savefig(figures / "protein_coverage_vs_n_lcr.png", dpi=220)
    plt.close(fig)

    cont_frames = [continuous_summary(per_protein, "coverage")]
    top_signatures = d["signature"].value_counts().head(8).index.tolist()
    frac_input = d.assign(signature_for_fraction=np.where(
        d["signature"].isin(top_signatures), d["signature"], "other"))
    fractions = (frac_input.groupby(["uniprot_accession", "rna_class",
                                     "uniprot_length", "signature_for_fraction"],
                                    dropna=False)["length"].sum()
                 .rename("signature_residues").reset_index())
    fractions["property_fraction"] = (fractions["signature_residues"]
                                      / fractions["uniprot_length"])
    for signature, x in fractions.groupby("signature_for_fraction"):
        cont_frames.append(
            continuous_summary(x, "property_fraction")
            .assign(metric=f"property_fraction__{signature}"))
    pivoted = fractions.pivot_table(index=["uniprot_accession", "rna_class"],
                                    columns="signature_for_fraction",
                                    values="property_fraction", aggfunc="first")
    pivoted = pivoted.reset_index()
    plot_fraction_panels(
        pivoted, [c for c in pivoted.columns
                  if c not in {"uniprot_accession", "rna_class"}],
        lambda m: f"property_fraction: {m}",
        f"Protein-level signature fraction distributions: {method}",
        figures / "protein_signature_fraction_histograms.png")

    return {"categorical_protein": pd.concat(cat_frames, ignore_index=True),
            "continuous_protein": pd.concat(cont_frames, ignore_index=True)}

def qc_frame(d: pd.DataFrame, method: str) -> pd.DataFrame:
    overview = pd.DataFrame({
        "section": "overview", "method": [method], "n_lcr": [len(d)],
        "n_proteins": [d["uniprot_accession"].nunique()],
        "fraction_matched_uniprot_length": [d["uniprot_length"].notna().mean()],
    })
    per_class = (d.groupby("rna_class")["uniprot_accession"].nunique()
                 .rename("n_proteins").reset_index()
                 .assign(section="proteins_per_class"))
    return pd.concat([overview, per_class], ignore_index=True)

# ---------------- main ----------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--annotations", default="lcr_analyses/pc_properties/lcr_annotations.xlsx")
    parser.add_argument("--measurements", default="lcr_analyses/pc_properties/lcr_features.xlsx")
    parser.add_argument("--master", default="datasets/combined_rbp_pfam38_rbpdb_modomics_uniprot.xlsx")
    parser.add_argument("--outdir", default="lcr_analyses/results")
    args = parser.parse_args()

    out = Path(args.outdir)
    methods_root = out / "methods"
    methods_root.mkdir(parents=True, exist_ok=True)
    joined = read_and_join(args)

    manifest = []
    for method, d in joined.groupby("method", sort=True):
        d = d.copy()
        d.name = method  # type: ignore # used in figure titles
        tables = methods_root / safe_name(method) / "tables"  # type: ignore
        figures = methods_root / safe_name(method) / "figures"  # type: ignore
        tables.mkdir(parents=True, exist_ok=True)
        figures.mkdir(parents=True, exist_ok=True)
        for old_png in figures.glob("*.png"):
            old_png.unlink()

        sheets = {"qc": qc_frame(d, method)}  # type: ignore
        sheets.update(analyze_lcr_level(d, figures))
        sheets.update(analyze_protein_level(d, figures))

        # one workbook per method; FDR within each analysis block of a sheet
        workbook = tables / f"{safe_name(method)}_results.xlsx"  # type: ignore
        with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
            for name, df in sheets.items():
                if name == "qc":
                    df.to_excel(writer, sheet_name=name, index=False)
                elif "categorical" in name:
                    order_cols(bh_adjust_per_block(df, "analysis")).to_excel(
                        writer, sheet_name=name, index=False)
                else:
                    order_cols(bh_adjust(df)).to_excel(
                        writer, sheet_name=name, index=False)
        manifest.append({"method": method, "workbook": str(workbook),
                         "n_lcr": len(d),
                         "n_proteins": d["uniprot_accession"].nunique()})
        print(f"Completed {method}: {workbook}")

    pd.DataFrame(manifest).to_csv(out / "method_manifest.csv", index=False)
    print(f"Done. Results under: {methods_root.resolve()}")

if __name__ == "__main__":
    main()
