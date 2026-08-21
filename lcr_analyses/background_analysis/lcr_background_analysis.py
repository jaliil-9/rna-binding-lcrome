#!/usr/bin/env python3
"""Phase 2.5 (v2.2): RBP LCRs vs. the proteome background.

Pipeline: load universe -> join LCRs -> per method { Mode A (RBP vs background),
Mode B (per RNA class vs background) } -> one Excel workbook per method ->
summary: headline numbers + the Phase 2.4.1 validation grid.
Output manual: results_guide_phase2_5.md.

Design invariants:
- Reference: the non-RBP background (rung 2). Rung 1 (full proteome) is
  descriptive only -- the proteome contains the RBPs and cannot be tested
  against itself.
- Statistics: Fisher exact + log2 odds ratio (binary); Mann-Whitney +
  Cliff's delta (continuous). BH-FDR within block x method; alpha = 0.05;
  effect sizes lead; n < 5 -> descriptive_only. When either arm has zero
  carriers the log2 OR is uninformative (Haldane sign artifact) -> NaN;
  the counts still describe the row.
- Length confounding: incidence also reported within background-length
  quintiles.
- Signature sheets are split by measure (rates vs. occupancy); empty
  columns and zero-carrier occupancy rows are dropped per sheet.
- Column order fixed by order_cols(): identity -> counts -> rates ->
  effect size -> statistics -> flags.
- Validation grid: internal (2.4.1, class vs. other RBPs) and external
  (2.5, class vs. background) log2 OR side by side.

v2.2 review fixes:
- load_features now restores canonical metric names (frac_A, frac_polar, ...)
  -- v1/v2.1 normalized names never mapped back, so composition blocks
  silently covered only fcr/ncpr/length.
- Sheets drop all-empty columns; occupancy rows with zero carriers dropped.
- Zero-carrier binary rows: log2_or = NaN (counts still reported).

Inputs
------
--annotations  ensembl_lcr_annotations.xlsx  (proteome-run annotation output)
--features     ensembl_lcr_features.xlsx     (composition/distribution/co_occurrence)
--background   finalized_background.xlsx     (accession, sequence, length)
--labels       finalized_background_labels.tsv
--classes      rbp_rna_classification.xlsx   (per-class sheets)
--patch        finalized_label_patch.tsv
--internal-dir Phase 2.4.1 v2 methods dir (per-method *_results.xlsx)
--outdir       results dir (methods/ + summary/ created below)
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests

AA = list("ACDEFGHIKLMNPQRSTVWY")
PROPERTY_METRICS = [
    "frac_polar", "frac_hydrophobic", "frac_strong_hydro", "frac_aromatic",
    "frac_disorder", "frac_positive", "frac_negative", "fcr", "ncpr", "frac_GS",
]
COMPOSITION_METRICS = [f"frac_{a}" for a in AA] + PROPERTY_METRICS
ARCHITECTURE_METRICS = ["coverage", "n_lcr", "mean_lcr_length"]
MIN_N = 5

METHOD_ALIASES = {
    "cast": "CAST", "seg": "SEG", "segstrict": "SEG",
    "segintermediate": "SEG_intermediate", "flps": "FLPS", "flpsstrict": "FLPS",
    "lcrfinder": "LCRFinder", "alcor": "AlcoR",
}
# normalized (lowercase, alnum-only) -> canonical metric name
METRIC_ALIASES = {
    "fracpolar": "frac_polar", "frachydrophobic": "frac_hydrophobic",
    "fracstronghydro": "frac_strong_hydro", "fracaromatic": "frac_aromatic",
    "fracdisorder": "frac_disorder", "fracpositive": "frac_positive",
    "fracnegative": "frac_negative", "fcr": "fcr", "ncpr": "ncpr",
    "fracgs": "frac_GS",
    "nrgmotifs": "n_RG_motifs", "nsrmotifs": "n_SR_motifs",
    "ngsmotifs": "n_GS_motifs", "coocpositivearomatic": "cooc_positive_aromatic",
}
RAW_KEY_COLUMNS = {"proteinid", "sourcemethod", "start", "end", "length"}

BINARY_COLS = ["block", "group", "signature", "length_bin",
               "n_pos_group", "n_tot_group", "n_pos_ref", "n_tot_ref",
               "rate_group", "rate_ref", "log2_or", "p_value", "p_fdr",
               "descriptive_only"]
CONT_COLS = ["block", "group", "metric", "signature",
             "n_group", "n_ref", "median_group", "median_ref",
             "mean_group", "mean_ref", "sd_group", "sd_ref",
             "cliffs_delta", "p_value", "p_fdr", "descriptive_only"]

def order_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Canonical column order; applied after empty columns are dropped."""
    preferred = BINARY_COLS if "log2_or" in df.columns else CONT_COLS
    cols = [c for c in preferred if c in df.columns]
    return df[cols + [c for c in df.columns if c not in cols]]

# ---------------- small helpers ----------------

def normalized_name(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).strip().lower())

def canonical_method(x: str) -> str:
    return METHOD_ALIASES.get(normalized_name(x), str(x).strip())

def canonical_metric_name(normalized: str) -> str:
    """Map a normalized feature column name back to canonical form."""
    if re.fullmatch(r"frac[a-z]", normalized):
        return f"frac_{normalized[-1].upper()}"
    return METRIC_ALIASES.get(normalized, normalized)

def extract_accession(x: str) -> str:
    text = str(x).strip()
    if "|" in text:
        parts = text.split("|")
        if parts[0] in {"sp", "tr"} and len(parts) > 1:
            return parts[1].strip()
        return parts[0].strip()
    return text

def bh_adjust(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["p_fdr"] = np.nan
    ok = out["p_value"].notna()
    if ok.any():
        out.loc[ok, "p_fdr"] = multipletests(out.loc[ok, "p_value"],
                                             method="fdr_bh")[1]
    return out

def bh_adjust_per_block(df: pd.DataFrame, block_col: str = "block") -> pd.DataFrame:
    """FDR scope = one analysis block, even when blocks share a sheet."""
    return df.groupby(block_col, group_keys=False)[df.columns].apply(bh_adjust)

# ---------------- statistics ----------------

def compare_binary(a_pos: int, a_tot: int, b_pos: int, b_tot: int) -> dict:
    """2x2 carrier table: group vs. reference. Fisher p + Haldane log2 OR.
    log2 OR is NaN when either arm has zero carriers (sign is an artifact
    of unequal group sizes under the Haldane correction)."""
    a_neg, b_neg = a_tot - a_pos, b_tot - b_pos
    try:
        p = float(fisher_exact([[a_pos, a_neg], [b_pos, b_neg]])[1])  # type: ignore
    except ValueError:
        p = np.nan
    if a_pos == 0 or b_pos == 0:
        log2_or = np.nan
    else:
        or_h = ((a_pos + 0.5) * (b_neg + 0.5)) / ((b_pos + 0.5) * (a_neg + 0.5))
        log2_or = float(np.log2(or_h))
    return {
        "n_pos_group": a_pos, "n_tot_group": a_tot,
        "n_pos_ref": b_pos, "n_tot_ref": b_tot,
        "rate_group": a_pos / a_tot if a_tot else np.nan,
        "rate_ref": b_pos / b_tot if b_tot else np.nan,
        "log2_or": log2_or, "p_value": p,
        "descriptive_only": bool(min(a_pos, a_neg, b_pos, b_neg) < MIN_N),
    }

def compare_continuous(x: pd.Series, y: pd.Series) -> dict:
    """Two distributions: per-arm n/median/mean/SD + MW p + Cliff's delta."""
    x = pd.to_numeric(x, errors="coerce").dropna()
    y = pd.to_numeric(y, errors="coerce").dropna()
    out = {
        "n_group": len(x), "n_ref": len(y),
        "median_group": x.median() if len(x) else np.nan,
        "median_ref": y.median() if len(y) else np.nan,
        "mean_group": x.mean() if len(x) else np.nan,
        "mean_ref": y.mean() if len(y) else np.nan,
        "sd_group": x.std() if len(x) else np.nan,
        "sd_ref": y.std() if len(y) else np.nan,
        "descriptive_only": bool(min(len(x), len(y)) < MIN_N),
        "p_value": np.nan, "cliffs_delta": np.nan,
    }
    if len(x) >= 2 and len(y) >= 2:
        try:
            u = mannwhitneyu(x, y, alternative="two-sided")
            out["p_value"] = float(u.pvalue)  # type: ignore
            out["cliffs_delta"] = float(2 * u.statistic / (len(x) * len(y)) - 1)  # type: ignore
        except ValueError:
            pass
    return out

# ---------------- loaders ----------------

def load_annotations(path: str) -> pd.DataFrame:
    df = pd.read_excel(path)
    df.columns = [normalized_name(c) for c in df.columns]
    out = pd.DataFrame({
        "uniprot_accession": df["proteinid"].map(extract_accession),
        "method": df["sourcemethod"].map(canonical_method),
        "start": pd.to_numeric(df["start"], errors="coerce").astype("Int64"),
        "end": pd.to_numeric(df["end"], errors="coerce").astype("Int64"),
        "length": pd.to_numeric(df["length"], errors="coerce"),
        "signature": df["primaryphysicochemicalannotation"].fillna("missing").astype(str),
    })
    return out.dropna(subset=["start", "end", "length"]).copy()

def load_features(path: str) -> pd.DataFrame:
    """Merge feature sheets on interval keys; metric columns renamed back to
    canonical form (frac_A, frac_polar, ...). Method stays in the keys --
    dropping it caused cross-method coordinate collisions in v1."""
    sheets = pd.read_excel(path, sheet_name=None)
    frames = []
    for name, raw in sheets.items():
        raw.columns = [normalized_name(c) for c in raw.columns]
        base = pd.DataFrame({
            "uniprot_accession": raw["proteinid"].map(extract_accession),
            "method": raw["sourcemethod"].map(canonical_method),
            "start": pd.to_numeric(raw["start"], errors="coerce").astype("Int64"),
            "end": pd.to_numeric(raw["end"], errors="coerce").astype("Int64"),
        })
        metrics = raw.drop(columns=[c for c in RAW_KEY_COLUMNS if c in raw.columns])
        metrics.columns = [canonical_metric_name(c) for c in metrics.columns]
        metrics = metrics.loc[:, ~metrics.columns.duplicated()]
        frames.append(pd.concat([base, metrics], axis=1))
    merged = frames[0]
    for f in frames[1:]:
        merged = merged.merge(f, on=["uniprot_accession", "method", "start", "end"],
                              how="outer")
    return merged

def load_universe(background: str, labels: str, patch: str, classes: str):
    """Protein universe: length, RBP label, and (for RBPs) RNA class per
    background protein. Class accessions map through the label patch."""
    bg = pd.read_excel(background)
    bg.columns = [normalized_name(c) for c in bg.columns]
    uni = bg[["uniprotaccession", "uniprotlength"]].copy()
    uni.columns = ["uniprot_accession", "uniprot_length"]
    uni["uniprot_accession"] = uni["uniprot_accession"].astype(str).str.strip()
    uni["uniprot_length"] = pd.to_numeric(uni["uniprot_length"], errors="coerce")
    uni = uni.drop_duplicates("uniprot_accession")

    lab = pd.read_csv(labels, sep="\t")
    lab["is_rbp"] = lab["is_rbp"].astype(str).str.lower().isin({"true", "1", "yes"})
    uni = uni.merge(lab[["uniprot_accession", "is_rbp", "label_source"]],
                    on="uniprot_accession", how="left")
    uni["is_rbp"] = uni["is_rbp"].fillna(False)
    uni["label_source"] = uni["label_source"].fillna("")

    cls = pd.concat(pd.read_excel(classes, sheet_name=None).values(),
                    ignore_index=True)
    cls.columns = [normalized_name(c) for c in cls.columns]
    cls = cls[["uniprotaccession", "rnaprimaryclass"]].dropna(subset=["uniprotaccession"])
    cls["uniprotaccession"] = cls["uniprotaccession"].astype(str).str.strip()
    conflicts = cls.groupby("uniprotaccession")["rnaprimaryclass"].nunique().gt(1).sum()
    if conflicts:
        print(f"WARNING: {conflicts} accessions with conflicting classes; keeping first.")
    cls = cls.drop_duplicates("uniprotaccession")

    patch_df = pd.read_csv(patch, sep="\t")
    acc_map = dict(zip(patch_df["uniprot_accession"].astype(str),
                       patch_df["covering_accession"].astype(str)))
    universe_accs = set(uni["uniprot_accession"])

    def map_acc(a: str):
        if a in universe_accs:
            return a
        cov = acc_map.get(a)
        return cov if cov in universe_accs else None

    cls["universe_acc"] = cls["uniprotaccession"].map(map_acc)
    n_mapped = cls["universe_acc"].notna().sum()
    print(f"classes: {len(cls)} classified, {n_mapped} mapped to background "
          f"({len(cls) - n_mapped} unmapped -> excluded from Mode B)")
    cls = cls.dropna(subset=["universe_acc"])
    uni["rna_class"] = uni["uniprot_accession"].map(
        dict(zip(cls["universe_acc"], cls["rnaprimaryclass"])))
    return uni

# ---------------- analysis blocks ----------------

def protein_table(lcr: pd.DataFrame, uni: pd.DataFrame) -> pd.DataFrame:
    """Per-protein architecture; non-carriers included with n_lcr = 0."""
    agg = (lcr.groupby("uniprot_accession")
           .agg(n_lcr=("length", "size"), lcr_residues=("length", "sum"),
                mean_lcr_length=("length", "mean"))
           .reset_index())
    p = uni.merge(agg, on="uniprot_accession", how="left")
    p[["n_lcr", "lcr_residues"]] = p[["n_lcr", "lcr_residues"]].fillna(0)
    p["coverage"] = p["lcr_residues"] / p["uniprot_length"]
    return p

def block_incidence(p: pd.DataFrame, group_col: str, ref_mask: pd.Series,
                    group_values: list, block: str) -> pd.DataFrame:
    """Carrier rate per group vs. reference, overall and per length quintile."""
    rows = []
    carriers = p["n_lcr"] > 0
    ref_pos, ref_tot = int(carriers[ref_mask].sum()), int(ref_mask.sum())
    for val in group_values:
        mask = p[group_col] == val
        r = compare_binary(int(carriers[mask].sum()), int(mask.sum()),
                           ref_pos, ref_tot)
        r.update({"block": block, "group": val, "length_bin": ""})
        rows.append(r)
    try:  # length-stratified: quintiles of the reference length distribution
        quint = pd.qcut(p.loc[ref_mask, "uniprot_length"], 5, duplicates="drop")
        bins = pd.IntervalIndex(quint.cat.categories)  # type: ignore
        p_len_bin = pd.cut(p["uniprot_length"], bins)
        for b in bins:
            in_bin = p_len_bin == b
            ref_b = ref_mask & in_bin
            rp, rt = int(carriers[ref_b].sum()), int(ref_b.sum())
            for val in group_values:
                mask = (p[group_col] == val) & in_bin
                if mask.sum() == 0:
                    continue
                r = compare_binary(int(carriers[mask].sum()), int(mask.sum()),
                                   rp, rt)
                r.update({"block": f"{block}_by_length", "group": val,
                          "length_bin": str(b)})
                rows.append(r)
    except ValueError:
        pass
    return pd.DataFrame(rows)

def block_architecture(p: pd.DataFrame, group_col: str, ref_mask: pd.Series,
                       group_values: list, block: str) -> pd.DataFrame:
    """Coverage / n_LCR / mean LCR length among carriers, group vs. reference."""
    rows = []
    pc = p[p["n_lcr"] > 0]
    ref_c = pc[ref_mask.loc[pc.index]]
    for val in group_values:
        g = pc[pc[group_col] == val]
        for metric in ARCHITECTURE_METRICS:
            r = compare_continuous(g[metric], ref_c[metric])
            r.update({"block": block, "group": val, "metric": metric})
            rows.append(r)
    return pd.DataFrame(rows)

def block_composition(lcr: pd.DataFrame, group_col: str, group_values: list,
                      block: str) -> pd.DataFrame:
    """LCR-level composition: group LCRs vs. non-RBP background LCRs."""
    rows = []
    ref = lcr[~lcr["is_rbp"]]
    metrics = [m for m in COMPOSITION_METRICS + ["length"] if m in lcr.columns]
    missing = [m for m in COMPOSITION_METRICS if m not in lcr.columns]
    if missing:
        print(f"  WARNING: {len(missing)} composition metrics missing from "
              f"features: {missing[:5]}...")
    for val in group_values:
        g = lcr[lcr[group_col] == val]
        for metric in metrics:
            r = compare_continuous(g[metric], ref[metric])
            r.update({"block": block, "group": val, "metric": metric})
            rows.append(r)
    return pd.DataFrame(rows)

def block_signatures(lcr: pd.DataFrame, p: pd.DataFrame, group_col: str,
                     group_values: list, block: str) -> pd.DataFrame:
    """Per signature: carrier rate (Fisher) + occupancy among carriers (MW).
    One frame with a `measure` column; the caller splits it into sheets."""
    rows = []
    ref_mask_uni = ~p["is_rbp"]
    for sig in sorted(lcr["signature"].dropna().unique()):
        carriers = set(lcr.loc[lcr["signature"] == sig, "uniprot_accession"])  # type: ignore
        has = p["uniprot_accession"].isin(carriers)
        ref_pos = int(has[ref_mask_uni].sum())
        ref_tot = int(ref_mask_uni.sum())
        for val in group_values:
            mask = p[group_col] == val
            if mask.sum() == 0:
                continue
            r = compare_binary(int(has[mask].sum()), int(mask.sum()),
                               ref_pos, ref_tot)
            r.update({"block": block, "group": val, "signature": sig,
                      "measure": "carrier_rate"})
            rows.append(r)
        occ = (lcr[lcr["signature"] == sig].groupby("uniprot_accession")["length"]
               .sum().rename("sig_residues").reset_index())
        po = p.merge(occ, on="uniprot_accession", how="inner")
        po["occupancy"] = po["sig_residues"] / po["uniprot_length"]
        ref_occ = po[~po["is_rbp"]]
        for val in group_values:
            g = po[po[group_col] == val]
            if len(g) == 0:
                continue  # zero carriers in this arm: the rates sheet says it
            r = compare_continuous(g["occupancy"], ref_occ["occupancy"])
            r.update({"block": block, "group": val, "signature": sig,
                      "measure": "occupancy_among_carriers",
                      "metric": "occupancy"})
            rows.append(r)
    return pd.DataFrame(rows)

# ---------------- validation grid ----------------

def load_internal(internal_dir: str) -> pd.DataFrame:
    """Phase 2.4.1 v2 internal enrichments: per-method workbooks, sheet
    categorical_protein, signature_presence rows (class vs. other RBPs)."""
    frames = []
    for f in sorted(Path(internal_dir).glob("*/tables/*_results.xlsx")):
        method = canonical_method(f.stem.replace("_results", ""))
        df = pd.read_excel(f, sheet_name="categorical_protein")
        df = df[df["analysis"] == "signature_presence"]
        df["method"] = method
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"no 2.4.1 workbooks under {internal_dir}")
    internal = pd.concat(frames, ignore_index=True)
    return internal[["method", "rna_class", "feature_value", "log2_or",
                     "p_fdr"]].rename(
        columns={"feature_value": "signature", "log2_or": "internal_log2_or",
                 "p_fdr": "internal_fdr"})

def build_validation_grid(b3_all: pd.DataFrame, internal_dir: str) -> pd.DataFrame:
    """Internal vs. external log2 OR per signature x class x method + verdict.
    Same metric both sides; only the reference differs."""
    internal = load_internal(internal_dir)
    ext = b3_all[b3_all["measure"] == "carrier_rate"][
        ["method", "group", "signature", "log2_or", "p_fdr",
         "n_pos_group", "n_tot_group"]].copy()
    ext.columns = ["method", "rna_class", "signature", "external_log2_or",
                   "external_fdr", "n_carriers", "n_class"]
    grid = internal.merge(ext, on=["method", "rna_class", "signature"],
                          how="inner")

    def verdict(r) -> str:
        if pd.isna(r["external_fdr"]) or r["n_carriers"] < MIN_N:
            return "external_underpowered"
        int_sig = pd.notna(r["internal_fdr"]) and r["internal_fdr"] < 0.05
        int_dir = np.sign(r["internal_log2_or"]) if pd.notna(r["internal_log2_or"]) else 0
        ext_sig = r["external_fdr"] < 0.05
        ext_dir = np.sign(r["external_log2_or"]) if pd.notna(r["external_log2_or"]) else 0
        int_enriched = int_sig and int_dir > 0
        if int_enriched and ext_sig and ext_dir > 0:
            return "externally_validated"
        if int_enriched and not ext_sig:
            return "rbp_generic"
        if int_enriched and ext_sig and ext_dir < 0:
            return "internal_only_artifact"
        if not int_enriched and ext_sig and ext_dir > 0:
            return "new_external_signal"
        if int_sig and int_dir < 0 and ext_sig and ext_dir < 0:
            return "validated_depletion"
        return "concordant_null"

    grid["verdict"] = grid.apply(verdict, axis=1)
    return grid

# ---------------- main ----------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotations",
                    default=r"lcr_analyses\background_analysis\ensembl_lcr_annotations.xlsx")
    ap.add_argument("--features",
                    default=r"lcr_analyses\background_analysis\ensembl_lcr_features.xlsx")
    ap.add_argument("--background",
                    default=r"datasets\rbps_census\ensembl_v116\finalized_background.xlsx")
    ap.add_argument("--labels",
                    default=r"datasets\rbps_census\ensembl_v116\finalized_background_labels.tsv")
    ap.add_argument("--classes",
                    default=r"rbp_superclasses\rbp_rna_classification.xlsx")
    ap.add_argument("--patch",
                    default=r"datasets\rbps_census\ensembl_v116\finalized_label_patch.tsv")
    ap.add_argument("--internal-dir", default=r"lcr_analyses\quantitative_analysis\results_v2\methods")
    ap.add_argument("--outdir", default=r"lcr_analyses\background_analysis\results_v2")
    args = ap.parse_args()

    out = Path(args.outdir)
    summary_tables = out / "summary" / "tables"
    summary_tables.mkdir(parents=True, exist_ok=True)

    print("Loading universe...")
    uni = load_universe(args.background, args.labels, args.patch, args.classes)
    print(f"universe: {len(uni)} proteins, {int(uni['is_rbp'].sum())} RBP, "
          f"{uni['rna_class'].notna().sum()} class-labelled")

    print("Loading annotations + features...")
    ann = load_annotations(args.annotations)
    feat = load_features(args.features)
    lcr = ann.merge(feat, on=["uniprot_accession", "method", "start", "end"],
                    how="left")
    assert len(lcr) == len(ann), f"merge inflated rows: {len(lcr)} vs {len(ann)}"
    lcr = lcr.merge(uni[["uniprot_accession", "is_rbp", "rna_class",
                         "uniprot_length"]], on="uniprot_accession", how="left")
    print(f"LCRs: {len(lcr)}; length match: {lcr['uniprot_length'].notna().mean():.1%}")

    all_b3, headline = [], []
    for method, d in lcr.groupby("method", sort=True):
        tables = out / "methods" / str(method) / "tables"
        tables.mkdir(parents=True, exist_ok=True)
        p = protein_table(d, uni)
        rbp_mask = p["is_rbp"]
        headline.append({"method": method, "n_lcr": len(d),
                         "n_rbp_carriers": int((p["n_lcr"] > 0)[rbp_mask].sum()),
                         "n_rbp_total": int(rbp_mask.sum()),
                         "n_bg_carriers": int((p["n_lcr"] > 0)[~rbp_mask].sum()),
                         "n_bg_total": int((~rbp_mask).sum())})

        ref_mask = ~p["is_rbp"]
        classes = sorted(p["rna_class"].dropna().unique())

        a_sig = block_signatures(d, p, "is_rbp", [True], "A3_signature")
        b_sig = block_signatures(d, p, "rna_class", classes, "B3_signature")

        # sheets are homogeneous: binary and continuous columns never mix;
        # all-empty columns dropped before ordering
        sheets = {
            "A_incidence": block_incidence(p, "is_rbp", ref_mask, [True],
                                           "A1_incidence"),
            "A_architecture": block_architecture(p, "is_rbp", ref_mask, [True],
                                                 "A1_architecture"),
            "A_composition": block_composition(d, "is_rbp", [True],
                                               "A2_composition"),
            "A_signature_rates": a_sig[a_sig["measure"] == "carrier_rate"],
            "A_signature_occupancy": a_sig[a_sig["measure"] != "carrier_rate"],
            "B_incidence": block_incidence(p, "rna_class", ref_mask, classes,
                                           "B1_incidence"),
            "B_composition": block_composition(
                d[d["rna_class"].notna() | ~d["is_rbp"]], "rna_class", classes,
                "B2_composition"),
            "B_signature_rates": b_sig[b_sig["measure"] == "carrier_rate"],
            "B_signature_occupancy": b_sig[b_sig["measure"] != "carrier_rate"],
        }
        workbook = tables / f"{method}_results.xlsx"
        with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
            for name, df in sheets.items():
                df = df.dropna(axis=1, how="all")
                if len(df):
                    order_cols(bh_adjust_per_block(df)).to_excel(
                        writer, sheet_name=name, index=False)
        b3 = bh_adjust_per_block(b_sig)
        b3["method"] = method
        all_b3.append(b3)
        print(f"done {method}: {len(d)} LCRs -> {workbook}")

    pd.DataFrame(headline).to_csv(summary_tables / "headline_numbers.csv",
                                  index=False)
    print("Building validation grid...")
    grid = build_validation_grid(pd.concat(all_b3, ignore_index=True),
                                 args.internal_dir)
    grid.to_csv(summary_tables / "validation_grid.csv", index=False)
    consensus = (grid.groupby(["rna_class", "signature"])["verdict"]
                 .agg(lambda s: s.value_counts().to_dict()).reset_index())
    consensus.to_csv(summary_tables / "validation_consensus.csv", index=False)
    print(grid["verdict"].value_counts().to_string())
    print(f"\nAll outputs under: {out.resolve()}")

if __name__ == "__main__":
    main()
