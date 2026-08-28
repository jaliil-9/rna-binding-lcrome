"""
LCR clustering — Analysis B: carrier phenotype (LCR carriers only, per method).

For each method:
  - Restrict to proteins with n_lcr > 0.
  - Build protein-level features:
      * architecture: n_lcr, coverage, mean_lcr_length
      * position presence
      * signature presence
      * LCR-residue-weighted physicochemical fractions
  - Run PAM, hierarchical, HDBSCAN.
  - Save cluster labels and summaries.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import RobustScaler
import hdbscan
import kmedoids
from scipy.spatial.distance import squareform


# ---------------- configuration ----------------

METHOD_ALIASES = {
    "cast": "CAST",
    "seg": "SEG",
    "segintermediate": "SEG_intermediate",
    "flps": "FLPS",
    "lcrfinder": "LCRFinder",
    "alcor": "AlcoR",
}

POSITION_CLASSES = [
    "domainintrinsic",
    "domainedge",
    "domainadjacent",
    "interdomainlinker",
    "distalterminal",
]

SIGNATURES = [
    "rg_rgg_repeat_region",
    "sr_rs_repeat_region",
    "gs_rich_neutral_region",
    "basic_enriched_region",
    "arginine_rich_rna_contact_region",
    "acidic_region",
    "mixed_charge_polyampholyte",
    "polar_linker",
    "hydrophobic_region",
    "aromatic_sticker_region",
    "aromatic_aggregation_prone_region",
    "cation_pi_rich_neighborhood",
    "disorder_rich_spacer_region",
    "serine_rich_region",
    "glutamine_rich_region",
    "glycine_rich_region",
    "proline_rich_disordered_region",
    "unmapped_physicochemical_properties",
]

PROPERTY_METRICS = [
    "frac_polar",
    "frac_hydrophobic",
    "frac_strong_hydro",
    "frac_aromatic",
    "frac_disorder",
    "frac_positive",
    "frac_negative",
    "fcr",
    "ncpr",
    "frac_GS",
]

K_VALUES = [2, 3, 4, 5, 6, 7, 8, 9]

MIN_CLUSTER_SIZE_ABS = 10
MIN_CLUSTER_SIZE_FRAC = 0.03


# ---------------- helpers ----------------

def normalized_name(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).strip().lower())

def canonical_method(x: str) -> str:
    return METHOD_ALIASES.get(normalized_name(x), str(x).strip())

def extract_accession(x: str) -> str:
    text = str(x).strip()
    if "|" in text:
        return text.split("|", 1)[0].strip()
    return text

def load_annotations(path: str) -> pd.DataFrame:
    raw = pd.read_excel(path)
    raw.columns = [normalized_name(c) for c in raw.columns]
    out = pd.DataFrame({
        "uniprot_accession": raw["proteinid"].map(extract_accession),
        "method": raw["sourcemethod"].map(canonical_method),
        "start": pd.to_numeric(raw["start"], errors="coerce").astype("Int64"),
        "end": pd.to_numeric(raw["end"], errors="coerce").astype("Int64"),
        "length": pd.to_numeric(raw["length"], errors="coerce"),
        "rna_class": raw["rnatargetsuperclass"].fillna("missing").astype(str),
        "position": raw["domainpositionclass"].fillna("unclassifiednopfam").astype(str),
        "signature": raw["primaryphysicochemicalannotation"].fillna("unmapped_physicochemical_properties").astype(str),
    })
    return out.dropna(subset=["start", "end", "length"]).copy()

def load_features(path: str) -> pd.DataFrame:
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
        metrics = raw.drop(columns=[c for c in {"proteinid", "sourcemethod", "start", "end", "length"} if c in raw.columns])
        # keep only known property metrics
        metrics = metrics[[c for c in metrics.columns if c in PROPERTY_METRICS]]
        frames.append(pd.concat([base, metrics], axis=1))
    merged = frames[0]
    for f in frames[1:]:
        merged = merged.merge(f, on=["uniprot_accession", "method", "start", "end"], how="outer")
    return merged

def load_master(path: str) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name="Combined")
    raw.columns = [normalized_name(c) for c in raw.columns]
    acc_col = len_col = None
    for c in raw.columns:
        if c in {"uniprot_accession", "uniprotaccession"}:
            acc_col = c
        if c in {"uniprot_length", "uniprotlength"}:
            len_col = c
    if acc_col is None or len_col is None:
        raise ValueError("Master file must contain accession and length columns")
    out = raw[[acc_col, len_col]].copy()
    out.columns = ["uniprot_accession", "uniprot_length"]
    out["uniprot_accession"] = out["uniprot_accession"].astype(str).str.strip()
    out["uniprot_length"] = pd.to_numeric(out["uniprot_length"], errors="coerce")
    return out.drop_duplicates("uniprot_accession")

def load_rna_classes(path: str) -> pd.DataFrame:
    df = pd.read_excel(path)
    df.columns = [normalized_name(c) for c in df.columns]
    out = df[["uniprotaccession", "rnaprimaryclass"]].dropna(subset=["uniprotaccession"]).copy()
    out["uniprot_accession"] = out["uniprotaccession"].astype(str).str.strip()
    out = out.drop_duplicates("uniprot_accession")
    return out[["uniprot_accession", "rnaprimaryclass"]]

def build_carrier_table(
    ann: pd.DataFrame,
    feat: pd.DataFrame,
    master: pd.DataFrame,
) -> pd.DataFrame:
    """Aggregate LCR-level data to protein level for carriers only."""
    # architecture
    arch = (
        ann.groupby("uniprot_accession")
        .agg(
            n_lcr=("length", "size"),
            lcr_residues=("length", "sum"),
            mean_lcr_length=("length", "mean"),
        )
        .reset_index()
    )
    p = master.merge(arch, on="uniprot_accession", how="inner")
    p = p[p["n_lcr"] > 0].copy()
    p["coverage"] = p["lcr_residues"] / p["uniprot_length"]

    # position presence
    ann_pos = ann[ann["position"] != "unclassifiednopfam"].copy()
    pos_present = (
        ann_pos.groupby(["uniprot_accession", "position"])
        .size()
        .reset_index(name="count")
    )
    pos_present["present"] = 1
    pos_piv = pos_present.pivot_table(
        index="uniprot_accession",
        columns="position",
        values="present",
        fill_value=0,
    )
    for col in POSITION_CLASSES:
        if col not in pos_piv.columns:
            pos_piv[col] = 0
    pos_piv = pos_piv[POSITION_CLASSES]
    p = p.merge(pos_piv, on="uniprot_accession", how="left")
    p[POSITION_CLASSES] = p[POSITION_CLASSES].fillna(0).astype(int)

    # signature presence
    sig_present = (
        ann.groupby(["uniprot_accession", "signature"])
        .size()
        .reset_index(name="count")
    )
    sig_present["present"] = 1
    sig_piv = sig_present.pivot_table(
        index="uniprot_accession",
        columns="signature",
        values="present",
        fill_value=0,
    )
    for sig in SIGNATURES:
        if sig not in sig_piv.columns:
            sig_piv[sig] = 0
    sig_piv = sig_piv[SIGNATURES]
    p = p.merge(sig_piv, on="uniprot_accession", how="left")
    p[SIGNATURES] = p[SIGNATURES].fillna(0).astype(int)

    # LCR-residue-weighted physicochemical fractions
    merged = ann.merge(
        feat,
        on=["uniprot_accession", "method", "start", "end"],
        how="left",
        validate="one_to_many",
    )
    merged = merged[merged["uniprot_accession"].isin(p["uniprot_accession"])].copy()
    for metric in PROPERTY_METRICS:
        if metric not in merged.columns:
            merged[metric] = np.nan
    # weight by LCR length
    merged["weight"] = merged["length"]
    agg = (
        merged.groupby("uniprot_accession")[PROPERTY_METRICS + ["weight"]]
        .apply(lambda g: (g[PROPERTY_METRICS].multiply(g["weight"], axis=0)).sum() / g["weight"].sum())
        .reset_index()
    )
    p = p.merge(agg, on="uniprot_accession", how="left")
    for m in PROPERTY_METRICS:
        p[m] = pd.to_numeric(p[m], errors="coerce")

    return p

def gower_mixed_distance(
    df: pd.DataFrame,
    binary_cols: List[str],
    continuous_cols: List[str],
) -> np.ndarray:
    n = df.shape[0]
    dist = np.zeros((n, n), dtype=float)
    w_bin = 0.5 / max(len(binary_cols), 1)
    w_con = 0.5 / max(len(continuous_cols), 1)

    if binary_cols:
        Xb = df[binary_cols].to_numpy()
        db = pairwise_distances(Xb, metric="hamming")
        dist += w_bin * db

    if continuous_cols:
        X = df[continuous_cols].to_numpy()
        scaler = RobustScaler()
        Xs = scaler.fit_transform(X)
        dc = pairwise_distances(Xs, metric="manhattan")
        dist += w_con * dc

    return dist

def run_pam(distance: np.ndarray, k: int, random_state: int = 0) -> Tuple[np.ndarray, object]:
    km = kmedoids.KMedoids(n_clusters=k, metric="precomputed", random_state=random_state, init="build")
    labels = km.fit_predict(distance)
    return labels, km

def run_hierarchical(
    distance: np.ndarray,
    k: int,
) -> np.ndarray:
    from scipy.cluster.hierarchy import linkage, fcluster

    condensed_distance = squareform(distance, checks=False)
    Z = linkage(condensed_distance, method="average")

    return fcluster(Z, t=k, criterion="maxclust")

def run_hdbscan(
    X: np.ndarray,
    min_cluster_size: int,
    min_samples: int | None = None,
    metric: str = "euclidean",
) -> Tuple[np.ndarray, np.ndarray]:
    if min_samples is None:
        min_samples = min_cluster_size
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        metric=metric,
        cluster_selection_method="eom",
    )
    clusterer.fit(X)
    labels = clusterer.labels_
    membership = clusterer.probabilities_ if hasattr(clusterer, "probabilities_") else None
    return labels, membership # type: ignore

def summarize_clusters(
    df: pd.DataFrame,
    labels: np.ndarray,
    rna_col: str | None = None,
) -> pd.DataFrame:
    out = df.copy()
    out["cluster"] = labels
    summ = []
    for c in sorted(set(labels)):
        if c == -1:
            continue
        sub = out[out["cluster"] == c]
        row = {
            "cluster": int(c),
            "size": int(len(sub)),
        }
        # architecture
        for col in ["n_lcr", "coverage", "mean_lcr_length"]:
            if col in sub.columns:
                row[f"median_{col}"] = float(sub[col].median()) # type: ignore
        # position
        for pos in POSITION_CLASSES:
            if pos in sub.columns:
                row[f"prop_{pos}"] = float(sub[pos].mean()) # type: ignore
        # signature
        for sig in SIGNATURES:
            if sig in sub.columns:
                row[f"prop_sig_{sig}"] = float(sub[sig].mean()) # type: ignore
        # composition
        for m in PROPERTY_METRICS:
            if m in sub.columns:
                row[f"median_{m}"] = float(sub[m].median()) # type: ignore
        # RNA class
        if rna_col and rna_col in sub.columns:
            class_counts = sub[rna_col].value_counts().to_dict()
            for cls, cnt in class_counts.items():
                row[f"class_{cls}"] = int(cnt)
        summ.append(row)
    return pd.DataFrame(summ)

# ---------------- main ----------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annotations", default="lcr_analyses/pc_properties/lcr_annotations.xlsx")
    ap.add_argument("--features", default="lcr_analyses/pc_properties/lcr_features.xlsx")
    ap.add_argument("--master", default="datasets/combined_rbp_pfam38_rbpdb_modomics_uniprot.xlsx")
    ap.add_argument("--rna-classes", default="rbp_superclasses/rbp_rna_classification.xlsx")
    ap.add_argument("--outdir", default="lcr_analyses/clustering/carriers")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    ann = load_annotations(args.annotations)
    feat = load_features(args.features)
    master = load_master(args.master)
    rna = load_rna_classes(args.rna_classes)

    rbps = master[master["uniprot_accession"].isin(rna["uniprot_accession"])].copy()
    print(f"RBP universe: {len(rbps)} proteins")

    ann = ann[ann["uniprot_accession"].isin(rbps["uniprot_accession"])].copy()
    feat = feat[feat["uniprot_accession"].isin(rbps["uniprot_accession"])].copy()

    results = []

    for method, d in ann.groupby("method", sort=True):
        print(f"\n=== Method: {method} ===")
        method_dir = out / str(method)
        method_dir.mkdir(parents=True, exist_ok=True)

        d_feat = feat[feat["method"] == method].copy()
        p = build_carrier_table(d, d_feat, rbps)
        p = p.merge(rna, on="uniprot_accession", how="left")

        # feature matrix
        binary_cols = POSITION_CLASSES + SIGNATURES
        continuous_cols = ["n_lcr", "coverage", "mean_lcr_length"] + PROPERTY_METRICS

        X = p[binary_cols + continuous_cols].copy().astype(float)

        # distance
        print("Computing Gower-style distance...")
        dist = gower_mixed_distance(X, binary_cols=binary_cols, continuous_cols=continuous_cols)

        # save results
        out_df = p[
            ["uniprot_accession", "uniprot_length", "rnaprimaryclass"]
        ].copy()
        out_df["method"] = method

        pam_summaries = {}
        hierarchical_summaries = {}

        for k in K_VALUES:
            print(f"Running PAM (k={k})...")
            labels_pam, _ = run_pam(
                dist,
                k=k,
                random_state=0,
            )
            out_df[f"cluster_pam_k{k}"] = labels_pam

            pam_summaries[f"k{k}"] = summarize_clusters(
                p,
                labels_pam,
                rna_col="rnaprimaryclass",
            )

            print(f"Running hierarchical clustering (k={k})...")
            labels_hc = run_hierarchical(
                dist,
                k=k,
            )
            out_df[f"cluster_hierarchical_k{k}"] = labels_hc

            hierarchical_summaries[f"k{k}"] = summarize_clusters(
                p,
                labels_hc,
                rna_col="rnaprimaryclass",
            )

            results.append({
                "method": method,
                "analysis": "global",  # use "carriers" in carrier script
                "algorithm": "pam",
                "k": k,
                "n_proteins": len(p),
            })
            results.append({
                "method": method,
                "analysis": "global",  # use "carriers" in carrier script
                "algorithm": "hierarchical",
                "k": k,
                "n_proteins": len(p),
            })

        with pd.ExcelWriter(
            method_dir / "cluster_summaries_pam.xlsx",
            engine="openpyxl",
        ) as writer:
            for sheet_name, summary_df in pam_summaries.items():
                summary_df.to_excel(
                    writer,
                    sheet_name=sheet_name,
                    index=False,
                )

        with pd.ExcelWriter(
            method_dir / "cluster_summaries_hierarchical.xlsx",
            engine="openpyxl",
        ) as writer:
            for sheet_name, summary_df in hierarchical_summaries.items():
                summary_df.to_excel(
                    writer,
                    sheet_name=sheet_name,
                    index=False,
                )

        out_df.to_csv(method_dir / "protein_clusters.csv", index=False)

        # HDBSCAN
        print("Running HDBSCAN...")
        scaler = RobustScaler()
        Xs = scaler.fit_transform(X)
        min_cs = max(MIN_CLUSTER_SIZE_ABS, int(MIN_CLUSTER_SIZE_FRAC * len(Xs)))
        min_s = max(5, min_cs // 2)
        labels_hdbscan, membership = run_hdbscan(
            Xs,
            min_cluster_size=min_cs,
            min_samples=min_s,
            metric="euclidean",
        )

        # save
        out_df["cluster_hdbscan"] = labels_hdbscan
        if membership is not None:
            out_df["hdbscan_confidence"] = membership
        else:
            out_df["hdbscan_confidence"] = np.nan

        out_df.to_csv(method_dir / "protein_clusters.csv", index=False)


    pd.DataFrame(results).to_csv(out / "method_manifest.csv", index=False)
    print(f"\nDone. Outputs under: {out.resolve()}")


if __name__ == "__main__":
    main()