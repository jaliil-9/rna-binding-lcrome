"""
LCR clustering — Analysis A: global LCR phenotype (all RBPs, per method).

For each LCR detection method independently:
  - Build a protein-level feature matrix including non-carriers.
  - Run three clusterers: PAM/k-medoids, hierarchical, HDBSCAN.
  - Save cluster labels, medoids, and per-cluster summaries.

No cross-method consensus, no RNA-class usage in clustering.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator
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

# Domain-position classes observed in annotations
POSITION_CLASSES = [
    "domain_intrinsic",
    "domain_edge",
    "domain_adjacent",
    "interdomain_linker",
    "distal_terminal",
]

# Signature catalogue (primary_physicochemical_annotation values)
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

K_VALUES = [2, 3, 4, 5, 6]

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
        "position": raw["domainpositionclass"].fillna("position_no_pfam").astype(str),
        "signature": raw["primaryphysicochemicalannotation"].fillna("unmapped_physicochemical_properties").astype(str),
    })
    return out.dropna(subset=["start", "end", "length"]).copy()

def load_master(path: str) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name="Combined")
    raw.columns = [normalized_name(c) for c in raw.columns]
    acc_col = None
    len_col = None
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

def build_protein_table(ann: pd.DataFrame, master: pd.DataFrame) -> pd.DataFrame:
    """Per-protein LCR architecture; non-carriers included with n_lcr = 0."""
    agg = (
        ann.groupby("uniprot_accession")
        .agg(
            n_lcr=("length", "size"),
            lcr_residues=("length", "sum"),
        )
        .reset_index()
    )
    p = master.merge(agg, on="uniprot_accession", how="left")
    p["n_lcr"] = p["n_lcr"].fillna(0).astype(int)
    p["lcr_residues"] = p["lcr_residues"].fillna(0.0)
    p["coverage"] = p["lcr_residues"] / p["uniprot_length"]
    return p

def encode_position(ann: pd.DataFrame, proteins: pd.DataFrame) -> pd.DataFrame:
    """Protein-level binary presence of each domain-position class."""
    ann = ann[ann["position"] != "position_no_pfam"].copy()
    present = (
        ann.groupby(["uniprot_accession", "position"])
        .size()
        .reset_index(name="count")
    )
    present["present"] = 1
    piv = present.pivot_table(
        index="uniprot_accession",
        columns="position",
        values="present",
        fill_value=0,
    )
    # ensure all expected columns exist
    for col in POSITION_CLASSES:
        if col not in piv.columns:
            piv[col] = 0
    piv = piv[POSITION_CLASSES]
    out = proteins.merge(piv, on="uniprot_accession", how="left")
    out[POSITION_CLASSES] = out[POSITION_CLASSES].fillna(0).astype(int)
    return out

def encode_signature(ann: pd.DataFrame, proteins: pd.DataFrame) -> pd.DataFrame:
    """Protein-level binary presence of each primary signature."""
    present = (
        ann.groupby(["uniprot_accession", "signature"])
        .size()
        .reset_index(name="count")
    )
    present["present"] = 1
    piv = present.pivot_table(
        index="uniprot_accession",
        columns="signature",
        values="present",
        fill_value=0,
    )
    for sig in SIGNATURES:
        if sig not in piv.columns:
            piv[sig] = 0
    piv = piv[SIGNATURES]
    out = proteins.merge(piv, on="uniprot_accession", how="left")
    out[SIGNATURES] = out[SIGNATURES].fillna(0).astype(int)
    return out

def build_multiplicity_category(n_lcr: pd.Series) -> pd.Series:
    """Map n_lcr to '0', '1', '2', '3+' categories."""
    def cat(x: int) -> str:
        if x == 0:
            return "0"
        if x == 1:
            return "1"
        if x == 2:
            return "2"
        return "3+"
    return n_lcr.apply(cat)

def gower_mixed_distance(
    df: pd.DataFrame,
    binary_cols: List[str],
    categorical_cols: List[str],
    continuous_cols: List[str],
) -> np.ndarray:
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

def run_pam(
    distance: np.ndarray,
    k: int,
    random_state: int = 0,
) -> Tuple[np.ndarray, BaseEstimator]:
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
        # architecture summaries
        if "coverage" in sub.columns:
            row["median_coverage"] = float(sub["coverage"].median()) # type: ignore
        if "n_lcr" in sub.columns:
            row["median_n_lcr"] = float(sub["n_lcr"].median()) # type: ignore
        # position signatures
        for pos in POSITION_CLASSES:
            if pos in sub.columns:
                row[f"prop_{pos}"] = float(sub[pos].mean()) # type: ignore
        # signature presence
        for sig in SIGNATURES:
            if sig in sub.columns:
                row[f"prop_sig_{sig}"] = float(sub[sig].mean()) # type: ignore
        # RNA-class composition if available
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
    ap.add_argument("--master", default="datasets/combined_rbp_pfam38_rbpdb_modomics_uniprot.xlsx")
    ap.add_argument("--rna-classes", default="rbp_superclasses/rbp_rna_classification.xlsx")
    ap.add_argument("--outdir", default="lcr_analyses/clustering/global")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    ann = load_annotations(args.annotations)
    master = load_master(args.master)
    rna = load_rna_classes(args.rna_classes)

    # restrict to RBPs only (present in master and with RNA class)
    rbps = master[master["uniprot_accession"].isin(rna["uniprot_accession"])].copy()
    print(f"RBP universe: {len(rbps)} proteins")

    ann = ann[ann["uniprot_accession"].isin(rbps["uniprot_accession"])].copy()

    results = []

    for method, d in ann.groupby("method", sort=True):
        print(f"\n=== Method: {method} ===")
        method_dir = out / str(method)
        method_dir.mkdir(parents=True, exist_ok=True)

        # protein table
        p = build_protein_table(d, rbps)
        p = p.merge(rna, on="uniprot_accession", how="left")

        # encode features
        p = encode_position(d, p)
        p = encode_signature(d, p)
        p["multiplicity_cat"] = build_multiplicity_category(p["n_lcr"])

        # feature matrix for clustering
        binary_cols = POSITION_CLASSES + SIGNATURES
        categorical_cols = ["multiplicity_cat"]
        continuous_cols = ["coverage"]

        X_cat = p[binary_cols + categorical_cols + continuous_cols].copy()

        # distance matrix
        print("Computing Gower-style distance...")
        dist = gower_mixed_distance(
            X_cat,
            binary_cols=binary_cols,
            categorical_cols=categorical_cols,
            continuous_cols=continuous_cols,
        )

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
                "analysis": "global",
                "algorithm": "pam",
                "k": k,
                "n_proteins": len(p),
            })
            results.append({
                "method": method,
                "analysis": "global",
                "algorithm": "hierarchical",
                "k": k,
                "n_proteins": len(p),
            })
        
        out_df.to_csv(method_dir / "protein_clusters.csv", index=False)    
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

        # HDBSCAN on scaled continuous + binary
        print("Running HDBSCAN...")
        X_num = p[binary_cols + continuous_cols].astype(float).to_numpy()
        scaler = RobustScaler()
        Xs = scaler.fit_transform(X_num)
        min_cs = max(MIN_CLUSTER_SIZE_ABS, int(MIN_CLUSTER_SIZE_FRAC * len(Xs)))
        min_s = max(5, min_cs // 2)
        labels_hdbscan, membership = run_hdbscan(
            Xs,
            min_cluster_size=min_cs,
            min_samples=min_s,
            metric="euclidean",
        )

        # save results
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