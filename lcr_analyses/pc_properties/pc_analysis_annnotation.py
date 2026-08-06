#!/usr/bin/env python3
"""
LCR physicochemical feature extraction and behavior annotation (v1).

Inputs:
    - lcr_methods_combined.xlsx  (sheet: "all_combined")
    - aa-physicochemical-properties-2.csv

Outputs:
    - lcr_features.xlsx
        * composition
        * distribution
        * co_occurrence
    - lcr_annotations.xlsx
        * metadata + RNA-target/domain-position (if available)
        * physicochemical behavior annotations and evidence

This script assumes both input files are in the current working directory.
"""

import pandas as pd
import numpy as np
from collections import defaultdict
import argparse

# =========================
# Configuration
# =========================

LCR_FILE = r"rbp_lcrs\lcr_methods_combined.xlsx"
AA_PROP_FILE = r"lcr_analyses\pc_properties\aa-physicochemical-properties.csv"
LCR_SHEET = "all_results"

OUTPUT_FEATURES = "lcr_features.xlsx"
OUTPUT_ANNOTATIONS = "lcr_annotations.xlsx"

# Co-occurrence window size
WINDOW_SIZE = 5

# Behavior-rule thresholds (can be tuned later)
THRESHOLDS = {
    "f_plus_min": 0.25,
    "ncpr_basic_min": 0.10,
    "f_minus_min": 0.25,
    "ncpr_acidic_max": -0.10,
    "fcr_polyampholyte_min": 0.30,
    "ncpr_polyampholyte_abs_max": 0.05,
    "f_polar_min": 0.50,
    "f_hydro_max_for_polar_linker": 0.20,
    "f_hydro_min": 0.40,
    "f_polar_mixed_min": 0.30,
    "f_polar_mixed_max": 0.50,
    "f_hydro_mixed_min": 0.30,
    "f_hydro_mixed_max": 0.50,
    "f_aromatic_min": 0.08,
    "f_small_min": 0.50,
    "f_small_max_for_bulky": 0.20,
    "cooc_top_quantile": 0.75,
}

# =========================
# Load data
# =========================

def load_aa_properties(filepath):
    """Load AA property CSV and return a dict keyed by single-letter code."""
    aa_df = pd.read_csv(filepath)
    aa_df.columns = aa_df.columns.str.strip()
    # Extract single-letter code from "AA" column like "A (Ala)"
    aa_df["AA_letter"] = aa_df["AA"].str.extract(r"^\s*([A-Z])\s*\(")[0]

    aa_lookup = {}
    for _, row in aa_df.iterrows():
        aa = row["AA_letter"]
        if pd.isna(aa):
            continue
        aa_lookup[aa] = {
            "charge": row["Charge (7.4)"],  # Pos, Neg, Neu, Neu*
            "polar": row["Polar?"] == "Yes",
            "hydrophobic": row["Hydrophobic?"] == "Yes",
            "aromatic": row["Aromatic?"] == "Yes",
            "small": row["Small?"] == "Yes",
        }
    return aa_lookup


def load_lcrs(filepath, sheet_name):
    """Load LCR table."""
    df = pd.read_excel(filepath, sheet_name=sheet_name)
    # Ensure key columns exist
    required_cols = ["protein_id", "method", "start", "end", "length", "sequence"]
    for c in required_cols:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")
    return df


# =========================
# Sequence → property tracks
# =========================

def encode_sequence(seq, aa_lookup):
    """
    For a given protein sequence, return:
      - dict of property tracks (lists of 0/1 or charge states)
      - length
    """
    seq = seq.upper()
    L = len(seq)

    # Initialize tracks
    polar = []
    hydro = []
    aromatic = []
    small = []
    charge = []  # +1, -1, 0

    for res in seq:
        if res not in aa_lookup:
            # Unknown residue (e.g. X); treat as neutral, non-property
            polar.append(0)
            hydro.append(0)
            aromatic.append(0)
            small.append(0)
            charge.append(0)
            continue

        props = aa_lookup[res]
        polar.append(1 if props["polar"] else 0)
        hydro.append(1 if props["hydrophobic"] else 0)
        aromatic.append(1 if props["aromatic"] else 0)
        small.append(1 if props["small"] else 0)

        ch = props["charge"]
        if ch == "Pos":
            charge.append(1)
        elif ch == "Neg":
            charge.append(-1)
        else:
            # Neu, Neu*, anything else
            charge.append(0)

    tracks = {
        "polar": np.array(polar, dtype=int),
        "hydrophobic": np.array(hydro, dtype=int),
        "aromatic": np.array(aromatic, dtype=int),
        "small": np.array(small, dtype=int),
        "charged": np.array([1 if c != 0 else 0 for c in charge], dtype=int),
        "charge_signed": np.array(charge, dtype=int),
    }
    return tracks, L


# =========================
# Composition features
# =========================

def compute_composition_features(seq, aa_lookup):
    """Return a dict of composition features for one sequence."""
    seq = seq.upper()
    L = len(seq)
    if L == 0:
        # Return empty features
        return {
            "length": 0,
            **{f"frac_{aa}": 0.0 for aa in "ACDEFGHIKLMNPQRSTVWY"},
            "frac_polar": 0.0,
            "frac_hydrophobic": 0.0,
            "frac_aromatic": 0.0,
            "frac_small": 0.0,
            "frac_positive": 0.0,
            "frac_negative": 0.0,
            "fcr": 0.0,
            "ncpr": 0.0,
        }

    # AA counts
    aa_counts = defaultdict(int)
    for res in seq:
        if res in "ACDEFGHIKLMNPQRSTVWY":
            aa_counts[res] += 1

    features = {"length": L}
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        features[f"frac_{aa}"] = aa_counts[aa] / L # type: ignore

    # Property fractions
    n_polar = 0
    n_hydro = 0
    n_arom = 0
    n_small = 0
    n_pos = 0
    n_neg = 0

    for res in seq:
        if res not in aa_lookup:
            continue
        p = aa_lookup[res]
        if p["polar"]:
            n_polar += 1
        if p["hydrophobic"]:
            n_hydro += 1
        if p["aromatic"]:
            n_arom += 1
        if p["small"]:
            n_small += 1
        ch = p["charge"]
        if ch == "Pos":
            n_pos += 1
        elif ch == "Neg":
            n_neg += 1

    features["frac_polar"] = n_polar / L # type: ignore
    features["frac_hydrophobic"] = n_hydro / L # type: ignore
    features["frac_aromatic"] = n_arom / L # type: ignore
    features["frac_small"] = n_small / L # type: ignore 
    features["frac_positive"] = n_pos / L # type: ignore 
    features["frac_negative"] = n_neg / L   # type: ignore
    features["fcr"] = (n_pos + n_neg) / L # type: ignore
    features["ncpr"] = (n_pos - n_neg) / L # type: ignore

    return features


# =========================
# Distribution features
# =========================

def runs_and_gaps(binary_array):
    """
    Given a 1D binary numpy array, compute:
      - n_runs
      - mean_run_length
      - max_run_length
      - mean_gap
      - cv_gap
      - distribution_label (compact, dispersed, periodic-ish, insufficient)
    """
    if binary_array.sum() == 0:
        return {
            "n_runs": 0,
            "mean_run_length": np.nan,
            "max_run_length": 0,
            "mean_gap": np.nan,
            "cv_gap": np.nan,
            "label": "insufficient",
        }

    # Find runs
    diff = np.diff(np.concatenate(([0], binary_array, [0])))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    run_lengths = ends - starts
    n_runs = len(run_lengths)
    mean_run = float(np.mean(run_lengths))
    max_run = int(np.max(run_lengths))

    # Gaps between runs
    if n_runs <= 1:
        mean_gap = np.nan
        cv_gap = np.nan
        label = "compact" if max_run / len(binary_array) > 0.4 else "dispersed"
    else:
        # gap between run i and i+1: start_{i+1} - end_i
        gaps = starts[1:] - ends[:-1]
        mean_gap = float(np.mean(gaps))
        if mean_gap == 0:
            cv_gap = 0.0
        else:
            cv_gap = float(np.std(gaps) / mean_gap) if mean_gap > 0 else np.nan

        # Simple heuristic for label
        frac_max = max_run / len(binary_array)
        if n_runs <= 3 and frac_max > 0.35:
            label = "compact"
        elif cv_gap < 0.4 and n_runs >= 4:
            label = "periodic-ish"
        else:
            label = "dispersed"

    return {
        "n_runs": n_runs,
        "mean_run_length": mean_run,
        "max_run_length": max_run,
        "mean_gap": mean_gap,
        "cv_gap": cv_gap,
        "label": label,
    }


def compute_distribution_features(seq, aa_lookup):
    """Return distribution features for each property track."""
    tracks, L = encode_sequence(seq, aa_lookup)

    features = {}
    prop_names = ["polar", "hydrophobic", "aromatic", "small", "charged"]

    for prop in prop_names:
        arr = tracks[prop]
        stats = runs_and_gaps(arr)
        for k, v in stats.items():
            features[f"{prop}_{k}"] = v

    return features


# =========================
# Co-occurrence features
# =========================

def compute_cooccurrence_features(seq, aa_lookup, window=WINDOW_SIZE):
    """
    Compute normalized co-occurrence scores for property pairs within a sliding window.
    Returns dict of scores c_{X,Y}.
    """
    tracks, L = encode_sequence(seq, aa_lookup)
    if L < window or L == 0:
        # Not enough data
        pairs = [
            ("positive", "aromatic"),
            ("hydrophobic", "hydrophobic"),
            ("polar", "polar"),
            ("polar", "hydrophobic"),
            ("positive", "negative"),
        ]
        return {f"cooc_{x}_{y}": 0.0 for x, y in pairs}

    # Build property boolean arrays
    pos = tracks["charge_signed"] > 0
    neg = tracks["charge_signed"] < 0
    arom = tracks["aromatic"].astype(bool)
    hydro = tracks["hydrophobic"].astype(bool)
    polar = tracks["polar"].astype(bool)

    # Define pairs
    pair_arrays = {
        ("positive", "aromatic"): (pos, arom),
        ("hydrophobic", "hydrophobic"): (hydro, hydro),
        ("polar", "polar"): (polar, polar),
        ("polar", "hydrophobic"): (polar, hydro),
        ("positive", "negative"): (pos, neg),
    }

    scores = {}
    n_windows = L - window + 1

    for (x, y), (arr_x, arr_y) in pair_arrays.items():
        count = 0
        for i in range(n_windows):
            w_x = arr_x[i : i + window]
            w_y = arr_y[i : i + window]
            if w_x.any() and w_y.any():
                count += 1
        scores[f"cooc_{x}_{y}"] = count / n_windows if n_windows > 0 else 0.0

    return scores


# =========================
# Behavior annotation rules
# =========================

def assign_behavior_labels(comp, distr, cooc, cooc_quantiles=None):
    """
    Given composition, distribution, and co-occurrence features,
    assign behavior labels according to the taxonomy.

    cooc_quantiles: dict mapping cooc feature name -> quantile thresholds (optional).
    Returns:
      - primary_label (str or None)
      - secondary_labels (list of str)
      - evidence (list of strings describing triggered rules)
    """
    labels = []
    evidence = []

    f_pos = comp["frac_positive"]
    f_neg = comp["frac_negative"]
    ncpr = comp["ncpr"]
    fcr = comp["fcr"]

    f_polar = comp["frac_polar"]
    f_hydro = comp["frac_hydrophobic"]
    f_arom = comp["frac_aromatic"]
    f_small = comp["frac_small"]

    charge_label = distr.get("charged_label", "insufficient")
    polar_label = distr.get("polar_label", "insufficient")
    hydro_label = distr.get("hydrophobic_label", "insufficient")
    arom_label = distr.get("aromatic_label", "insufficient")
    small_label = distr.get("small_label", "insufficient")

    # A. Basic-enriched patch
    if (
        f_pos >= THRESHOLDS["f_plus_min"]
        and ncpr >= THRESHOLDS["ncpr_basic_min"]
        and charge_label == "compact"
    ):
        labels.append("basic_enriched_patch")
        evidence.append(
            f"Basic patch: f+={f_pos:.2f}, NCPR={ncpr:.2f}, charge_dist={charge_label}"
        )

    # B. Acidic patch
    if (
        f_neg >= THRESHOLDS["f_minus_min"]
        and ncpr <= THRESHOLDS["ncpr_acidic_max"]
        and charge_label == "compact"
    ):
        labels.append("acidic_patch")
        evidence.append(
            f"Acidic patch: f-={f_neg:.2f}, NCPR={ncpr:.2f}, charge_dist={charge_label}"
        )

    # C. Mixed-charge polyampholyte
    if (
        fcr >= THRESHOLDS["fcr_polyampholyte_min"]
        and abs(ncpr) <= THRESHOLDS["ncpr_polyampholyte_abs_max"]
        and charge_label in ["dispersed", "periodic-ish"]
    ):
        labels.append("mixed_charge_polyampholyte")
        evidence.append(
            f"Polyampholyte: FCR={fcr:.2f}, |NCPR|={abs(ncpr):.2f}, charge_dist={charge_label}"
        )

    # D. Polar linker
    if (
        f_polar >= THRESHOLDS["f_polar_min"]
        and f_hydro <= THRESHOLDS["f_hydro_max_for_polar_linker"]
        and polar_label == "dispersed"
    ):
        labels.append("polar_linker")
        evidence.append(
            f"Polar linker: f_polar={f_polar:.2f}, f_hydro={f_hydro:.2f}, polar_dist={polar_label}"
        )

    # E. Hydrophobic patch
    if f_hydro >= THRESHOLDS["f_hydro_min"] and hydro_label == "compact":
        labels.append("hydrophobic_patch")
        evidence.append(
            f"Hydrophobic patch: f_hydro={f_hydro:.2f}, hydro_dist={hydro_label}"
        )

    # F. Mixed polar–hydrophobic segment
    if (
        (THRESHOLDS["f_polar_mixed_min"] <= f_polar <= THRESHOLDS["f_polar_mixed_max"])
        and (THRESHOLDS["f_hydro_mixed_min"] <= f_hydro <= THRESHOLDS["f_hydro_mixed_max"])
        and hydro_label in ["dispersed", "periodic-ish"]
        and polar_label in ["dispersed", "periodic-ish"]
    ):
        labels.append("mixed_polar_hydrophobic_segment")
        evidence.append(
            f"Mixed polar-hydro: f_polar={f_polar:.2f}, f_hydro={f_hydro:.2f}, "
            f"polar_dist={polar_label}, hydro_dist={hydro_label}"
        )

    # G. Aromatic-enriched patch
    if f_arom >= THRESHOLDS["f_aromatic_min"] and arom_label == "compact":
        labels.append("aromatic_enriched_patch")
        evidence.append(
            f"Aromatic patch: f_arom={f_arom:.2f}, arom_dist={arom_label}"
        )

    # H. Cation–π rich neighborhood
    cooc_key = "cooc_positive_aromatic"
    if cooc_quantiles is not None and cooc_key in cooc_quantiles:
        q = cooc[cooc_key]
        if q >= THRESHOLDS["cooc_top_quantile"]:
            labels.append("cation_pi_rich_neighborhood")
            evidence.append(
                f"Cation–π rich: cooc_positive_aromatic={q:.2f} (top {int((1-THRESHOLDS['cooc_top_quantile'])*100)}%)"
            )
    else:
        # Without quantiles, skip this rule in v1
        pass

    # I. Small-rich flexible segment
    if f_small >= THRESHOLDS["f_small_min"] and small_label == "dispersed":
        labels.append("small_rich_flexible_segment")
        evidence.append(
            f"Small-rich: f_small={f_small:.2f}, small_dist={small_label}"
        )

    # J. Bulky segment
    if f_small <= THRESHOLDS["f_small_max_for_bulky"]:
        labels.append("bulky_segment")
        evidence.append(f"Bulky segment: f_small={f_small:.2f}")

    # K. Chemically neutral linker (fallback if no other label)
    if len(labels) == 0:
        # Check that nothing is extremely enriched and distributions are dispersed
        extreme = (
            (f_polar > 0.6) or (f_hydro > 0.5) or (f_arom > 0.15) or
            (f_small > 0.6) or (f_pos > 0.3) or (f_neg > 0.3)
        )
        if not extreme:
            labels.append("chemically_neutral_linker")
            evidence.append("No strong compositional bias; all distributions dispersed.")

    if len(labels) == 0:
        # Safety fallback
        labels.append("unclassified")
        evidence.append("No behavior rule matched.")

    primary = labels[0]
    secondary = labels[1:] if len(labels) > 1 else []
    return primary, secondary, evidence


# =========================
# Main pipeline
# =========================

def main():
    parser = argparse.ArgumentParser(
        description="Compute LCR physicochemical features and behavior annotations."
    )
    parser.add_argument("--lcr_file", default=LCR_FILE, help="LCR Excel file")
    parser.add_argument("--aa_file", default=AA_PROP_FILE, help="AA properties CSV")
    parser.add_argument("--output_features", default=OUTPUT_FEATURES)
    parser.add_argument("--output_annotations", default=OUTPUT_ANNOTATIONS)
    args = parser.parse_args()

    print("Loading amino-acid properties...")
    aa_lookup = load_aa_properties(args.aa_file)

    print("Loading LCRs...")
    lcrs = load_lcrs(args.lcr_file, LCR_SHEET)

    print(f"Processing {len(lcrs)} LCRs...")

    # Pre-allocate feature dicts
    comp_rows = []
    distr_rows = []
    cooc_rows = []
    annot_rows = []

    # First pass: compute all features
    all_cooc = []

    for idx, row in lcrs.iterrows():
        seq = row["sequence"]
        if not isinstance(seq, str):
            seq = ""

        # Composition
        comp = compute_composition_features(seq, aa_lookup)
        comp_row = {
            "protein_id": row["protein_id"],
            "method": row["method"],
            "start": row["start"],
            "end": row["end"],
            "length": row["length"],
            **comp,
        }
        comp_rows.append(comp_row)

        # Distribution
        distr = compute_distribution_features(seq, aa_lookup)
        distr_row = {
            "protein_id": row["protein_id"],
            "method": row["method"],
            "start": row["start"],
            "end": row["end"],
            "length": row["length"],
            **distr,
        }
        distr_rows.append(distr_row)

        # Co-occurrence
        cooc = compute_cooccurrence_features(seq, aa_lookup)
        cooc_row = {
            "protein_id": row["protein_id"],
            "method": row["method"],
            "start": row["start"],
            "end": row["end"],
            "length": row["length"],
            **cooc,
        }
        cooc_rows.append(cooc_row)
        all_cooc.append(cooc)

    comp_df = pd.DataFrame(comp_rows)
    distr_df = pd.DataFrame(distr_rows)
    cooc_df = pd.DataFrame(cooc_rows)

    # Compute co-occurrence quantiles for cation–π rule
    cooc_quantiles = {}
    for key in ["cooc_positive_aromatic", "cooc_hydrophobic_hydrophobic",
                "cooc_polar_polar", "cooc_polar_hydrophobic", "cooc_positive_negative"]:
        if key in cooc_df.columns:
            vals = cooc_df[key].dropna()
            if len(vals) > 0:
                cooc_quantiles[key] = np.quantile(vals, [0.25, 0.5, 0.75])

    # Second pass: assign behavior labels
    print("Assigning behavior annotations...")
    for idx, row in lcrs.iterrows():
        seq = row["sequence"]
        if not isinstance(seq, str):
            seq = ""

        comp = compute_composition_features(seq, aa_lookup)
        distr = compute_distribution_features(seq, aa_lookup)
        cooc = compute_cooccurrence_features(seq, aa_lookup)

        primary, secondary, evidence = assign_behavior_labels(comp, distr, cooc, cooc_quantiles)

        annot_row = {
            "protein_id": row["protein_id"],
            "method": row["method"],
            "start": row["start"],
            "end": row["end"],
            "length": row["length"],
            "sequence": seq,
            "description": row.get("description", ""),
            "primary_physicochemical_annotation": primary,
            "secondary_physicochemical_annotations": ";".join(secondary) if secondary else "",
            "annotation_evidence": " | ".join(evidence),
        }
        annot_rows.append(annot_row)

    annot_df = pd.DataFrame(annot_rows)

    # Save outputs
    print("Saving feature workbook...")
    with pd.ExcelWriter(args.output_features, engine="openpyxl") as writer:
        comp_df.to_excel(writer, sheet_name="composition", index=False)
        distr_df.to_excel(writer, sheet_name="distribution", index=False)
        cooc_df.to_excel(writer, sheet_name="co_occurrence", index=False)

    print("Saving annotation workbook...")
    annot_df.to_excel(args.output_annotations, index=False)

    print("Done.")
    print(f"Features written to: {args.output_features}")
    print(f"Annotations written to: {args.output_annotations}")


if __name__ == "__main__":
    main()