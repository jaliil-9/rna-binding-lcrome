#!/usr/bin/env python3
"""
LCR physicochemical feature extraction and behavior annotation (v2).

Changes from v1:
  - "Small?" property replaced with "Disorder?" (TOP-IDP disorder-promoting set)
  - RG/RGG and SR/RS motif counting added (regex-based)
  - Single-AA enrichment rules added (Ser, Gln, Gly, Pro)
  - Aromatic patch rule corrected: dispersed aromatics = LLPS stickers,
    compact aromatics = aggregation-prone (Martin et al. 2020)
  - Bulky segment rule dropped (not literature-supported)
  - Co-occurrence window made length-adaptive; polar-polar and
    hydrophobic-hydrophobic pairs dropped (no added signal)
  - Distribution labels simplified: compact / dispersed / insufficient
    (dropped fragile "periodic-ish")
  - Output schema includes rna_target_superclass and
    domain_position_class placeholder columns

Inputs:
    - lcr_methods_combined.xlsx  (sheet: "all_results")
    - aa-physicochemical-properties-v2.csv

Outputs:
    - lcr_features_v2.xlsx    (composition / distribution / co_occurrence sheets)
    - lcr_annotations_v2.xlsx (metadata + behavior annotations + evidence)

Literature basis:
    Das & Pappu (2013) PNAS 110:13392  -- charge patterning, diagram of states
    Holehouse et al. (2017) Biophys J  -- localCIDER / CIDER toolkit
    Martin et al. (2020) Science 367:694 -- aromatic valence & patterning
    Holehouse et al. (2021) Biochemistry -- stickers & spacers framework
    Chong, Vernon & Forman-Kay (2018) JMB -- RGG/RG motifs in RNA binding
    Shepard & Hertel (2009) Genome Biol -- SR/RS domain characterization
"""

import pandas as pd
import numpy as np
from collections import defaultdict
import argparse
import re


# =========================
# Configuration
# =========================

LCR_FILE = r"rbp_lcrs\lcr_methods_combined.xlsx"
AA_PROP_FILE = r"lcr_analyses\pc_properties\aa-physicochemical-properties-v2.csv"
LCR_SHEET = "all_results"

OUTPUT_FEATURES = "lcr_features_v2.xlsx"
OUTPUT_ANNOTATIONS = "lcr_annotations_v2.xlsx"

# Behavior-rule thresholds (v2)
THRESHOLDS = {
    # Charge-driven
    "f_plus_min": 0.25,
    "ncpr_basic_min": 0.10,
    "frac_R_min": 0.12,
    "f_minus_min": 0.25,
    "ncpr_acidic_max": -0.10,
    "fcr_polyampholyte_min": 0.30,
    "ncpr_polyampholyte_abs_max": 0.05,
    # Polarity / hydropathy
    "f_polar_min": 0.50,
    "f_hydro_max_for_polar_linker": 0.20,
    "f_hydro_min": 0.40,
    "f_polar_mixed_min": 0.25,
    "f_polar_mixed_max": 0.55,
    "f_hydro_mixed_min": 0.25,
    "f_hydro_mixed_max": 0.55,
    # Aromatic / stickers
    "f_aromatic_min": 0.08,
    # Disorder
    "f_disorder_min": 0.60,
    # Single-AA enrichment
    "frac_S_min": 0.30,
    "frac_Q_min": 0.20,
    "frac_G_min": 0.25,
    "frac_P_min": 0.15,
    # Motif counting
    "min_RG_repeats": 3,
    "min_SR_repeats": 3,
    # Co-occurrence
    "cooc_top_quantile": 0.75,
}


# =========================
# Load data
# =========================

def load_aa_properties(filepath):
    """Load AA property CSV (v2 with Disorder? column)."""
    aa_df = pd.read_csv(filepath)
    aa_df.columns = aa_df.columns.str.strip()
    aa_df["AA_letter"] = aa_df["AA"].str.extract(r"^\s*([A-Z])\s*\(")[0]

    aa_lookup = {}
    for _, row in aa_df.iterrows():
        aa = row["AA_letter"]
        if pd.isna(aa):
            continue
        aa_lookup[aa] = {
            "charge": row["Charge (7.4)"],
            "polar": row["Polar?"] == "Yes",
            "hydrophobic": row["Hydrophobic?"] == "Yes",
            "aromatic": row["Aromatic?"] == "Yes",
            "disorder": row["Disorder?"] == "Yes",
        }
    return aa_lookup


def load_lcrs(filepath, sheet_name):
    """Load LCR table."""
    df = pd.read_excel(filepath, sheet_name=sheet_name)
    required_cols = ["protein_id", "method", "start", "end", "length", "sequence"]
    for c in required_cols:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")
    return df


# =========================
# Sequence → property tracks
# =========================

def encode_sequence(seq, aa_lookup):
    """Return dict of binary property tracks and signed charge track."""
    seq = seq.upper()
    L = len(seq)

    polar, hydro, aromatic, disorder, charge = [], [], [], [], []

    for res in seq:
        if res not in aa_lookup:
            polar.append(0)
            hydro.append(0)
            aromatic.append(0)
            disorder.append(0)
            charge.append(0)
            continue

        props = aa_lookup[res]
        polar.append(1 if props["polar"] else 0)
        hydro.append(1 if props["hydrophobic"] else 0)
        aromatic.append(1 if props["aromatic"] else 0)
        disorder.append(1 if props["disorder"] else 0)

        ch = props["charge"]
        if ch == "Pos":
            charge.append(1)
        elif ch == "Neg":
            charge.append(-1)
        else:
            charge.append(0)

    tracks = {
        "polar": np.array(polar, dtype=int),
        "hydrophobic": np.array(hydro, dtype=int),
        "aromatic": np.array(aromatic, dtype=int),
        "disorder": np.array(disorder, dtype=int),
        "charged": np.array([1 if c != 0 else 0 for c in charge], dtype=int),
        "charge_signed": np.array(charge, dtype=int),
    }
    return tracks, L


# =========================
# Composition features
# =========================

def count_motifs(seq):
    """Count RG/RGG and SR/RS motif occurrences."""
    rg_matches = re.findall(r"(?:RG{1,2}|RGG)", seq.upper())
    sr_matches = re.findall(r"[ST][RK]", seq.upper())
    return len(rg_matches), len(sr_matches)


def compute_composition_features(seq, aa_lookup):
    """Composition features for one sequence."""
    seq = seq.upper()
    L = len(seq)

    empty = {
        "length": 0,
        **{f"frac_{aa}": 0.0 for aa in "ACDEFGHIKLMNPQRSTVWY"},
        "frac_polar": 0.0,
        "frac_hydrophobic": 0.0,
        "frac_aromatic": 0.0,
        "frac_disorder": 0.0,
        "frac_positive": 0.0,
        "frac_negative": 0.0,
        "fcr": 0.0,
        "ncpr": 0.0,
        "n_RG_motifs": 0,
        "n_SR_motifs": 0,
        "frac_R": 0.0,
        "frac_K": 0.0,
    }
    if L == 0:
        return empty

    aa_counts = defaultdict(int)
    for res in seq:
        if res in "ACDEFGHIKLMNPQRSTVWY":
            aa_counts[res] += 1

    features = {"length": L}
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        features[f"frac_{aa}"] = aa_counts[aa] / L # type: ignore

    n_polar = n_hydro = n_arom = n_disorder = n_pos = n_neg = 0

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
        if p["disorder"]:
            n_disorder += 1
        ch = p["charge"]
        if ch == "Pos":
            n_pos += 1
        elif ch == "Neg":
            n_neg += 1

    features["frac_polar"] = n_polar / L # type: ignore
    features["frac_hydrophobic"] = n_hydro / L # type: ignore
    features["frac_aromatic"] = n_arom / L # type: ignore
    features["frac_disorder"] = n_disorder / L # type: ignore
    features["frac_positive"] = n_pos / L # type: ignore
    features["frac_negative"] = n_neg / L # type: ignore
    features["fcr"] = (n_pos + n_neg) / L # type: ignore
    features["ncpr"] = (n_pos - n_neg) / L # type: ignore
    features["frac_R"] = features["frac_R"]  # already set above
    features["frac_K"] = features["frac_K"]

    n_rg, n_sr = count_motifs(seq)
    features["n_RG_motifs"] = n_rg
    features["n_SR_motifs"] = n_sr

    return features


# =========================
# Distribution features
# =========================

def runs_and_gaps(binary_array):
    """Run/gap statistics; simplified labels: compact / dispersed / insufficient."""
    if binary_array.sum() == 0:
        return {
            "n_runs": 0,
            "mean_run_length": np.nan,
            "max_run_length": 0,
            "mean_gap": np.nan,
            "cv_gap": np.nan,
            "label": "insufficient",
        }

    diff = np.diff(np.concatenate(([0], binary_array, [0])))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    run_lengths = ends - starts
    n_runs = len(run_lengths)
    mean_run = float(np.mean(run_lengths))
    max_run = int(np.max(run_lengths))

    if n_runs <= 1:
        mean_gap = np.nan
        cv_gap = np.nan
        label = "compact" if max_run / len(binary_array) > 0.4 else "dispersed"
    else:
        gaps = starts[1:] - ends[:-1]
        mean_gap = float(np.mean(gaps))
        cv_gap = float(np.std(gaps) / mean_gap) if mean_gap > 0 else 0.0

        frac_max = max_run / len(binary_array)
        if n_runs <= 3 and frac_max > 0.35:
            label = "compact"
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
    """Distribution features for each property track."""
    tracks, L = encode_sequence(seq, aa_lookup)
    prop_names = ["polar", "hydrophobic", "aromatic", "disorder", "charged"]

    features = {}
    for prop in prop_names:
        stats = runs_and_gaps(tracks[prop])
        for k, v in stats.items():
            features[f"{prop}_{k}"] = v

    return features


# =========================
# Co-occurrence features
# =========================

def compute_cooccurrence_features(seq, aa_lookup):
    """
    Cation-pi co-occurrence score with length-adaptive window.
    Only the positive-aromatic pair is retained in v2 (the others added
    no signal beyond what distribution features already capture).
    """
    tracks, L = encode_sequence(seq, aa_lookup)

    if L == 0:
        return {"cooc_positive_aromatic": 0.0}

    window = max(5, L // 4)

    pos = tracks["charge_signed"] > 0
    arom = tracks["aromatic"].astype(bool)

    if L < window:
        return {"cooc_positive_aromatic": 0.0}

    n_windows = L - window + 1
    count = 0
    for i in range(n_windows):
        if pos[i : i + window].any() and arom[i : i + window].any():
            count += 1

    return {"cooc_positive_aromatic": count / n_windows}


# =========================
# Behavior annotation rules
# =========================

def assign_behavior_labels(comp, distr, cooc, cooc_quantile=None):
    """
    Assign behavior labels based on composition, distribution, and
    co-occurrence features. Rule order matters; first match = primary.

    Literature basis for each rule is noted inline.
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
    f_disorder = comp["frac_disorder"]
    f_S = comp["frac_S"]
    f_Q = comp["frac_Q"]
    f_G = comp["frac_G"]
    f_P = comp["frac_P"]
    f_R = comp["frac_R"]
    n_rg = comp["n_RG_motifs"]
    n_sr = comp["n_SR_motifs"]

    charge_label = distr.get("charged_label", "insufficient")
    polar_label = distr.get("polar_label", "insufficient")
    hydro_label = distr.get("hydrophobic_label", "insufficient")
    arom_label = distr.get("aromatic_label", "insufficient")
    disorder_label = distr.get("disorder_label", "insufficient")

    # ---- Motif-driven signatures (highest specificity first) ----

    # 1. RG/RGG repeat region
    # Chong, Vernon & Forman-Kay (2018); >1000 human RBPs, G4-RNA binding
    if n_rg >= THRESHOLDS["min_RG_repeats"]:
        labels.append("rg_rgg_repeat_region")
        evidence.append(
            f"RG/RGG repeats: {n_rg} motifs; RNA-binding and phase-separation driver"
        )

    # 2. SR/RS repeat region
    # Shepard & Hertel (2009); phosphorylation-gated RNA interaction switch
    if n_sr >= THRESHOLDS["min_SR_repeats"]:
        labels.append("sr_rs_repeat_region")
        evidence.append(
            f"SR/RS repeats: {n_sr} motifs; phosphorylation-sensitive RNA-binding"
        )

    # ---- Charge-driven signatures ----

    # 3. Basic patch (Arg-enriched)
    # Das & Pappu (2013) strong polyelectrolyte; Arg bidentate H-bonds to RNA
    if (
        f_pos >= THRESHOLDS["f_plus_min"]
        and ncpr >= THRESHOLDS["ncpr_basic_min"]
        and charge_label == "compact"
    ):
        if f_R >= THRESHOLDS["frac_R_min"]:
            labels.append("arginine_rich_rna_contact_patch")
            evidence.append(
                f"Arg-rich patch: f_R={f_R:.2f}, f+={f_pos:.2f}, NCPR={ncpr:.2f}"
            )
        else:
            labels.append("basic_enriched_patch")
            evidence.append(
                f"Basic patch: f+={f_pos:.2f}, NCPR={ncpr:.2f}, charge_dist={charge_label}"
            )

    # 4. Acidic patch
    if (
        f_neg >= THRESHOLDS["f_minus_min"]
        and ncpr <= THRESHOLDS["ncpr_acidic_max"]
        and charge_label == "compact"
    ):
        labels.append("acidic_patch")
        evidence.append(
            f"Acidic patch: f-={f_neg:.2f}, NCPR={ncpr:.2f}, charge_dist={charge_label}"
        )

    # 5. Mixed-charge polyampholyte
    # Das & Pappu (2013) strong polyampholyte regime
    if (
        fcr >= THRESHOLDS["fcr_polyampholyte_min"]
        and abs(ncpr) <= THRESHOLDS["ncpr_polyampholyte_abs_max"]
        and charge_label == "dispersed"
    ):
        labels.append("mixed_charge_polyampholyte")
        evidence.append(
            f"Polyampholyte: FCR={fcr:.2f}, |NCPR|={abs(ncpr):.2f}, "
            f"charge_dist={charge_label}"
        )

    # ---- Polarity / hydropathy signatures ----

    # 6. Polar linker (spacer signature)
    # Stickers-and-spacers: polar residues act as spacers between stickers
    if (
        f_polar >= THRESHOLDS["f_polar_min"]
        and f_hydro <= THRESHOLDS["f_hydro_max_for_polar_linker"]
        and polar_label == "dispersed"
    ):
        labels.append("polar_linker")
        evidence.append(
            f"Polar linker: f_polar={f_polar:.2f}, f_hydro={f_hydro:.2f}"
        )

    # 7. Hydrophobic patch (aggregation-prone)
    if f_hydro >= THRESHOLDS["f_hydro_min"] and hydro_label == "compact":
        labels.append("hydrophobic_patch")
        evidence.append(
            f"Hydrophobic patch: f_hydro={f_hydro:.2f}, hydro_dist={hydro_label}"
        )

    # 8. Mixed polar-hydrophobic segment
    if (
        (THRESHOLDS["f_polar_mixed_min"] <= f_polar <= THRESHOLDS["f_polar_mixed_max"])
        and (THRESHOLDS["f_hydro_mixed_min"] <= f_hydro <= THRESHOLDS["f_hydro_mixed_max"])
    ):
        labels.append("mixed_polar_hydrophobic_segment")
        evidence.append(
            f"Mixed polar-hydro: f_polar={f_polar:.2f}, f_hydro={f_hydro:.2f}"
        )

    # ---- Aromatic / sticker signatures ----

    # 9a. Aromatic sticker region (dispersed aromatics → LLPS)
    # Martin et al. (2020): uniform aromatic distribution drives phase separation
    if f_arom >= THRESHOLDS["f_aromatic_min"] and arom_label == "dispersed":
        labels.append("aromatic_sticker_region")
        evidence.append(
            f"Aromatic stickers: f_arom={f_arom:.2f}, arom_dist={arom_label}; "
            f"uniform distribution supports LLPS (Martin 2020)"
        )

    # 9b. Aromatic aggregation-prone patch (compact aromatics)
    # Martin et al. (2020): clustered aromatics drive aggregation, not LLPS
    if f_arom >= THRESHOLDS["f_aromatic_min"] and arom_label == "compact":
        labels.append("aromatic_aggregation_prone_patch")
        evidence.append(
            f"Aromatic cluster: f_arom={f_arom:.2f}, arom_dist={arom_label}; "
            f"compact patterning suggests aggregation risk"
        )

    # 10. Cation-pi rich neighborhood
    if cooc_quantile is not None:
        q = cooc["cooc_positive_aromatic"]
        if q >= cooc_quantile:
            labels.append("cation_pi_rich_neighborhood")
            evidence.append(
                f"Cation-pi rich: cooc={q:.2f} (top 25% dataset threshold={cooc_quantile:.2f})"
            )

    # ---- Disorder signature ----

    # 11. Disorder-rich spacer region
    # TOP-IDP disorder-promoting set; CIDER/LLPS literature
    if f_disorder >= THRESHOLDS["f_disorder_min"]:
        labels.append("disorder_rich_spacer_region")
        evidence.append(
            f"Disorder-rich: f_disorder={f_disorder:.2f}, disorder_dist={disorder_label}"
        )

    # ---- Single-AA enrichment signatures ----

    # 12. Serine-rich region
    if f_S >= THRESHOLDS["frac_S_min"]:
        labels.append("serine_rich_region")
        evidence.append(f"Serine-rich: f_S={f_S:.2f}")

    # 13. Glutamine-rich region
    if f_Q >= THRESHOLDS["frac_Q_min"]:
        labels.append("glutamine_rich_region")
        evidence.append(f"Glutamine-rich: f_Q={f_Q:.2f}")

    # 14. Glycine-rich region
    if f_G >= THRESHOLDS["frac_G_min"]:
        labels.append("glycine_rich_region")
        evidence.append(f"Glycine-rich: f_G={f_G:.2f}")

    # 15. Proline-rich disordered region
    if f_P >= THRESHOLDS["frac_P_min"]:
        labels.append("proline_rich_disordered_region")
        evidence.append(f"Proline-rich: f_P={f_P:.2f}")

    # ---- Fallback ----

    # 16. Chemically neutral linker
    if len(labels) == 0:
        extreme = (
            (f_polar > 0.6) or (f_hydro > 0.5) or (f_arom > 0.15)
            or (f_disorder > 0.7) or (f_pos > 0.3) or (f_neg > 0.3)
        )
        if not extreme:
            labels.append("chemically_neutral_linker")
            evidence.append("No strong compositional bias; generic flexible connector")

    if len(labels) == 0:
        labels.append("unclassified")
        evidence.append("No behavior rule matched")

    primary = labels[0]
    secondary = labels[1:] if len(labels) > 1 else []
    return primary, secondary, evidence


# =========================
# Main pipeline
# =========================

def main():
    parser = argparse.ArgumentParser(
        description="LCR physicochemical features and behavior annotations (v2)."
    )
    parser.add_argument("--lcr_file", default=LCR_FILE)
    parser.add_argument("--aa_file", default=AA_PROP_FILE)
    parser.add_argument("--output_features", default=OUTPUT_FEATURES)
    parser.add_argument("--output_annotations", default=OUTPUT_ANNOTATIONS)
    args = parser.parse_args()

    print("Loading amino-acid properties (v2)...")
    aa_lookup = load_aa_properties(args.aa_file)

    print("Loading LCRs...")
    lcrs = load_lcrs(args.lcr_file, LCR_SHEET)
    print(f"Processing {len(lcrs)} LCRs...")

    comp_rows = []
    distr_rows = []
    cooc_rows = []
    annot_rows = []

    all_cooc = []

    for idx, row in lcrs.iterrows():
        seq = row["sequence"]
        if not isinstance(seq, str):
            seq = ""

        base = {
            "protein_id": row["protein_id"],
            "method": row["method"],
            "start": row["start"],
            "end": row["end"],
            "length": row["length"],
        }

        comp = compute_composition_features(seq, aa_lookup)
        comp_rows.append({**base, **comp})

        distr = compute_distribution_features(seq, aa_lookup)
        distr_rows.append({**base, **distr})

        cooc = compute_cooccurrence_features(seq, aa_lookup)
        cooc_rows.append({**base, **cooc})
        all_cooc.append(cooc)

    comp_df = pd.DataFrame(comp_rows)
    distr_df = pd.DataFrame(distr_rows)
    cooc_df = pd.DataFrame(cooc_rows)

    # Compute 75th percentile threshold for cation-pi co-occurrence
    cooc_vals = cooc_df["cooc_positive_aromatic"].dropna()
    cooc_quantile = (
        float(np.quantile(cooc_vals, THRESHOLDS["cooc_top_quantile"]))
        if len(cooc_vals) > 0
        else None
    )

    print("Assigning behavior annotations...")
    for idx, row in lcrs.iterrows():
        seq = row["sequence"]
        if not isinstance(seq, str):
            seq = ""

        comp = compute_composition_features(seq, aa_lookup)
        distr = compute_distribution_features(seq, aa_lookup)
        cooc = compute_cooccurrence_features(seq, aa_lookup)

        primary, secondary, evidence = assign_behavior_labels(
            comp, distr, cooc, cooc_quantile
        )

        annot_rows.append({
            "protein_id": row["protein_id"],
            "method": row["method"],
            "start": row["start"],
            "end": row["end"],
            "length": row["length"],
            "sequence": seq,
            "description": row.get("description", ""),
            "rna_target_superclass": row.get("rna_target_superclass", ""),
            "domain_position_class": row.get("domain_position_class", ""),
            "primary_physicochemical_annotation": primary,
            "secondary_physicochemical_annotations": (
                ";".join(secondary) if secondary else ""
            ),
            "annotation_evidence": " | ".join(evidence),
        })

    annot_df = pd.DataFrame(annot_rows)

    print("Saving feature workbook...")
    with pd.ExcelWriter(args.output_features, engine="openpyxl") as writer:
        comp_df.to_excel(writer, sheet_name="composition", index=False)
        distr_df.to_excel(writer, sheet_name="distribution", index=False)
        cooc_df.to_excel(writer, sheet_name="co_occurrence", index=False)

    print("Saving annotation workbook...")
    annot_df.to_excel(args.output_annotations, index=False)

    print(f"Done.\nFeatures -> {args.output_features}")
    print(f"Annotations -> {args.output_annotations}")


if __name__ == "__main__":
    main()
