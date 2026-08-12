#!/usr/bin/env python3
"""
LCR physicochemical feature extraction and behavior annotation (v2.2).

Changes from v2.1 (rule-logic fixes from output review):
  - Basic/acidic rules: dropped the charge-distribution ("compact")
    requirement. Neutral interrupters split charge tracks into many runs,
    causing poly-E tracts (f- ~0.8) and Arg-rich segments to be missed.
    NCPR sign+magnitude already encodes one-sign dominance at LCR scales.
  - hydrophobic_region now keyed on strong-hydrophobe fraction
    (V, I, L, M, F, W, Y) instead of the Taylor/Jalview hydrophobic set,
    which includes A, G, T and produced false "aggregation-prone" labels
    on Ala/Thr/Gly homopolymers (and contradictory co-annotations with
    gs_rich_neutral_region).
  - polar_linker now requires low charge (fcr <= 0.30): spacers in the
    stickers-and-spacers sense are UNCHARGED polar residues. Previously
    poly-E tracts (f_polar ~0.9) and Arg-rich RP segments were mislabeled
    as polar linkers.
  - gs_rich_neutral_region fraction path now requires BOTH components
    (min(frac_G, frac_S) >= 0.10) so pure poly-S tracts fall through to
    serine_rich_region. The GS-repeat-cluster path is unchanged (a GS|SG
    cluster contains both by definition).

Inputs:
    - lcr_methods_combined.xlsx  (sheet: "all_results")
    - aa-physicochemical-properties-v2.csv

Outputs:
    - lcr_features_v2.xlsx    (composition / distribution / co_occurrence sheets)
    - lcr_annotations_v2.xlsx (metadata + behavior annotations + evidence)

Literature basis:
    Das & Pappu (2013) PNAS 110:13392  -- FCR/NCPR, charge patterning
    Holehouse et al. (2017) Biophys J  -- localCIDER / CIDER toolkit
    Martin et al. (2020) Science 367:694 -- aromatic valence & patterning
    Holehouse et al. (2021) Biochemistry -- stickers & spacers, GS spacers
    Chong, Vernon & Forman-Kay (2018) JMB -- RGG/RG motifs in RNA binding
    Gerstberger et al. (2014) Nat Rev Genet -- RG/RGG repeat definition
    Xiang et al. (2013) Structure -- RS domain phosphorylation switch
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

OUTPUT_FEATURES = "lcr_features_v2_2.xlsx"
OUTPUT_ANNOTATIONS = "lcr_annotations_v2_2.xlsx"

# Consolidated thresholds (to be calibrated in a dedicated phase)
THRESHOLDS = {
    # Charge
    "f_charge_min": 0.25,        # enrichment of one charge sign
    "ncpr_min": 0.10,            # net charge bias for basic/acidic
    "frac_R_min": 0.12,          # Arg enrichment within basic regions
    "fcr_min": 0.30,             # polyampholyte (Das & Pappu 2013)
    "ncpr_neutral_max": 0.05,    # near-neutral net charge
    # Polarity / hydropathy
    "f_polar_min": 0.50,
    "fcr_polar_linker_max": 0.30,  # polar spacers are uncharged (v2.2)
    "f_hydro_low": 0.20,         # max hydrophobicity for polar linker
    "frac_strong_hydro_min": 0.30,  # V,I,L,M,F,W,Y fraction (v2.2)
    # Aromatic
    "f_arom_min": 0.08,
    # Disorder
    "f_disorder_min": 0.60,
    # Single-AA enrichment
    "frac_S_min": 0.30,
    "frac_Q_min": 0.20,
    "frac_G_min": 0.25,
    "frac_P_min": 0.15,
    # GS-rich neutral signature
    "f_GS_min": 0.40,            # combined Gly+Ser fraction
    "frac_GS_each_min": 0.10,    # both G and S must be present (v2.2)
    # Repeat detection (Gerstberger et al. 2014)
    "min_repeats": 3,
    "max_repeat_spacing": 10,
    # Co-occurrence
    "cooc_top_quantile": 0.75,
}

# Strong hydrophobes used for the aggregation-prone rule (v2.2)
STRONG_HYDRO = "VILMFWY"

# Motif regexes
RG_REGEX = r"RG{1,2}|RPR"   # RG, RGG, and Arg-Pro-Arg repeats
RS_REGEX = r"RS|SR"          # strict: Arg required (fixes SK/TK false positives)
GS_REGEX = r"GS|SG"


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
# Repeat detection
# =========================

def find_repeat_cluster(seq, motif_regex, min_repeats, max_spacing):
    """
    Robust repeat-region detection.

    Finds all non-overlapping motif matches, then checks for a cluster of
    >= min_repeats motifs where consecutive matches are spaced
    <= max_spacing residues apart (start-to-start), following the
    Gerstberger et al. (2014) RG/RGG definition.

    Returns:
        n_motifs   -- total non-overlapping matches in the sequence
        n_cluster  -- size of the largest qualifying cluster (0 if none)
        qualifies   -- True if any cluster meets the criterion
    """
    starts = [m.start() for m in re.finditer(motif_regex, seq)]
    n_motifs = len(starts)
    if n_motifs < min_repeats:
        return n_motifs, 0, False

    best = run = 1
    for i in range(1, n_motifs):
        if starts[i] - starts[i - 1] <= max_spacing:
            run += 1
            best = max(best, run)
        else:
            run = 1

    qualifies = best >= min_repeats
    return n_motifs, best if qualifies else 0, qualifies


def compute_repeat_features(seq):
    """Repeat/motif features for RG-like, RS-like, and GS-like patterns."""
    seq = seq.upper()

    n_rg, rg_cluster, rg_ok = find_repeat_cluster(
        seq, RG_REGEX, THRESHOLDS["min_repeats"], THRESHOLDS["max_repeat_spacing"]
    )
    n_sr, sr_cluster, sr_ok = find_repeat_cluster(
        seq, RS_REGEX, THRESHOLDS["min_repeats"], THRESHOLDS["max_repeat_spacing"]
    )
    n_gs, gs_cluster, gs_ok = find_repeat_cluster(
        seq, GS_REGEX, THRESHOLDS["min_repeats"], THRESHOLDS["max_repeat_spacing"]
    )

    return {
        "n_RG_motifs": n_rg,
        "RG_repeat_cluster": rg_cluster,
        "RG_repeat_region": rg_ok,
        "n_SR_motifs": n_sr,
        "SR_repeat_cluster": sr_cluster,
        "SR_repeat_region": sr_ok,
        "n_GS_motifs": n_gs,
        "GS_repeat_cluster": gs_cluster,
        "GS_repeat_region": gs_ok,
    }


# =========================
# Composition features
# =========================

def compute_composition_features(seq, aa_lookup):
    """Composition features for one sequence."""
    seq = seq.upper()
    L = len(seq)

    empty = {
        "length": 0,
        **{f"frac_{aa}": 0.0 for aa in "ACDEFGHIKLMNPQRSTVWY"},
        "frac_polar": 0.0,
        "frac_hydrophobic": 0.0,
        "frac_strong_hydro": 0.0,
        "frac_aromatic": 0.0,
        "frac_disorder": 0.0,
        "frac_positive": 0.0,
        "frac_negative": 0.0,
        "fcr": 0.0,
        "ncpr": 0.0,
        "frac_GS": 0.0,
        **compute_repeat_features(""),
    }
    if L == 0:
        return empty

    aa_counts = defaultdict(int)
    for res in seq:
        if res in "ACDEFGHIKLMNPQRSTVWY":
            aa_counts[res] += 1

    features = {"length": L}
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        features[f"frac_{aa}"] = aa_counts[aa] / L

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

    features["frac_polar"] = n_polar / L
    features["frac_hydrophobic"] = n_hydro / L
    features["frac_strong_hydro"] = sum(aa_counts[a] for a in STRONG_HYDRO) / L
    features["frac_aromatic"] = n_arom / L
    features["frac_disorder"] = n_disorder / L
    features["frac_positive"] = n_pos / L
    features["frac_negative"] = n_neg / L
    features["fcr"] = (n_pos + n_neg) / L
    features["ncpr"] = (n_pos - n_neg) / L
    features["frac_GS"] = (aa_counts["G"] + aa_counts["S"]) / L

    features.update(compute_repeat_features(seq))

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
    """Cation-pi co-occurrence score with length-adaptive window."""
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
    Assign behavior labels. Rule order matters; first match = primary.

    Signatures (v2.2):
      Repeat:    rg_rgg_repeat_region, sr_rs_repeat_region,
                 gs_rich_neutral_region
      Charge:    basic_enriched_region, arginine_rich_rna_contact_region,
                 acidic_region, mixed_charge_polyampholyte
      Hydropathy: polar_linker, hydrophobic_region
      Aromatic:  aromatic_sticker_region, aromatic_aggregation_prone_region,
                 cation_pi_rich_neighborhood
      Disorder:  disorder_rich_spacer_region
      Single-AA: serine_rich_region, glutamine_rich_region,
                 glycine_rich_region, proline_rich_disordered_region
      Fallback:  unmapped_physicochemical_properties
    """
    labels = []
    evidence = []

    f_pos = comp["frac_positive"]
    f_neg = comp["frac_negative"]
    ncpr = comp["ncpr"]
    fcr = comp["fcr"]
    f_polar = comp["frac_polar"]
    f_hydro = comp["frac_hydrophobic"]
    f_sh = comp["frac_strong_hydro"]
    f_arom = comp["frac_aromatic"]
    f_disorder = comp["frac_disorder"]
    f_GS = comp["frac_GS"]

    charge_label = distr.get("charged_label", "insufficient")
    polar_label = distr.get("polar_label", "insufficient")
    hydro_label = distr.get("hydrophobic_label", "insufficient")
    arom_label = distr.get("aromatic_label", "insufficient")

    # ---- Repeat-driven signatures (highest specificity first) ----

    # 1. RG/RGG/RPR repeat region
    # Chong et al. (2018); Gerstberger et al. (2014) repeat definition
    if comp["RG_repeat_region"]:
        labels.append("rg_rgg_repeat_region")
        evidence.append(
            f"RG/RGG/RPR repeats: cluster of {comp['RG_repeat_cluster']} motifs "
            f"(spacing <= {THRESHOLDS['max_repeat_spacing']} aa); "
            f"RNA-binding and phase-separation driver"
        )

    # 2. SR/RS repeat region
    # Xiang et al. (2013): phosphorylation-gated RNA interaction switch
    if comp["SR_repeat_region"]:
        labels.append("sr_rs_repeat_region")
        evidence.append(
            f"RS/SR repeats: cluster of {comp['SR_repeat_cluster']} motifs; "
            f"phosphorylation-sensitive RNA-binding"
        )

    # 3. GS-rich neutral region
    # Holehouse et al. (2021): GS repeats ~ ideal Gaussian chains.
    # Fraction path requires BOTH G and S present (v2.2): pure poly-S
    # tracts fall through to serine_rich_region.
    gs_cluster_ok = comp["GS_repeat_region"]
    gs_fraction_ok = (
        f_GS >= THRESHOLDS["f_GS_min"]
        and min(comp["frac_G"], comp["frac_S"]) >= THRESHOLDS["frac_GS_each_min"]
    )
    if gs_cluster_ok or gs_fraction_ok:
        labels.append("gs_rich_neutral_region")
        evidence.append(
            f"GS-rich: f_GS={f_GS:.2f}, GS-repeat cluster={comp['GS_repeat_cluster']}; "
            f"chemically neutral spacer (~ideal chain)"
        )

    # ---- Charge-driven signatures ----

    # 4. Basic region (Arg-enriched gets specific label)
    # Das & Pappu (2013); Chong et al. (2018) Arg-RNA mechanism.
    # v2.2: no charge-distribution requirement -- neutral interrupters
    # split charge tracks into runs and caused false negatives.
    if f_pos >= THRESHOLDS["f_charge_min"] and ncpr >= THRESHOLDS["ncpr_min"]:
        if comp["frac_R"] >= THRESHOLDS["frac_R_min"]:
            labels.append("arginine_rich_rna_contact_region")
            evidence.append(
                f"Basic Arg-rich region: f_R={comp['frac_R']:.2f}, "
                f"f+={f_pos:.2f}, NCPR={ncpr:.2f}; RNA-contact"
            )
        else:
            labels.append("basic_enriched_region")
            evidence.append(
                f"Basic region: f+={f_pos:.2f}, NCPR={ncpr:.2f}; "
                f"expanded-coil polyelectrolyte"
            )

    # 5. Acidic region (v2.2: no charge-distribution requirement)
    if f_neg >= THRESHOLDS["f_charge_min"] and ncpr <= -THRESHOLDS["ncpr_min"]:
        labels.append("acidic_region")
        evidence.append(
            f"Acidic region: f-={f_neg:.2f}, NCPR={ncpr:.2f}; functional module"
        )

    # 6. Mixed-charge polyampholyte
    # Das & Pappu (2013) strong polyampholyte regime
    if (
        fcr >= THRESHOLDS["fcr_min"]
        and abs(ncpr) <= THRESHOLDS["ncpr_neutral_max"]
        and charge_label == "dispersed"
    ):
        labels.append("mixed_charge_polyampholyte")
        evidence.append(
            f"Polyampholyte: FCR={fcr:.2f}, |NCPR|={abs(ncpr):.2f}; "
            f"salt-tunable chain dimensions"
        )

    # ---- Polarity / hydropathy signatures ----

    # 7. Polar linker (uncharged polar spacer)
    # v2.2: FCR cap -- stickers-and-spacers polar spacers are uncharged;
    # previously poly-E tracts and Arg-rich segments were mislabeled.
    if (
        f_polar >= THRESHOLDS["f_polar_min"]
        and f_hydro <= THRESHOLDS["f_hydro_low"]
        and fcr <= THRESHOLDS["fcr_polar_linker_max"]
        and polar_label == "dispersed"
    ):
        labels.append("polar_linker")
        evidence.append(
            f"Polar linker: f_polar={f_polar:.2f}, FCR={fcr:.2f}; uncharged spacer"
        )

    # 8. Hydrophobic region (aggregation-prone, atypical for LCRs)
    # v2.2: keyed on strong hydrophobes (VILMFWY); Taylor/Jalview set
    # includes A, G, T and produced false positives on homopolymers.
    if f_sh >= THRESHOLDS["frac_strong_hydro_min"] and hydro_label == "compact":
        labels.append("hydrophobic_region")
        evidence.append(
            f"Hydrophobic region: f_strong_hydro={f_sh:.2f}, compact; "
            f"aggregation-prone"
        )

    # ---- Aromatic / sticker signatures ----

    # 9a. Aromatic sticker region (dispersed -> LLPS)
    # Martin et al. (2020)
    if f_arom >= THRESHOLDS["f_arom_min"] and arom_label == "dispersed":
        labels.append("aromatic_sticker_region")
        evidence.append(
            f"Aromatic stickers: f_arom={f_arom:.2f}, dispersed; "
            f"promotes LLPS, inhibits aggregation"
        )

    # 9b. Aromatic aggregation-prone region (compact)
    if f_arom >= THRESHOLDS["f_arom_min"] and arom_label == "compact":
        labels.append("aromatic_aggregation_prone_region")
        evidence.append(
            f"Aromatic cluster: f_arom={f_arom:.2f}, compact; aggregation risk"
        )

    # 10. Cation-pi rich neighborhood
    if cooc_quantile is not None:
        q = cooc["cooc_positive_aromatic"]
        if q >= cooc_quantile:
            labels.append("cation_pi_rich_neighborhood")
            evidence.append(
                f"Cation-pi rich: cooc={q:.2f} (dataset top-25% cutoff "
                f"{cooc_quantile:.2f}); ligand/nucleic-acid binding"
            )

    # ---- Disorder signature ----

    # 11. Disorder-rich spacer region
    if f_disorder >= THRESHOLDS["f_disorder_min"]:
        labels.append("disorder_rich_spacer_region")
        evidence.append(f"Disorder-promoting region: f_disorder={f_disorder:.2f}")

    # ---- Single-AA enrichment signatures ----

    # 12. Serine-rich region
    if comp["frac_S"] >= THRESHOLDS["frac_S_min"]:
        labels.append("serine_rich_region")
        evidence.append(f"Serine-rich: f_S={comp['frac_S']:.2f}")

    # 13. Glutamine-rich region
    if comp["frac_Q"] >= THRESHOLDS["frac_Q_min"]:
        labels.append("glutamine_rich_region")
        evidence.append(f"Glutamine-rich: f_Q={comp['frac_Q']:.2f}")

    # 14. Glycine-rich region
    if comp["frac_G"] >= THRESHOLDS["frac_G_min"]:
        labels.append("glycine_rich_region")
        evidence.append(f"Glycine-rich: f_G={comp['frac_G']:.2f}")

    # 15. Proline-rich disordered region
    if comp["frac_P"] >= THRESHOLDS["frac_P_min"]:
        labels.append("proline_rich_disordered_region")
        evidence.append(f"Proline-rich: f_P={comp['frac_P']:.2f}")

    # ---- Fallback ----

    if len(labels) == 0:
        labels.append("unmapped_physicochemical_properties")
        evidence.append("No behavior signature matched")

    primary = labels[0]
    secondary = labels[1:] if len(labels) > 1 else []
    return primary, secondary, evidence


# =========================
# Main pipeline
# =========================

def main():
    parser = argparse.ArgumentParser(
        description="LCR physicochemical features and behavior annotations (v2.2)."
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

    comp_df = pd.DataFrame(comp_rows)
    distr_df = pd.DataFrame(distr_rows)
    cooc_df = pd.DataFrame(cooc_rows)

    # 75th percentile threshold for cation-pi co-occurrence
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
