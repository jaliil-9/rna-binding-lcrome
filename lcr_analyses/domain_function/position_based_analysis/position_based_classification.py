#!/usr/bin/env python3
"""Classify method-specific LCRs by their positions relative to retained Pfam domains.

Inputs
------
1. lcr_methods_combined.xlsx, sheet: all_results
2. pfam_lcr_overlap.xlsx, sheet: Pfam_hits

The script treats every LCR call independently; it does not merge calls or use
multi-method support. All retained Pfam domains are assumed RNA-related.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

# ---- Classification settings (1-based, inclusive protein coordinates) ----
INTRINSIC_MIN_FRACTION = 0.70
LINKER_MAX_FRACTION = 0.30
ADJACENT_MAX_DISTANCE = 50
EDGE_MARGIN = 10


def accession(value: str) -> str:
    """Return the accession before the first '|' so input identifiers match."""
    return str(value).split("|")[0]


def overlap_length(a_start, a_end, b_start, b_end) -> int:
    """Number of shared residues between two inclusive coordinate intervals."""
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def interval_distance(a_start, a_end, b_start, b_end) -> int:
    """Residues between two non-overlapping intervals; zero if they touch/overlap."""
    if a_end < b_start:
        return b_start - a_end - 1
    if b_end < a_start:
        return a_start - b_end - 1
    return 0


def is_domain_edge(lcr_start, lcr_end, domain_start, domain_end, overlap) -> bool:
    """True when an overlapping LCR crosses or lies within EDGE_MARGIN of a domain end."""
    if overlap == 0:
        return False
    crosses_boundary = (lcr_start < domain_start <= lcr_end) or (lcr_start <= domain_end < lcr_end)
    near_boundary = (
        abs(lcr_start - domain_start) <= EDGE_MARGIN
        or abs(lcr_end - domain_end) <= EDGE_MARGIN
    )
    return crosses_boundary or near_boundary


def classify_lcr(lcr, protein_domains: pd.DataFrame) -> dict:
    """Return one primary class and all applicable position-based class labels."""
    start, end = int(lcr.lcr_start), int(lcr.lcr_end)
    lcr_length = end - start + 1

    if protein_domains.empty:
        return {
            "primary_class": "unclassified_no_pfam",
            "class_labels": "unclassified_no_pfam",
            "max_lcr_domain_fraction": 0.0,
            "nearest_domain_distance": pd.NA,
            "nearest_pfam_accession": pd.NA,
            "nearest_pfam_name": pd.NA,
            "flanking_left_pfam": pd.NA,
            "flanking_right_pfam": pd.NA,
        }

    domains = protein_domains.copy()
    domains["overlap_length"] = domains.apply(
        lambda d: overlap_length(start, end, d.domain_start, d.domain_end), axis=1
    )
    domains["lcr_domain_fraction"] = domains["overlap_length"] / lcr_length
    domains["distance"] = domains.apply(
        lambda d: interval_distance(start, end, d.domain_start, d.domain_end), axis=1
    )
    domains["is_edge"] = domains.apply(
        lambda d: is_domain_edge(start, end, d.domain_start, d.domain_end, d.overlap_length), axis=1
    )

    max_fraction = domains.lcr_domain_fraction.max()
    nearest = domains.sort_values(["distance", "i_evalue"], na_position="last").iloc[0]
    left = domains[domains.domain_end < start].sort_values("domain_end").tail(1)
    right = domains[domains.domain_start > end].sort_values("domain_start").head(1)

    labels = []
    if (domains.lcr_domain_fraction >= INTRINSIC_MIN_FRACTION).any():
        labels.append("domain_intrinsic")
    if domains.is_edge.any():
        labels.append("domain_edge")
    if max_fraction <= LINKER_MAX_FRACTION and not left.empty and not right.empty:
        labels.append("interdomain_linker")
    if max_fraction == 0 and nearest.distance <= ADJACENT_MAX_DISTANCE:
        labels.append("domain_adjacent")
    if max_fraction == 0 and nearest.distance > ADJACENT_MAX_DISTANCE:
        if end < domains.domain_start.min() or start > domains.domain_end.max():
            labels.append("distal_terminal")

    # One reporting class; class_labels retains all valid classes.
    priority = ["domain_edge", "domain_intrinsic", "interdomain_linker", "domain_adjacent", "distal_terminal"]
    primary = next((label for label in priority if label in labels), "unclassified_partial_overlap")

    return {
        "primary_class": primary,
        "class_labels": ";".join(labels) if labels else "unclassified_partial_overlap",
        "max_lcr_domain_fraction": round(max_fraction, 4),
        "nearest_domain_distance": int(nearest.distance),
        "nearest_pfam_accession": nearest.pfam_accession,
        "nearest_pfam_name": nearest.pfam_name,
        "flanking_left_pfam": left.iloc[0].pfam_accession if not left.empty else pd.NA,
        "flanking_right_pfam": right.iloc[0].pfam_accession if not right.empty else pd.NA,
    }


def main(lcr_file: Path, pfam_file: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    lcrs = pd.read_excel(lcr_file, sheet_name="all_results")
    domains = pd.read_excel(pfam_file, sheet_name="Pfam_hits")

    lcrs = lcrs.rename(columns={"start": "lcr_start", "end": "lcr_end", "length": "lcr_length"})
    lcrs["protein_accession"] = lcrs["protein_id"].map(accession)
    domains["protein_accession"] = domains["protein_id"].map(accession)
    domains = domains.drop_duplicates(["protein_accession", "pfam_accession", "domain_start", "domain_end"])

    by_protein = {key: group for key, group in domains.groupby("protein_accession", sort=False)}
    results = []
    relationships = []

    for lcr in lcrs.itertuples(index=False):
        protein_domains = by_protein.get(lcr.protein_accession, domains.iloc[0:0])
        result = classify_lcr(lcr, protein_domains)
        results.append(result)

        # Keep only true overlaps in the relationship table; adjacency is in the LCR table.
        for domain in protein_domains.itertuples(index=False):
            ov = overlap_length(lcr.lcr_start, lcr.lcr_end, domain.domain_start, domain.domain_end)
            if ov:
                relationships.append({
                    "protein_accession": lcr.protein_accession,
                    "method": lcr.method,
                    "lcr_start": lcr.lcr_start,
                    "lcr_end": lcr.lcr_end,
                    "pfam_accession": domain.pfam_accession,
                    "pfam_name": domain.pfam_name,
                    "domain_start": domain.domain_start,
                    "domain_end": domain.domain_end,
                    "overlap_length": ov,
                    "lcr_domain_fraction": round(ov / lcr.lcr_length, 4), # type: ignore
                    "domain_lcr_fraction": round(ov / (domain.domain_end - domain.domain_start + 1), 4), # type: ignore
                    "is_domain_edge": is_domain_edge(lcr.lcr_start, lcr.lcr_end, domain.domain_start, domain.domain_end, ov),
                })

    classified = pd.concat([lcrs.reset_index(drop=True), pd.DataFrame(results)], axis=1)
    relation_table = pd.DataFrame(relationships)
    summary = (classified.groupby(["method", "primary_class"], dropna=False)
               .size().reset_index(name="n_lcrs")
               .sort_values(["method", "n_lcrs"], ascending=[True, False]))

    classified.to_csv(output_dir / "lcr_position_classes.csv", index=False)
    relation_table.to_csv(output_dir / "lcr_domain_relationships.csv", index=False)
    summary.to_csv(output_dir / "lcr_position_class_summary.csv", index=False)
    print(f"Classified {len(classified):,} method-specific LCRs")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Position-based LCR classification relative to Pfam domains")
    parser.add_argument("--lcr-file", type=Path, default=Path(r"rbp_lcrs\lcr_methods_combined.xlsx"))
    parser.add_argument("--pfam-file", type=Path, default=Path(r"rbp_lcrs\pfam_lcr_overlap.xlsx"))
    parser.add_argument("--output-dir", type=Path, default=Path(r"lcr_analyses\domain_function"))
    args = parser.parse_args()
    main(args.lcr_file, args.pfam_file, args.output_dir)
