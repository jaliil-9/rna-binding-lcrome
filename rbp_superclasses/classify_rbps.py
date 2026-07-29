#!/usr/bin/env python3
"""Assign Gerstberger-style RNA target superclasses to a UniProt RBP workbook.

Usage
-----
python classify_rbps.py \
  --input rbps_input.xlsx \
  --gerstberger gerstberger_s3.xlsx \
  --output rbp_rna_classification.xlsx

Inputs
------
1) UniProt workbook (first sheet, or --input-sheet) with columns:
   Entry; Entry Name; Protein names; Gene Ontology (biological process);
   Gene Ontology (molecular function); Subcellular location [CC]; Domain [CC]
2) Gerstberger Supplementary Table S3 workbook/CSV with columns:
   gene name; protein id; consensus RNA target; putative RNA target;
   supporting evidence (# pubmed ID)

Output
------
An XLSX file containing all_combined plus one sheet per primary RNA class.
Only these six columns are written: uniprot_accession, rna_primary_class,
rna_secondary_class, confidence, evidence_source, evidence_matched.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

CLASSES = [
    "ribosomal protein", "mRNA", "tRNA", "pre-rRNA", "snRNA",
    "snoRNA", "ncRNA", "diverse", "unknown",
]

# Target-specific terms only. Generic RBDs (RRM, KH, helicase, zinc finger,
# G-patch, etc.) intentionally do not assign a target class.
TERMS = {
    "ribosomal protein": [
        r"\bribosomal protein\b", r"\bstructural constituent of ribosome\b",
        r"\bmitochondrial ribosomal\b", r"\bcytosolic ribosomal\b",
        r"\bsmall ribosomal subunit\b", r"\blarge ribosomal subunit\b",
    ],
    "mRNA": [
        r"\bmessenger rna\b", r"\bmRNA\b", r"\bpre-mRNA\b",
        r"\bmrna splicing\b", r"\brna splicing\b", r"\bspliceosomal\b",
        r"\bspliceosome\b", r"\bhnrnp\b", r"\bpolyadenylation\b",
        r"\bpoly\(a\)\b", r"\bmrna export\b", r"\bmrna localization\b",
        r"\bmrna decay\b", r"\bmrna stability\b", r"\bexon junction\b",
        r"\bcap-binding\b", r"\btranslation initiation\b",
        r"\btranslation regulation\b",
    ],
    "tRNA": [
        r"\btRNA\b", r"\btransfer rna\b", r"\baminoacyl-tRNA\b",
        r"\baminoacyl trna\b", r"\banticodon\b", r"\btRNA splicing\b",
        r"\btRNA processing\b", r"\btRNA modification\b",
        r"\btRNA methyltransferase\b",
    ],
    "pre-rRNA": [
        r"\bpre-rRNA\b", r"\bpre rRNA\b", r"\brRNA processing\b",
        r"\brRNA maturation\b", r"\brRNA modification\b", r"\b18S rRNA\b",
        r"\b28S rRNA\b", r"\b5\.8S rRNA\b", r"\b5S rRNA\b",
        r"\bribosome biogenesis\b", r"\bsmall subunit processome\b",
        r"\bnucleolar rna processing\b",
    ],
    "snRNA": [
        r"\bsnRNA\b", r"\bsmall nuclear rna\b", r"\bsnRNP\b",
        r"\bsmall nuclear ribonucleoprotein\b", r"\bU1 snRNA\b",
        r"\bU2 snRNA\b", r"\bU4 snRNA\b", r"\bU5 snRNA\b",
        r"\bU6 snRNA\b", r"\bSMN complex\b",
    ],
    "snoRNA": [
        r"\bsnoRNA\b", r"\bsmall nucleolar rna\b", r"\bsnoRNP\b",
        r"\bsmall nucleolar ribonucleoprotein\b", r"\bscaRNA\b",
        r"\bsmall cajal body-specific rna\b", r"\bbox c/d\b", r"\bbox h/aca\b",
    ],
    "ncRNA": [
        r"\bmiRNA\b", r"\bmicroRNA\b", r"\bpiRNA\b", r"\bPIWI\b",
        r"\blncRNA\b", r"\blong non-coding rna\b", r"\blong noncoding rna\b",
        r"\b7SK\b", r"\b7SL\b", r"\bY RNA\b", r"\bvault rna\b",
        r"\btelomerase rna\b", r"\bRNase P\b", r"\bRNase MRP\b",
        r"\bmicroprocessor\b", r"\bargonaute\b",
    ],
    "diverse": [
        r"\brna exosome\b", r"\bribonuclease\b", r"\brna exonuclease\b",
        r"\brna endonuclease\b", r"\bgeneral rna turnover\b",
        r"\brna surveillance\b", r"\brna degradation\b", r"\brna:DNA hybrid\b",
        r"\brna-dna hybrid\b",
    ],
}

UNIPROT_COLUMNS = {
    "accession": "uniprot_accession",
    "entry_name": "uniprot_entry_name",
    "protein_name": "uniprot_protein_name",
    "go_bp": "go_biological_process",
    "go_mf": "go_molecular_function",
    "location": "subcellular_location",
    "domain": "uniprot_domains",
}


def text(value) -> str:
    return "" if pd.isna(value) else str(value).strip()


def accession(value) -> str:
    """Normalise a UniProt accession and remove an optional isoform suffix."""
    return text(value).upper().split("-")[0]


def gene_from_entry_name(value) -> str:
    """TSR3_HUMAN -> TSR3; used only as a pragmatic S3 fallback key."""
    return text(value).upper().split("_")[0]


def normalise_s3_class(value) -> str:
    """Map expected S3 labels to one controlled project vocabulary."""
    x = text(value).lower().replace("–", "-")
    x = re.sub(r"\s+", " ", x)
    aliases = {
        "ribosome": "ribosomal protein",
        "ribosomal proteins": "ribosomal protein",
        "ribosomal protein": "ribosomal protein",
        "mrna": "mRNA", "mrna-binding": "mRNA", "mrna-binding proteins": "mRNA",
        "trna": "tRNA", "trna-binding": "tRNA", "trna-binding proteins": "tRNA",
        "pre-rrna": "pre-rRNA", "pre-rrna-binding": "pre-rRNA",
        "rrna": "pre-rRNA", "rrna-binding": "pre-rRNA",
        "snrna": "snRNA", "snrna-binding": "snRNA",
        "snorna": "snoRNA", "snorna-binding": "snoRNA",
        "ncrna": "ncRNA", "ncrna-binding": "ncRNA",
        "diverse": "diverse", "diverse targets": "diverse",
        "unknown": "unknown", "unknown targets": "unknown",
    }
    return aliases.get(x, x)


def read_table(path: str, sheet: str | int = 0) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() == ".csv":
        return pd.read_csv(p)
    return pd.read_excel(p, sheet_name=sheet)


def require_columns(df: pd.DataFrame, columns: List[str], label: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing columns: {missing}\nFound: {df.columns.tolist()}")


def load_s3(path: str, sheet: str | int = 0) -> Tuple[Dict[str, dict], Dict[str, dict]]:
    """Load Gerstberger S3, automatically locating its sheet and header row."""
    required = [
        "gene name",
        "protein id",
        "consensus rna target",
        "putative rna target",
        "supporting evidence (# pubmed id)",
    ]
    p = Path(path)
    s3 = pd.read_excel(
            p,
            sheet_name="RBP table"
        )

    # Standardize headers after the correct header row is found.
    s3.columns = [text(c).strip().lower() for c in s3.columns]
    require_columns(s3, required, "Gerstberger S3")

    protein_map = {}
    gene_map = {}

    for _, r in s3.iterrows():
        primary = normalise_s3_class(r["consensus rna target"])

        raw_primary = text(r["consensus rna target"])
        primary = normalise_s3_class(raw_primary)

        if primary not in CLASSES:
            continue

        secondary = normalise_s3_class(r["putative rna target"])

        if secondary == primary or secondary not in CLASSES:
            secondary = ""

        pubmed = text(r["supporting evidence (# pubmed id)"])

        record = {
            "rna_primary_class": primary,
            "rna_secondary_class": secondary,
            "evidence_matched": (
                f"Gerstberger S3 consensus RNA target: {primary}; "
                f"putative RNA target: {secondary or 'none'}; "
                f"PubMed: {pubmed or 'not listed'}"
            ),
        }

        pid = text(r["protein id"]).upper()
        gene = text(r["gene name"]).upper()

        if pid:
            protein_map[pid] = record

        if gene:
            gene_map[gene] = record

    return protein_map, gene_map


def matched_terms(value: str) -> Dict[str, List[str]]:
    hits = {}
    for cls, patterns in TERMS.items():
        found = [p for p in patterns if re.search(p, value, flags=re.I)]
        if found:
            hits[cls] = found
    return hits


def printable(pattern: str) -> str:
    return pattern.replace(r"\b", "").replace(r"\.", ".").replace(r"\(", "(").replace(r"\)", ")")


def classify_row(row: pd.Series, s3_protein: Dict[str, dict], s3_gene: Dict[str, dict]) -> dict:
    acc = accession(row.get(UNIPROT_COLUMNS["accession"], ""))

    # Direct S3 match. protein id is checked first; Entry Name-derived gene is fallback.
    possible_ids = [acc]
    if "ensembl_protein_id" in row.index:
        possible_ids.append(text(row["ensembl_protein_id"]).upper())
    for pid in possible_ids:
        if pid and pid in s3_protein:
            h = s3_protein[pid]
            return {"uniprot_accession": acc, **h, "confidence": "very high",
                    "evidence_source": "literature based (Gerstberger S3)"}

    gene = gene_from_entry_name(row.get(UNIPROT_COLUMNS["entry_name"], ""))
    if gene and gene in s3_gene:
        h = s3_gene[gene]
        return {"uniprot_accession": acc, **h, "confidence": "very high",
                "evidence_source": "literature based (Gerstberger S3)"}

    # Collect independent support by field from the UniProt workbook.
    support: Dict[str, List[Tuple[str, str, str]]] = {c: [] for c in CLASSES[:-1]}
    for field in ("protein_name", "go_bp", "go_mf", "location", "domain"):
        col = UNIPROT_COLUMNS[field]
        if col not in row.index:
            continue
        for cls, patterns in matched_terms(text(row[col])).items():
            for pat in patterns:
                support[cls].append((field, col, printable(pat)))
    support = {cls: items for cls, items in support.items() if items}

    if not support:
        return {"uniprot_accession": acc, "rna_primary_class": "unknown",
                "rna_secondary_class": "", "confidence": "unknown",
                "evidence_source": "", "evidence_matched": ""}

    # Class ranking = number of distinct UniProt fields, then number of matched terms.
    ranked = sorted(
        support,
        key=lambda c: (-len({x[0] for x in support[c]}), -len(support[c]), c),
    )
    primary = ranked[0]
    secondary = "; ".join(ranked[1:])
    hits = support[primary]
    fields = {x[0] for x in hits}

    # Requested, intentionally simple confidence scheme.
    if len(fields) >= 2:
        confidence = "high"
    elif fields == {"protein_name"}:
        # A target-specific name (e.g. tRNA methyltransferase) is high; otherwise low.
        target_words = ("trna", "rrna", "mrna", "snrna", "snorna", "mirna", "pirna", "ribosomal", "rna exosome")
        confidence = "high" if any(w in text(row.get(UNIPROT_COLUMNS["protein_name"], "")).lower() for w in target_words) else "low"
    elif fields.issubset({"go_bp", "go_mf", "location", "domain"}):
        confidence = "medium"
    else:
        confidence = "unknown"

    sources = []
    if {"go_bp", "go_mf"} & fields:
        sources.append("GO based")
    if {"protein_name", "location", "domain"} & fields:
        sources.append("UniProt annotation/metadata")

    return {
        "uniprot_accession": acc,
        "rna_primary_class": primary,
        "rna_secondary_class": secondary,
        "confidence": confidence,
        "evidence_source": "; ".join(sources),
        "evidence_matched": "; ".join(f"{col}: {term}" for _, col, term in hits),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Simple Gerstberger-style RBP RNA-target classifier")
    parser.add_argument("--input", required=True, help="UniProt-mapped input XLSX/CSV")
    parser.add_argument("--gerstberger", required=True, help="Original Gerstberger S3 XLSX/CSV")
    parser.add_argument("--output", default="rbp_rna_classification.xlsx", help="Output XLSX path")
    parser.add_argument("--input-sheet", default=0, help="Input sheet name or zero-based index")
    parser.add_argument("--s3-sheet", default=0, help="S3 sheet name or zero-based index")
    args = parser.parse_args()

    def sheet_arg(x):
        return int(x) if str(x).isdigit() else x

    inp = read_table(args.input, sheet_arg(args.input_sheet)).copy()
    require_columns(inp, [UNIPROT_COLUMNS["accession"]], "Input workbook")
    s3_protein, s3_gene = load_s3(args.gerstberger, sheet_arg(args.s3_sheet))

    result = pd.DataFrame([classify_row(r, s3_protein, s3_gene) for _, r in inp.iterrows()])
    out_cols = ["uniprot_accession", "rna_primary_class", "rna_secondary_class", "confidence", "evidence_source", "evidence_matched"]
    result = result[out_cols]

    with pd.ExcelWriter(args.output, engine="openpyxl") as writer:
        result.to_excel(writer, sheet_name="all_combined", index=False)
        for cls in CLASSES:
            result[result["rna_primary_class"] == cls].to_excel(
                writer, sheet_name=cls[:31], index=False
            )

    print(f"Wrote {len(result)} rows to {Path(args.output).resolve()}")
    print("Primary classes:\n", result["rna_primary_class"].value_counts())
    print("Confidence:\n", result["confidence"].value_counts())


if __name__ == "__main__":
    main()
