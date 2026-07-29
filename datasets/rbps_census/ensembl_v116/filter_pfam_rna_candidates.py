#!/usr/bin/env python3
"""Keep RNA-functional Pfam candidates in three evidence tiers.

This is a conservative name-based triage for the output of
find_new_pfam_rna_candidates.py. It deliberately excludes viral, prokaryotic
system, and DNA-centric domains. Inspect the review workbook before using its
contents as a final domain allow-list.
"""
import argparse
import re
import pandas as pd

# Order matters: exclusions are applied before positive rules.
EXCLUDE = re.compile(
    r"(?:flu_|flavi|rota|corona|\bcov|pico|pox|herpes|birna|bunya|arena|orbi|"
    r"seadorna|phage|\bt4\b|\bt7\b|baculo|mitovir|mycovirus|viral|virus|"
    r"capsid|matrix|hepatitis|ebola|mononeg|retro|rd?rp|rna[_-]?pol|v?rnap|"
    r"toxin|antitoxin|colicin|cdi|\bcas\d|crispr|restriction|\brnla\b|\bhica\b|"
    r"\bhig\b|\bsyme\b|\bbrnt\b|ghos|"
    r"\bdnab\b|primase|uvr|recq|recg|\brad\b|\brpa\b|dna[_-]?pol|helitron|"
    r"transpos|\brag1\b|chromodomain|\bchd\b|\bhth\b)", re.I)

DIRECT = re.compile(
    r"(?:\brhv\b|\bkh\b|rna.?bind|\brbd\b|\brrm\b|dsrm|dsrbd|"
    r"double.strand.*rna|\bpuf\b|\bpab\b|nab2|vts1|anticodon|"
    r"zf[-_]?ccch|zinc.*ccch|\bs1\b|\br3h\b|\bpu[a-z0-9_-]*\b)", re.I)

METABOLISM = re.compile(
    r"(?:rnase|rna[ _-]?se|ribonuc|exonuc|endonuc|\brnh\b|rnaseh|"
    r"helicase|\bdea[dhd]\b|\bdexq\b|ski2|pif1|sen1|\bxrn\b|dicer|drosha|"
    r"\btrm\b|trna|t.?rna|pseudour|methyltr|methyl.*rna|rna.*methyl|"
    r"rna.*cap|decap|polyapol|rna.?lig|rli[g]?|rnas[e]?j|ribonucleotide)", re.I)

RNP = re.compile(
    r"(?:ribosom|\bmrp[ls]\b|\brp[ols]\b|\brrn\b|snrnp|\bsnu\b|sf3|\bprp\b|"
    r"splice|\blsm\b|\bnop\b|\butp\b|\brra\b|\brpf\b|\brpp\b|exosc|mtr4|"
    r"brr2|cpsf|pcf11|clp1|mago|\bedc\b|larp|mterf|\bsrp\b|translation|"
    r"\beif\b|\bef[gt]\b|aminoacyl|synthetase|\bpeptidyl\b)", re.I)

def classify(name: str):
    """Return tier and a short reproducible rationale, or (None, None)."""
    name = str(name)
    if EXCLUDE.search(name):
        return None, None
    if DIRECT.search(name):
        return "direct_RNA_binding", "RNA-binding fold or RNA-binding module"
    if METABOLISM.search(name):
        return "RNA_metabolism", "RNA processing, modification, turnover, or remodelling"
    if RNP.search(name):
        return "RNP_component_RNA_pathway", "stable RNP, ribosome, splicing, or RNA-biogenesis component"
    return None, None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_xlsx", help="pfam38.2_new_rna_candidates.xlsx")
    ap.add_argument("-o", "--output", default="pfam38.2_rna_functional_tiers.xlsx")
    ap.add_argument("--include-go", action="store_true", help="also classify rows sourced from GO")
    args = ap.parse_args()

    df = pd.read_excel(args.input_xlsx)
    required = {"accession", "name", "source"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns: {', '.join(sorted(missing))}")

    work = df.copy() if args.include_go else df[df["source"].astype(str).str.lower().eq("metadata")].copy()
    assigned = work["name"].map(classify)
    work[["evidence_tier", "triage_rationale"]] = pd.DataFrame(assigned.tolist(), index=work.index)
    kept = work.dropna(subset=["evidence_tier"]).sort_values(["evidence_tier", "accession"])

    # Audit sheets make discarded candidates and unassigned candidates transparent.
    excluded = work[work["evidence_tier"].isna() & work["name"].map(lambda x: bool(EXCLUDE.search(str(x))))].copy()
    excluded["reason"] = "viral, prokaryotic defence/toxin, or DNA-centric"
    review = work[work["evidence_tier"].isna() & ~work.index.isin(excluded.index)].copy()
    review["reason"] = "ambiguous name; manual functional review required"

    with pd.ExcelWriter(args.output, engine="openpyxl") as xw:
        kept.to_excel(xw, sheet_name="retained_tiers", index=False)
        kept[kept.evidence_tier.eq("direct_RNA_binding")].to_excel(xw, sheet_name="direct_binding", index=False)
        kept[kept.evidence_tier.eq("RNA_metabolism")].to_excel(xw, sheet_name="RNA_metabolism", index=False)
        kept[kept.evidence_tier.eq("RNP_component_RNA_pathway")].to_excel(xw, sheet_name="RNP_component", index=False)
        excluded.to_excel(xw, sheet_name="excluded", index=False)
        review.to_excel(xw, sheet_name="manual_review", index=False)
    print(f"Retained {len(kept)} of {len(work)} examined rows: " + ", ".join(f"{k}={v}" for k,v in kept.evidence_tier.value_counts().items()))

if __name__ == "__main__":
    main()
