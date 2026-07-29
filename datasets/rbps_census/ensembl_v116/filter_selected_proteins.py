#!/usr/bin/env python3
import sys
from pathlib import Path
import pandas as pd

DOMTBLOUT = Path(sys.argv[1])

COLS = [
    "target_name", "target_accession", "target_len",
    "query_name", "query_accession", "query_len",
    "full_evalue", "full_score", "full_bias",
    "domain_number", "domain_of",
    "c_evalue", "i_evalue", "domain_score", "domain_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to",
    "env_from", "env_to", "acc", "description"
]

def read_domtblout(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split(maxsplit=22)
            if len(fields) == 23:
                rows.append(fields)

    df = pd.DataFrame(rows, columns=COLS)

    numeric = [
        "target_len", "query_len", "full_evalue", "full_score", "full_bias",
        "domain_number", "domain_of", "c_evalue", "i_evalue",
        "domain_score", "domain_bias", "hmm_from", "hmm_to",
        "ali_from", "ali_to", "env_from", "env_to", "acc"
    ]
    for col in numeric:
        df[col] = pd.to_numeric(df[col])

    df["hmm_coverage"] = (
        (df["hmm_to"] - df["hmm_from"] + 1) / df["target_len"]
    )

    return df

df = read_domtblout(DOMTBLOUT)

selected = pd.read_excel(
    "pfam38_census_comparison.xlsx",
    sheet_name="2_Deduplicated"
)

selected_ids = set(selected["human_sequence_id"])
df = df[df["query_name"].isin(selected_ids)].copy()

# Strict, final-set domain criteria
strict = df[
    (df["i_evalue"] <= 1e-5) &
    (df["hmm_coverage"] >= 0.50)
    ].copy()

# Useful but not final: retain for focused later inspection
borderline = df[
    (df["i_evalue"] <= 1e-3) &
    (df["hmm_coverage"] >= 0.30)
].copy()

borderline = borderline[
    ~borderline.index.isin(strict.index)
].copy()

def summarise(hits):
    if hits.empty:
        return pd.DataFrame()

    return (
        hits.sort_values(["query_name", "i_evalue", "domain_score"])
        .groupby("query_name", as_index=False)
        .agg(
            n_domains=("target_accession", "size"),
            pfam_accessions=("target_accession", lambda x: "; ".join(sorted(set(x)))),
            pfam_models=("target_name", lambda x: "; ".join(sorted(set(x)))),
            best_i_evalue=("i_evalue", "min"),
            best_domain_score=("domain_score", "max"),
            best_hmm_coverage=("hmm_coverage", "max"),
        )
    )


print(f"All domain hits:        {len(df)}")
print(f"Strict domain hits:     {len(strict)}")
print(f"Strict proteins:        {strict['query_name'].nunique()}")
print(f"Borderline proteins:    {borderline['query_name'].nunique()}")

strict_summary = summarise(strict)
borderline_summary = summarise(borderline)

# Domain-level TSVs
strict.to_csv("strict_domain_hits.tsv", sep="\t", index=False)
borderline.to_csv("borderline_domain_hits.tsv", sep="\t", index=False)

# Protein-level TSVs
strict_summary.to_csv(
    "strict_rna_domain_candidates.tsv", sep="\t", index=False
)
borderline_summary.to_csv(
    "borderline_rna_domain_candidates.tsv", sep="\t", index=False
)

# Protein-level Excel workbook: one sheet per confidence class
with pd.ExcelWriter(
    "rna_domain_protein_classes_v2.xlsx",
    engine="openpyxl"
) as writer:
    strict_summary.to_excel(writer, sheet_name="Strict", index=False)
    borderline_summary.to_excel(writer, sheet_name="Borderline", index=False)