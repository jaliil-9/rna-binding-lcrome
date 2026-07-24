#!/usr/bin/env python3
"""Build one gene-level table from the 2014 RBP census and new Pfam-38 candidates."""

import sys
import time
import pandas as pd
import requests

# Usage:
# python build_rbp_sheet.py 2014_census.xls new_candidates.xlsx output.xlsx

OLD_SHEET = "RBP table"
NEW_SHEET = "symbol_deduplicated_859"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FIELDS = [
    "accession", "id", "gene_names", "protein_name", "reviewed", "length", "sequence",
    "go_id", "go_p", "go_f", "cc_subcellular_location", "cc_domain",
]


def clean_gene(x):
    return str(x).strip().upper() if pd.notna(x) else ""


def join_unique(values):
    values = [str(v).strip() for v in values if pd.notna(v) and str(v).strip()]
    return "; ".join(dict.fromkeys(values))


def read_inputs(old_file, new_file):
    old = pd.read_excel(old_file, sheet_name=OLD_SHEET)
    old = old.rename(columns={
        "gene name": "gene_symbol", "description": "description_2014",
        "gene id": "ensembl_gene_id_2014", "protein id": "ensembl_protein_id_2014",
        "number of domains": "domain_count_2014", "domains[count]": "domains_2014",
    })
    old["gene_symbol"] = old["gene_symbol"].map(clean_gene)
    old = old[["gene_symbol", "description_2014", "ensembl_gene_id_2014",
               "ensembl_protein_id_2014", "domain_count_2014", "domains_2014"]]
    old["in_2014_census"] = "yes"

    new = pd.read_excel(new_file, sheet_name=NEW_SHEET)
    new = new.rename(columns={
        "gene_symbol": "gene_symbol", "ensembl_gene_id": "ensembl_gene_id_new",
        "human_sequence_id": "ensembl_protein_id_new", "protein_length": "protein_length_new",
        "fasta_description": "description_new", "number_of_domain_instances": "domain_count_new",
        "pfam_accessions": "pfam_accessions_new", "pfam_models": "pfam_domains_new",
        "full_sequence": "sequence_ensembl",
    })
    new["gene_symbol"] = new["gene_symbol"].map(clean_gene)
    wanted = ["gene_symbol", "description_new", "ensembl_gene_id_new", "ensembl_protein_id_new",
              "protein_length_new", "sequence_ensembl", "domain_count_new",
              "pfam_accessions_new", "pfam_domains_new"]
    new = new[[c for c in wanted if c in new.columns]]
    new["in_new_pfam38_set"] = "yes"
    return old, new


def uniprot_lookup(genes):
    """Retrieve one preferred human UniProtKB entry per gene symbol."""
    rows = []
    for start in range(0, len(genes), 100):
        batch = genes[start:start + 100]
        gene_query = " OR ".join(f"gene_exact:{g}" for g in batch)
        params = {"query": f"organism_id:9606 AND ({gene_query})", "format": "tsv",
                  "fields": ",".join(UNIPROT_FIELDS), "size": 500}
        response = requests.get(UNIPROT_URL, params=params, timeout=90)
        response.raise_for_status()
        rows.append(pd.read_csv(pd.io.common.StringIO(response.text), sep="\t"))
        time.sleep(0.2)
    if not rows:
        return pd.DataFrame()
    hits = pd.concat(rows, ignore_index=True)
    if hits.empty:
        return hits

    # Match each entry to every submitted symbol listed in UniProt's Gene Names field.
    hits["gene_symbol"] = hits["Gene Names"].fillna("").str.upper().str.split()
    hits = hits.explode("gene_symbol")
    hits = hits[hits["gene_symbol"].isin(genes)].copy()
    hits["reviewed_rank"] = hits["Reviewed"].eq("reviewed").astype(int)
    hits = hits.sort_values(["gene_symbol", "reviewed_rank", "Length"], ascending=[True, False, False])
    hits = hits.drop_duplicates("gene_symbol", keep="first")
    return hits.rename(columns={
        "Entry": "uniprot_accession", "Entry Name": "uniprot_entry_name",
        "Protein names": "uniprot_protein_name", "Reviewed": "uniprot_reviewed",
        "Length": "uniprot_length", "Sequence": "uniprot_sequence",
        "Gene Ontology IDs": "go_ids", "Gene Ontology (biological process)": "go_biological_process",
        "Gene Ontology (molecular function)": "go_molecular_function",
        "Subcellular location [CC]": "subcellular_location", "Domain [CC]": "uniprot_domains",
    }).drop(columns=["Gene Names", "reviewed_rank"], errors="ignore")


def main(old_file, new_file, output_file):
    old, new = read_inputs(old_file, new_file)
    table = old.merge(new, on="gene_symbol", how="outer")
    table["in_2014_census"] = table["in_2014_census"].fillna("no")
    table["in_new_pfam38_set"] = table["in_new_pfam38_set"].fillna("no")
    table["source"] = table.apply(lambda r: "2014; Pfam-38" if r.in_2014_census == "yes" and r.in_new_pfam38_set == "yes" else ("2014" if r.in_2014_census == "yes" else "Pfam-38"), axis=1)

    genes = table["gene_symbol"].dropna().unique().tolist()
    uni = uniprot_lookup(genes)
    table = table.merge(uni, on="gene_symbol", how="left")

    # Prefer current UniProt sequence/length; retain Ensembl data when UniProt has no match.
    table["protein_length"] = table["uniprot_length"].fillna(table.get("protein_length_new"))
    table["protein_sequence"] = table["uniprot_sequence"].fillna(table.get("sequence_ensembl"))
    domain_counts = table[["domain_count_2014", "domain_count_new"]].apply(
        pd.to_numeric, errors="coerce"
    )
    table["rna_domain_count"] = domain_counts.fillna(0).max(axis=1).astype(int)    
    table["rna_domains"] = table.apply(lambda r: join_unique([r.get("domains_2014"), r.get("pfam_domains_new")]), axis=1)
    columns = [
        "gene_symbol", "source", "in_2014_census", "in_new_pfam38_set",
        "uniprot_accession", "uniprot_entry_name", "uniprot_reviewed", "uniprot_protein_name",
        "ensembl_gene_id_2014", "ensembl_protein_id_2014", "ensembl_gene_id_new", "ensembl_protein_id_new",
        "protein_length", "protein_sequence", "go_ids", "go_biological_process", "go_molecular_function",
        "subcellular_location", "rna_domain_count", "rna_domains", "uniprot_domains",
        "domains_2014", "pfam_accessions_new", "pfam_domains_new",
        "description_2014", "description_new",
    ]
    columns = [c for c in columns if c in table.columns]
    table = table[columns].sort_values("gene_symbol")
    table.to_excel(output_file, index=False, sheet_name="combined_RBPs")
    print(f"Wrote {len(table):,} genes to {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit("Usage: python build_rbp_sheet.py 2014_census.xls new_candidates.xlsx output.xlsx")
    main(*sys.argv[1:])
