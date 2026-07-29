#!/usr/bin/env python3
"""Combine the 2014 RBP census with Ensembl v116 Pfam-38 strict/borderline candidates.

Usage:
  python build_rbp_sheet.py rna_binding_proteins.xls ensembl116_pfam38_filtered_proteins.xlsx combined_rbps.xlsx
"""
import sys
import time
from io import StringIO
import pandas as pd
import requests

OLD_SHEET = "RBP table"
CANDIDATE_SHEETS = ("Strict", "Borderline")
ENSEMBL_LOOKUP = "https://rest.ensembl.org/lookup/id"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FIELDS = [
    "accession", "id", "gene_names", "protein_name", "reviewed", "length", "sequence",
    "go_id", "go_p", "go_f", "cc_subcellular_location", "cc_domain",
]

def clean_text(x):
    return "" if pd.isna(x) else str(x).strip()

def clean_gene(x):
    return clean_text(x).upper()

def strip_version(x):
    return clean_text(x).split(".")[0]

def join_unique(values):
    seen, out = set(), []
    for value in values:
        value = clean_text(value)
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return "; ".join(out)

def get_col(df, name, default=""):
    return df[name] if name in df.columns else pd.Series(default, index=df.index)

def read_census(path):
    old = pd.read_excel(path, sheet_name=OLD_SHEET)
    old = old.rename(columns={
        "gene name": "gene_symbol", "description": "description_2014",
        "gene id": "ensembl_gene_id_2014", "protein id": "ensembl_protein_id_2014",
        "number of domains": "domain_count_2014", "domains[count]": "domains_2014",
    })
    required = ["gene_symbol", "ensembl_protein_id_2014"]
    missing = [c for c in required if c not in old.columns]
    if missing:
        raise ValueError(f"Census sheet is missing columns: {missing}")
    old["gene_symbol"] = old["gene_symbol"].map(clean_gene)
    old["ensp_key"] = old["ensembl_protein_id_2014"].map(strip_version)
    keep = ["gene_symbol", "ensp_key", "description_2014", "ensembl_gene_id_2014",
            "ensembl_protein_id_2014", "domain_count_2014", "domains_2014"]
    old = old[[c for c in keep if c in old.columns]].copy()
    old["in_2014_census"] = "yes"
    return old

def read_candidates(path):
    frames = []
    for sheet in CANDIDATE_SHEETS:
        df = pd.read_excel(path, sheet_name=sheet)
        required = {"query_name", "n_domains", "pfam_accessions", "pfam_models"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"{sheet} sheet is missing columns: {sorted(missing)}")
        df["candidate_class"] = sheet
        frames.append(df)
    df = pd.concat(frames, ignore_index=True).rename(columns={
        "query_name": "ensembl_protein_id_new", "n_domains": "domain_count_new",
        "pfam_accessions": "pfam_accessions_new", "pfam_models": "pfam_domains_new",
    })
    df["ensp_key"] = df["ensembl_protein_id_new"].map(strip_version)
    aggregations = {
        "ensembl_protein_id_new": "first", "domain_count_new": "max",
        "pfam_accessions_new": join_unique, "pfam_domains_new": join_unique,
        "candidate_class": join_unique,
    }
    for col, rule in [("best_i_evalue", "min"), ("best_domain_score", "max"),
                      ("best_hmm_coverage", "max"), ("best_bias_fraction", "min")]:
        if col in df.columns:
            aggregations[col] = rule
    df = df.groupby("ensp_key", as_index=False).agg(aggregations)
    df["in_new_pfam38_set"] = "yes"
    return df

def ensembl_metadata(ensp_keys, batch_size=50, max_retries=5):
    """Resolve ENSP IDs to Ensembl metadata with retry/backoff."""
    rows = []
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
    }
    session = requests.Session()

    for start in range(0, len(ensp_keys), batch_size):
        ids = ensp_keys[start:start + batch_size]

        for attempt in range(max_retries):
            try:
                response = session.post(
                    ENSEMBL_LOOKUP,
                    json={"ids": ids},
                    headers=headers,
                    timeout=(15, 180),
                )

                if response.status_code == 429:
                    wait = float(response.headers.get(
                        "Retry-After", 2 ** (attempt + 1)
                    ))
                    print(f"Ensembl rate limit; waiting {wait:.0f} s")
                    time.sleep(wait)
                    continue

                response.raise_for_status()
                payload = response.json()
                break

            except (requests.exceptions.Timeout,
                    requests.exceptions.ConnectionError) as exc:
                if attempt == max_retries - 1:
                    raise RuntimeError(
                        f"Ensembl lookup failed for IDs {start}–"
                        f"{start + len(ids) - 1} after {max_retries} attempts"
                    ) from exc

                wait = 2 ** (attempt + 1)
                print(
                    f"Ensembl timeout for batch {start}–"
                    f"{start + len(ids) - 1}; retrying in {wait} s"
                )
                time.sleep(wait)

        else:
            raise RuntimeError(f"Could not retrieve Ensembl batch starting at {start}")

        for ensp, record in payload.items():
            if not record:
                continue

            rows.append({
                "ensp_key": ensp,
                "ensembl_gene_id_new": record.get("Parent", ""),
                "gene_symbol_ensembl": record.get("display_name", ""),
                "protein_length_new": record.get("length"),
                "description_new": record.get("description", ""),
                "sequence_ensembl": record.get("seq", ""),
            })

        print(f"Ensembl: resolved {min(start + len(ids), len(ensp_keys))}/{len(ensp_keys)}")
        time.sleep(0.25)

    return pd.DataFrame(rows)

def uniprot_metadata(genes, batch_size=100):
    """Retrieve one preferred human UniProtKB entry per HGNC symbol."""
    frames, session = [], requests.Session()
    for start in range(0, len(genes), batch_size):
        batch = genes[start:start + batch_size]
        query = " OR ".join(f"gene_exact:{g}" for g in batch)
        params = {"query": f"organism_id:9606 AND ({query})", "format": "tsv",
                  "fields": ",".join(UNIPROT_FIELDS), "size": 500}
        response = session.get(UNIPROT_SEARCH, params=params, timeout=90)
        response.raise_for_status()
        frames.append(pd.read_csv(StringIO(response.text), sep="\t"))
        time.sleep(0.2)
    if not frames:
        return pd.DataFrame(columns=["gene_symbol"])
    hits = pd.concat(frames, ignore_index=True)
    if hits.empty:
        return pd.DataFrame(columns=["gene_symbol"])
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

def main(census_file, candidate_file, output_file):
    old = read_census(census_file)
    new = read_candidates(candidate_file)
    meta = ensembl_metadata(new["ensp_key"].tolist())
    new = new.merge(meta, on="ensp_key", how="left")

    # Do NOT discard candidate rows if Ensembl annotation is missing.
    new["gene_symbol"] = new["gene_symbol_ensembl"].fillna("").map(clean_gene)

    unmapped = new["gene_symbol"].eq("")
    print(f"Pfam candidates with unresolved Ensembl symbol: {unmapped.sum():,}")
    print(f"Pfam candidates resolved to a symbol: {(~unmapped).sum():,}")

    # Combine census and Pfam candidate rows; retain candidates even if
    # Ensembl metadata could not resolve an HGNC/gene symbol.
    table = pd.concat([old, new], ignore_index=True, sort=False)

    table["gene_symbol"] = get_col(table, "gene_symbol", "").map(clean_gene)
    table["ensp_key"] = get_col(table, "ensp_key", "").fillna("").astype(str)

    table["in_2014_census"] = (
        get_col(table, "in_2014_census", "no").replace("", "no").fillna("no")
    )
    table["in_new_pfam38_set"] = (
        get_col(table, "in_new_pfam38_set", "no").replace("", "no").fillna("no")
    )

    # Preserve unresolved new candidates with an ENSP-based placeholder.
    missing_symbol = table["gene_symbol"].eq("")
    table.loc[missing_symbol, "gene_symbol"] = (
        "UNMAPPED_" + table.loc[missing_symbol, "ensp_key"].replace("", "UNKNOWN")
    )
    gene_rows = []
    for gene, g in table.groupby("gene_symbol", sort=True):
        row = {"gene_symbol": gene}
        row["in_2014_census"] = "yes" if (g["in_2014_census"] == "yes").any() else "no"
        row["in_new_pfam38_set"] = "yes" if (g["in_new_pfam38_set"] == "yes").any() else "no"
        for col in ["ensembl_gene_id_2014", "ensembl_protein_id_2014", "description_2014",
                    "ensembl_gene_id_new", "ensembl_protein_id_new", "description_new",
                    "candidate_class", "pfam_accessions_new", "pfam_domains_new"]:
            if col in g: row[col] = join_unique(g[col])
        for col, rule in [("domain_count_2014", "max"), ("domain_count_new", "max"),
                          ("best_i_evalue", "min"), ("best_domain_score", "max"),
                          ("best_hmm_coverage", "max"), ("best_bias_fraction", "min"),
                          ("protein_length_new", "max")]:
            if col in g: row[col] = pd.to_numeric(g[col], errors="coerce").agg(rule)
        row["sequence_ensembl"] = join_unique(g.get("sequence_ensembl", []))
        gene_rows.append(row)
    result = pd.DataFrame(gene_rows)
    result["source"] = result.apply(lambda r: "2014; Pfam-38" if r.in_2014_census == "yes" and r.in_new_pfam38_set == "yes" else ("2014" if r.in_2014_census == "yes" else "Pfam-38"), axis=1)

    uni = uniprot_metadata(result["gene_symbol"].tolist())
    result = result.merge(uni, on="gene_symbol", how="left")
    result["protein_length"] = result["uniprot_length"].fillna(result.get("protein_length_new"))
    result["protein_sequence"] = result["uniprot_sequence"].fillna(result.get("sequence_ensembl"))
    result["rna_domain_count"] = result[[c for c in ["domain_count_2014", "domain_count_new"] if c in result]].fillna(0).max(axis=1).astype(int)
    result["rna_domains"] = result.apply(lambda r: join_unique([r.get("domains_2014", ""), r.get("pfam_domains_new", "")]), axis=1)

    first = ["gene_symbol", "source", "in_2014_census", "in_new_pfam38_set", "candidate_class",
             "uniprot_accession", "uniprot_entry_name", "uniprot_reviewed", "uniprot_protein_name",
             "ensembl_gene_id_2014", "ensembl_protein_id_2014", "ensembl_gene_id_new", "ensembl_protein_id_new",
             "protein_length", "protein_sequence", "go_ids", "go_biological_process", "go_molecular_function",
             "subcellular_location", "rna_domain_count", "rna_domains", "pfam_accessions_new", "pfam_domains_new",
             "best_i_evalue", "best_domain_score", "best_hmm_coverage", "best_bias_fraction", "uniprot_domains",
             "description_2014", "description_new"]
    result = result[[c for c in first if c in result.columns]].sort_values("gene_symbol")
    with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
        result.to_excel(writer, sheet_name="combined_RBPs", index=False)
    print(f"Wrote {len(result):,} genes to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit("Usage: python build_rbp_sheet.py 2014_census.xls filtered_candidates.xlsx output.xlsx")
    main(*sys.argv[1:])
