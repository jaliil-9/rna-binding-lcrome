#!/usr/bin/env python3
"""
Add RBPDB and MODOMICS protein entries to the combined Pfam27/Pfam38/UniProt
RBP table, skipping anything already present to avoid duplication.

Input:
    combined_rbp_pfam38_uniprot.xlsx
        - sheet: Combined
        - output of combine_rbps.py; key column "ensembl_protein_id",
          annotation column "uniprot_accession"

    RBPDB_v1_3_1_proteins_human_2012-11-21.csv
        - RBPDB "proteins, human" export, no header row, 13 columns:
          id, ensembl_gene_id, creation_date, update_date, gene_symbol,
          gene_description, species, taxid, domains, aliases, col11,
          col12, col13
        - "ensembl_gene_id" is an ENSG (gene), not an ENSP (protein), so
          RBPDB entries are mapped to UniProt via Ensembl *gene* ID rather
          than Ensembl *protein* ID.

    protein_sequences.fasta
        - MODOMICS protein sequences. Header format:
          > acronym:<name>|Uniprot:<accession or None>|gi:<...>|
            cog:<...>|enzyme_type:<...>|species:<species>
        - Multi-species; restricted here to species == "Homo sapiens" (to
          match the rest of this human RBP census) and to headers with a
          resolvable UniProt accession (some are "Uniprot:None").

Output:
    combined_rbp_pfam38_rbpdb_modomics_uniprot.xlsx
"""

from collections import defaultdict
from io import StringIO
import re
import time

import pandas as pd
import requests


# ---------- File names ----------

COMBINED_FILE = "datasets/rbps_census/combined_rbp_pfam38_uniprot.xlsx"
RBPDB_FILE = "datasets/rbpdb/RBPDB_v1.3.1_human_2012-11-21_CSV/RBPDB_v1.3.1_proteins_human_2012-11-21.csv"
MODOMICS_FILE = "datasets/modomics/protein_sequences.fasta"
OUTPUT_FILE = "datasets/rbps_census/combined_rbp_pfam38_rbpdb_modomics_uniprot.xlsx"
UNMAPPED_RBPDB_FILE = "unmapped_rbpdb_genes.xlsx"

UNIPROT_API = "https://rest.uniprot.org"


# ---------- Existing combined table ----------

def load_existing_combined():
    """Load the pfam27/pfam38 combined table and its known UniProt accessions."""

    combined = pd.read_excel(COMBINED_FILE, sheet_name="Combined")

    known_accessions = set(
        combined["uniprot_accession"].dropna().astype(str).str.strip().str.upper()
    )

    print(f"Existing combined rows:        {len(combined)}")
    print(f"Existing UniProt accessions:   {len(known_accessions)}")

    return combined, known_accessions


# ---------- RBPDB ----------

RBPDB_COLUMNS = [
    "id", "ensembl_gene_id", "creation_date", "update_date", "gene_symbol",
    "gene_description", "species", "taxid", "domains", "aliases",
    "col11", "col12", "col13",
]


def load_rbpdb_gene_ids():
    """Load versionless Ensembl gene IDs (ENSG) from the RBPDB human export."""

    rbpdb = pd.read_csv(RBPDB_FILE, header=None, names=RBPDB_COLUMNS)

    rbpdb = rbpdb[rbpdb["species"].astype(str).str.strip() == "Homo sapiens"]

    gene_ids = (
        rbpdb["ensembl_gene_id"]
        .dropna()
        .astype(str)
        .str.strip()
        .str.split(".").str[0]
    )
    gene_ids = sorted(set(gene_ids) - {""})

    print(f"RBPDB human genes:              {len(gene_ids)}")

    return gene_ids


# ---------- MODOMICS ----------

FASTA_HEADER_RE = re.compile(
    r"acronym:(?P<acronym>[^|]*)\|"
    r"Uniprot:(?P<uniprot>[^|]*)\|"
    r"gi:(?P<gi>[^|]*)\|"
    r"cog:(?P<cog>[^|]*)\|"
    r"enzyme_type:(?P<enzyme_type>[^|]*)\|"
    r"species:(?P<species>.*)"
)


def load_modomics_accessions():
    """
    Parse MODOMICS FASTA header lines and return the UniProt accessions of
    human entries that have a resolvable UniProt accession.
    """

    accessions = set()
    total_headers = 0
    human_headers = 0

    with open(MODOMICS_FILE) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue

            total_headers += 1
            match = FASTA_HEADER_RE.search(line.strip())
            if not match:
                continue

            if match.group("species").strip() != "Homo sapiens":
                continue
            human_headers += 1

            uniprot = match.group("uniprot").strip()
            if not uniprot or uniprot.lower() == "none":
                continue

            accessions.add(uniprot.upper())

    print(f"MODOMICS FASTA headers:         {total_headers}")
    print(f"MODOMICS human headers:         {human_headers}")
    print(f"MODOMICS human w/ UniProt ID:   {len(accessions)}")

    return accessions


# ---------- UniProt ID mapping (generic) ----------

def run_uniprot_id_mapping(ids, from_db, to_db):
    """
    Submit an ID list to UniProt's ID-mapping service and return
    {from_id: [to_id, ...]}. Works for any from/to database pair, e.g.
    "Ensembl" -> "UniProtKB" or "UniProtKB" -> "Ensembl_Protein".
    """

    if not ids:
        return {}

    response = requests.post(
        f"{UNIPROT_API}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
        timeout=60,
    )
    if not response.ok:
        raise RuntimeError(
            f"UniProt ID-mapping submission failed ({response.status_code}).\n"
            f"From: {from_db}\n"
            f"To: {to_db}\n"
            f"Response: {response.text}"
        )

    job_id = response.json()["jobId"]

    status_url = f"{UNIPROT_API}/idmapping/status/{job_id}"
    while True:
        response = requests.get(status_url, timeout=60)
        response.raise_for_status()
        status = response.json()

        if "jobStatus" not in status or status["jobStatus"] == "FINISHED":
            break
        if status["jobStatus"] in {"ERROR", "FAILED"}:
            raise RuntimeError(f"UniProt mapping failed: {status}")

        time.sleep(3)

    # UniProtKB targets carry full records and need the dedicated
    # "uniprotkb/results" endpoint; other targets (e.g. Ensembl_Protein)
    # use the plain "results" endpoint.
    results_endpoint = "uniprotkb/results" if to_db == "UniProtKB" else "results"
    url = f"{UNIPROT_API}/idmapping/{results_endpoint}/{job_id}"
    params = {"format": "json", "size": 500}

    mapping = defaultdict(list)
    while url:
        response = requests.get(url, params=params, timeout=60)
        response.raise_for_status()
        data = response.json()

        for result in data.get("results", []):
            to_value = result["to"]
            to_id = (
                to_value["primaryAccession"]
                if isinstance(to_value, dict)
                else str(to_value)
            )
            mapping[result["from"]].append(to_id)

        url = response.links.get("next", {}).get("url")
        params = None

    return mapping


# ---------- UniProt annotation download ----------

def chunks(items, size=200):
    """Yield consecutive chunks from a list."""
    for start in range(0, len(items), size):
        yield items[start:start + size]


def is_valid_uniprot_accession(accession):
    """
    Accept standard UniProtKB primary accessions, e.g. P12345, Q9Y6K9,
    A0A024R214. Excludes malformed mapping outputs.
    """
    accession = str(accession).strip().upper()

    pattern = (
        r"^(?:"
        r"[OPQ][0-9][A-Z0-9]{3}[0-9]"           # 6-character format
        r"|"
        r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9]"        # 6-character format
        r"|"
        r"[A-Z0-9]{10}"                         # 10-character format
        r")$"
    )

    return bool(re.fullmatch(pattern, accession))


def fetch_uniprot_batch(accessions, fields):
    """Download one accession batch from UniProt."""

    query = " OR ".join(f"accession:{accession}" for accession in accessions)

    response = requests.get(
        f"{UNIPROT_API}/uniprotkb/search",
        params={
            "query": f"({query})",
            "format": "tsv",
            "fields": ",".join(fields),
            "size": 500,
        },
        timeout=120,
    )

    if not response.ok:
        raise RuntimeError(
            f"UniProt request failed ({response.status_code}).\n"
            f"Accessions: {', '.join(accessions)}\n"
            f"Response: {response.text}"
        )

    return pd.read_csv(StringIO(response.text), sep="\t")


def download_uniprot_annotations(accessions):
    """
    Download UniProt annotations robustly.

    Uses batches of 50 accessions to avoid problematic long query URLs.
    If a batch fails, it tries every accession individually and reports
    accession(s) that UniProt rejects.
    """

    fields = [
        "accession",
        "id",
        "protein_name",
        "reviewed",
        "length",
        "sequence",
        "go_id",
        "go_p",
        "go_f",
        "cc_subcellular_location",
        "cc_domain",
    ]

    clean_accessions = sorted({
        str(accession).strip().upper()
        for accession in accessions
        if is_valid_uniprot_accession(accession)
    })

    invalid_accessions = sorted({
        str(accession).strip()
        for accession in accessions
        if not is_valid_uniprot_accession(accession)
    })

    if invalid_accessions:
        print(
            f"Skipped {len(invalid_accessions)} invalid accession(s): "
            f"{invalid_accessions[:10]}"
        )

    if not clean_accessions:
        return pd.DataFrame(columns=["uniprot_accession"])

    tables = []
    failed_accessions = []

    for accession_chunk in chunks(clean_accessions, size=50):
        try:
            tables.append(fetch_uniprot_batch(accession_chunk, fields))

        except RuntimeError as error:
            print(f"\nBatch failed; retrying accessions individually.\n{error}")

            for accession in accession_chunk:
                try:
                    tables.append(fetch_uniprot_batch([accession], fields))
                except RuntimeError:
                    failed_accessions.append(accession)

    if failed_accessions:
        print(
            f"\nUniProt could not retrieve {len(failed_accessions)} accession(s): "
            f"{failed_accessions}"
        )

    if not tables:
        raise RuntimeError("No UniProt annotations could be downloaded.")

    annotations = pd.concat(tables, ignore_index=True)

    annotations = annotations.rename(
        columns={
            "Entry": "uniprot_accession",
            "Entry Name": "uniprot_entry_name",
            "Protein names": "uniprot_protein_name",
            "Reviewed": "uniprot_reviewed",
            "Length": "uniprot_length",
            "Sequence": "uniprot_sequence",
            "Gene Ontology IDs": "go_ids",
            "Gene Ontology (biological process)": "go_biological_process",
            "Gene Ontology (molecular function)": "go_molecular_function",
            "Subcellular location [CC]": "subcellular_location",
            "Domain [CC]": "uniprot_domains",
        }
    )

    annotations["uniprot_accession"] = (
        annotations["uniprot_accession"].astype(str).str.upper()
    )

    return annotations


def select_best_candidate(candidate_accessions, annotation_by_accession):
    """Pick one annotated UniProt accession from a list of candidates,
    preferring a reviewed (Swiss-Prot) entry."""

    candidates = [
        accession for accession in candidate_accessions
        if accession in annotation_by_accession
    ]
    if not candidates:
        return None

    reviewed = [
        accession for accession in candidates
        if str(annotation_by_accession[accession]["uniprot_reviewed"]).lower()
        == "reviewed"
    ]

    return reviewed[0] if reviewed else candidates[0]


# ---------- Main ----------

def main():
    combined, known_accessions = load_existing_combined()

    # --- RBPDB: gene -> UniProt (candidate list) ---
    print("\nMapping RBPDB Ensembl gene IDs to UniProt...")
    rbpdb_gene_ids = load_rbpdb_gene_ids()
    gene_to_uniprot = run_uniprot_id_mapping(
        rbpdb_gene_ids, from_db="Ensembl", to_db="UniProtKB"
    )

    unmapped_genes = sorted(set(rbpdb_gene_ids) - set(gene_to_uniprot))
    print(f"RBPDB genes not mapped to UniProt: {len(unmapped_genes)}")

    rbpdb_candidate_accessions = {
        accession.upper()
        for accessions in gene_to_uniprot.values()
        for accession in accessions
    }

    # --- MODOMICS: accession already known from FASTA header ---
    print("\nParsing MODOMICS FASTA...")
    modomics_accessions = load_modomics_accessions()

    # --- Download annotations for every candidate accession up front, so
    #     "select the reviewed one" can be decided from real data. ---
    all_candidate_accessions = rbpdb_candidate_accessions | modomics_accessions
    print(f"\nCandidate accessions to annotate: {len(all_candidate_accessions)}")

    annotations = download_uniprot_annotations(all_candidate_accessions)
    annotation_by_accession = {
        row["uniprot_accession"]: row
        for row in annotations.drop_duplicates(subset="uniprot_accession").to_dict("records")
    }

    # --- RBPDB: resolve one accession per gene, then keep only the ones
    #     that are genuinely new (not already in the combined table). ---
    rbpdb_new_accessions = set()
    for gene_id, candidates in gene_to_uniprot.items():
        best = select_best_candidate(
            [a.upper() for a in candidates], annotation_by_accession
        )
        if best and best not in known_accessions:
            rbpdb_new_accessions.add(best)

    # --- MODOMICS: already resolved to a single accession per entry. ---
    modomics_new_accessions = {
        accession for accession in modomics_accessions
        if accession in annotation_by_accession and accession not in known_accessions
    }

    print(f"\nNew accessions from RBPDB:    {len(rbpdb_new_accessions)}")
    print(f"New accessions from MODOMICS: {len(modomics_new_accessions)}")

    new_accessions = rbpdb_new_accessions | modomics_new_accessions
    print(f"New accessions total (union): {len(new_accessions)}")

    if not new_accessions:
        print("\nNothing new to add; combined table already covers RBPDB and MODOMICS.")
        combined.to_excel(OUTPUT_FILE, sheet_name="Combined", index=False)
        return

    # --- Best-effort reverse mapping to an Ensembl protein ID, so new rows
    #     line up with the existing "ensembl_protein_id" column where possible. ---
    print("\nReverse-mapping new UniProt accessions to Ensembl protein IDs...")
    uniprot_to_ensp = run_uniprot_id_mapping(
        sorted(new_accessions),
        from_db="UniProtKB_AC-ID",
        to_db="Ensembl_Protein",
    )

    def source_for(accession):
        in_rbpdb = accession in rbpdb_new_accessions
        in_modomics = accession in modomics_new_accessions
        if in_rbpdb and in_modomics:
            return "rbpdb;modomics"
        return "rbpdb" if in_rbpdb else "modomics"

    def ensp_for(accession):
        ensp_candidates = uniprot_to_ensp.get(accession)
        if not ensp_candidates:
            return None
        return sorted({remove_ensembl_version(e) for e in ensp_candidates})[0] # type: ignore

    def remove_ensembl_version(protein_id):
        if pd.isna(protein_id):
            return None
        return str(protein_id).strip().split(".")[0]

    new_rows = []
    for accession in sorted(new_accessions):
        row = dict(annotation_by_accession[accession])
        row["ensembl_protein_id"] = ensp_for(accession)
        row["source"] = source_for(accession)
        new_rows.append(row)

    new_df = pd.DataFrame(new_rows)

    output_columns = [
        "ensembl_protein_id",
        "source",
        "uniprot_accession",
        "uniprot_entry_name",
        "uniprot_protein_name",
        "uniprot_reviewed",
        "uniprot_length",
        "uniprot_sequence",
        "go_ids",
        "go_biological_process",
        "go_molecular_function",
        "subcellular_location",
        "uniprot_domains",
    ]

    new_df = new_df.reindex(columns=output_columns)
    combined = combined.reindex(columns=output_columns)

    final = pd.concat([combined, new_df], ignore_index=True)

    print(f"\nRows before: {len(combined)}")
    print(f"Rows added:  {len(new_df)}")
    print(f"Rows after:  {len(final)}")

    final.to_excel(OUTPUT_FILE, sheet_name="Combined", index=False)
    print(f"\nOutput written to: {OUTPUT_FILE}")

    if unmapped_genes:
        pd.DataFrame({"ensembl_gene_id": unmapped_genes}).to_excel(
            UNMAPPED_RBPDB_FILE, index=False
        )
        print(f"Saved unmapped RBPDB genes to: {UNMAPPED_RBPDB_FILE}")


if __name__ == "__main__":
    main()
