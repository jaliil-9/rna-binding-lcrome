"""
Combine Pfam27 and Pfam38 protein sets, then annotate proteins from UniProt.

Input:
    ensembl116_pfam38_selected_proteins.xlsx
        - sheets: Strict, Borderline
        - protein column: query_name

    rna_binding_proteins.xls
        - sheet: RBP table
        - protein column: protein id

Output:
    combined_rbp_pfam38_uniprot.xlsx
"""

from collections import defaultdict
from io import StringIO
import time
import re

import pandas as pd
import requests


PFAM38_FILE = "datasets/rbps_census/ensembl_v116/ensembl116_pfam38_selected_proteins.xlsx"
RBP_FILE = "datasets/rbps_census/gerstberg/rna_binding_proteins.xls"
OUTPUT_FILE = "datasets/rbps_census/combined_rbp_pfam38_uniprot.xlsx"

UNIPROT_API = "https://rest.uniprot.org"


def remove_ensembl_version(protein_id):
    """Convert ENSP00000xxxxx.3 to ENSP00000xxxxx."""
    if pd.isna(protein_id):
        return None
    return str(protein_id).strip().split(".")[0]


def load_protein_sets():
    """Load protein IDs from both datasets and assign source labels."""

    strict = pd.read_excel(PFAM38_FILE, sheet_name="Strict")
    borderline = pd.read_excel(PFAM38_FILE, sheet_name="Borderline")

    pfam38_ids = pd.concat(
        [strict["query_name"], borderline["query_name"]],
        ignore_index=True,
    ).map(remove_ensembl_version)

    rbp = pd.read_excel(RBP_FILE, sheet_name="RBP table")

    pfam27_ids = rbp["protein id"].map(remove_ensembl_version)

    pfam38_ids = set(pfam38_ids.dropna())
    pfam27_ids = set(pfam27_ids.dropna())

    all_ids = sorted(pfam27_ids | pfam38_ids)

    source_by_id = {}
    for protein_id in all_ids:
        in_pfam27 = protein_id in pfam27_ids
        in_pfam38 = protein_id in pfam38_ids

        if in_pfam27 and in_pfam38:
            source_by_id[protein_id] = "pfam27;pfam38"
        elif in_pfam27:
            source_by_id[protein_id] = "pfam27"
        else:
            source_by_id[protein_id] = "pfam38"

    print(f"Pfam27/RBP proteins: {len(pfam27_ids)}")
    print(f"Pfam38 proteins:     {len(pfam38_ids)}")
    print(f"Shared proteins:     {len(pfam27_ids & pfam38_ids)}")
    print(f"Unique proteins:     {len(all_ids)}")

    return all_ids, source_by_id


# ---------- UniProt mapping ----------

def submit_uniprot_mapping(ensembl_ids):
    """Submit Ensembl protein IDs to UniProt's ID-mapping service."""

    response = requests.post(
        f"{UNIPROT_API}/idmapping/run",
        data={
            "from": "Ensembl_Protein",
            "to": "UniProtKB",
            "ids": ",".join(ensembl_ids),
        },
        timeout=60,
    )
    response.raise_for_status()

    return response.json()["jobId"]


def wait_for_mapping(job_id):
    """Wait until UniProt finishes the ID-mapping job."""

    url = f"{UNIPROT_API}/idmapping/status/{job_id}"

    while True:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        status = response.json()

        if "jobStatus" not in status:
            return

        if status["jobStatus"] == "FINISHED":
            return

        if status["jobStatus"] in {"ERROR", "FAILED"}:
            raise RuntimeError(f"UniProt mapping failed: {status}")

        time.sleep(3)


def get_mapping_results(job_id):
    """Return versionless Ensembl protein ID -> UniProt accessions."""

    url = f"{UNIPROT_API}/idmapping/uniprotkb/results/{job_id}"
    params = {"format": "json", "size": 500}

    ensembl_to_uniprot = defaultdict(list)

    while url:
        response = requests.get(url, params=params, timeout=60)
        response.raise_for_status()
        data = response.json()

        for result in data.get("results", []):
            ensembl_id = remove_ensembl_version(result["from"])
            uniprot_record = result["to"]

            if isinstance(uniprot_record, dict):
                accession = uniprot_record["primaryAccession"]
            else:
                accession = str(uniprot_record)

            ensembl_to_uniprot[ensembl_id].append(accession)

        url = response.links.get("next", {}).get("url")
        params = None

    return ensembl_to_uniprot


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
        r"[OPQ][0-9][A-Z0-9]{3}[0-9]"           
        r"|"
        r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9]"        
        r"|"
        r"[A-Z0-9]{10}"                         
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
    Download UniProt annotations.

    Using batches of 50 accessions to avoid problematic long query URLs.
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

    return annotations


def select_best_uniprot_match(ensembl_to_uniprot, annotations):
    """Select one annotated UniProt entry per Ensembl protein."""

    # Keep the UniProt accession as an ordinary column in every record.
    annotation_by_accession = {
        row["uniprot_accession"]: row
        for row in (
            annotations
            .drop_duplicates(subset="uniprot_accession")
            .to_dict("records")
        )
    }

    rows = []

    for ensembl_id, candidate_accessions in ensembl_to_uniprot.items():
        ensembl_id = remove_ensembl_version(ensembl_id)

        candidate_rows = [
            annotation_by_accession[accession]
            for accession in candidate_accessions
            if accession in annotation_by_accession
        ]

        if not candidate_rows:
            continue

        reviewed_rows = [
            row for row in candidate_rows
            if str(row["uniprot_reviewed"]).lower() == "reviewed"
        ]

        best_row = reviewed_rows[0] if reviewed_rows else candidate_rows[0]

        rows.append({
            "ensembl_protein_id": ensembl_id,
            **best_row,
        })

    return pd.DataFrame(rows)


def main():
    ensembl_ids, source_by_id = load_protein_sets()

    print("\nSubmitting IDs to UniProt...")
    job_id = submit_uniprot_mapping(ensembl_ids)

    print("Waiting for UniProt mapping...")
    wait_for_mapping(job_id)

    ensembl_to_uniprot = get_mapping_results(job_id)

    all_accessions = {
        accession
        for accessions in ensembl_to_uniprot.values()
        for accession in accessions
    }

    print(f"UniProt accessions retrieved: {len(all_accessions)}")
    annotations = download_uniprot_annotations(all_accessions)

    combined = pd.DataFrame(
        {
            "ensembl_protein_id": ensembl_ids,
            "source": [source_by_id[protein_id] for protein_id in ensembl_ids],
        }
    )

    submitted = set(ensembl_ids)
    mapped = set(ensembl_to_uniprot)

    print(f"Submitted but not mapped:   {len(submitted - mapped)}")
    print("Example unmapped IDs:", sorted(submitted - mapped)[:10])
    best_matches = select_best_uniprot_match(ensembl_to_uniprot, annotations)
    print("Best-match columns:", best_matches.columns.tolist())

    combined["ensembl_protein_id"] = (
        combined["ensembl_protein_id"]
        .astype(str)
        .str.strip()
        .str.split(".")
        .str[0]
    )

    best_matches["ensembl_protein_id"] = (
        best_matches["ensembl_protein_id"]
        .astype(str)
        .str.strip()
        .str.split(".")
        .str[0]
    )

    # Prevent duplicate mapping rows from multiplying rows in the final table.
    best_matches = best_matches.drop_duplicates(
        subset="ensembl_protein_id",
        keep="first",
    )

    shared_ids = set(combined["ensembl_protein_id"]) & set(
        best_matches["ensembl_protein_id"]
    )

    print(f"IDs shared by both merge tables: {len(shared_ids)}")
    print("Combined example:", combined["ensembl_protein_id"].head().tolist())
    print("Mapped example:", best_matches["ensembl_protein_id"].head().tolist())

    combined = combined.merge(
        best_matches,
        on="ensembl_protein_id",
        how="left",
        validate="one_to_one",
    )

    mapped_count = combined["uniprot_accession"].notna().sum()
    unmapped_count = combined["uniprot_accession"].isna().sum()

    print(f"Rows with UniProt accession: {mapped_count}")
    print(f"Rows without UniProt accession: {unmapped_count}")

    assert mapped_count > 1700, (
        "Merge failed: expected 1,790 UniProt-mapped proteins."
    )

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

    combined = combined.reindex(columns=output_columns)
    combined.to_excel(OUTPUT_FILE, sheet_name="Combined", index=False)

    unmapped = combined["uniprot_accession"].isna().sum()

    print(f"\nOutput written to: {OUTPUT_FILE}")
    print(f"Rows written:       {len(combined)}")
    print(f"Without UniProt map: {unmapped}")

    unmapped = combined.loc[
        combined["uniprot_accession"].isna(),
        ["ensembl_protein_id", "source"],
    ]

    unmapped.to_excel(
        "datasets/rbps_census/unmapped_ensembl_proteins.xlsx",
        index=False,
    )

    print("Saved unmapped proteins to: unmapped_ensembl_proteins.xlsx")


if __name__ == "__main__":
    main()