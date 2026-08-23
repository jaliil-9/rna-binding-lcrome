import csv
import io
import re
import sys
import time
from pathlib import Path

import requests

INPUT = Path(r"datasets\rbpdb\uniprot\rbpdb_human_current_uniprot_swissprot.csv")
OUTDIR = Path("rbpdb_uniprot_output")
OUTDIR.mkdir(exist_ok=True)

BASE = "https://rest.uniprot.org"
ENSG_RE = re.compile(r"ENSG\d{11}")
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "rbpdb-uniprot-fast/1.0 (research use)"})


def get(url, **kwargs):
    for attempt in range(5):
        try:
            r = SESSION.get(url, timeout=90, **kwargs)
            if r.status_code == 429:
                time.sleep(int(r.headers.get("Retry-After", "5")))
                continue
            r.raise_for_status()
            return r
        except requests.RequestException as exc:
            if attempt == 4:
                raise RuntimeError(f"Request failed: {url}\n{exc}") from exc
            time.sleep(2 ** attempt)


def post(url, **kwargs):
    for attempt in range(5):
        try:
            r = SESSION.post(url, timeout=90, **kwargs)
            if r.status_code == 429:
                time.sleep(int(r.headers.get("Retry-After", "5")))
                continue
            r.raise_for_status()
            return r
        except requests.RequestException as exc:
            if attempt == 4:
                raise RuntimeError(f"Request failed: {url}\n{exc}") from exc
            time.sleep(2 ** attempt)


def read_rbpdb(path):
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    try:
        dialect = csv.Sniffer().sniff(text[:5000], delimiters=",\t;")
    except csv.Error:
        dialect = csv.excel_tab

    rows = []
    for n, row in enumerate(csv.DictReader(text.splitlines(), dialect=dialect), 1):
        row = {f"rbpdb_{k.strip()}": (v or "").strip() for k, v in row.items()}
        hit = ENSG_RE.search(" ".join(row.values()))
        row["rbpdb_row_number"] = n
        row["ensembl_gene_id"] = hit.group(0) if hit else ""
        rows.append(row)
    return rows


def map_ensembl_to_swissprot(ensembl_ids):
    """Map ENSG gene IDs only to reviewed UniProtKB/Swiss-Prot."""
    response = post(
        f"{BASE}/idmapping/run",
        data={
            "from": "Ensembl",
            "to": "UniProtKB-Swiss-Prot",
            "ids": ",".join(ensembl_ids),
        },
    )
    job_id = response.json()["jobId"]
    print(f"UniProt mapping job: {job_id}", flush=True)

    while True:
        status = get(f"{BASE}/idmapping/status/{job_id}").json()
        if status.get("jobStatus") in {"RUNNING", "NEW"}:
            print("Waiting for UniProt mapping...", flush=True)
            time.sleep(3)
            continue
        if status.get("jobStatus"):
            raise RuntimeError(f"UniProt mapping failed: {status['jobStatus']}")
        break

    url = f"{BASE}/idmapping/results/{job_id}?format=tsv&size=500"
    mapping = {}

    while url:
        response = get(url)
        reader = csv.DictReader(io.StringIO(response.text), delimiter="\t")
        for record in reader:
            source = record.get("From", "")
            accession = record.get("To", "").split()[0]
            if source and accession:
                mapping.setdefault(source, []).append(accession)

        match = re.search(r'<([^>]+)>;\s*rel="next"', response.headers.get("Link", ""))
        url = match.group(1) if match else None

    return mapping


def fetch_metadata(accessions, batch_size=50):
    """
    Retrieve selected UniProt fields in small batches.
    batch_size=50 avoids excessive URL length.
    """
    fields = [
        "accession",
        "id",
        "reviewed",
        "protein_name",
        "gene_names",
        "organism_name",
        "organism_id",
        "length",
        "cc_function",
        "cc_subcellular_location",
        "go_id",
        "keyword",
        "ft_domain",
        "ft_region",
        "cc_ptm",
        "cc_disease",
        "cc_tissue_specificity",
        "sequence",
    ]

    all_rows = []

    for start in range(0, len(accessions), batch_size):
        batch = accessions[start:start + batch_size]

        # Compact valid UniProt query syntax.
        query = " OR ".join(f"accession:{acc}" for acc in batch)

        params = {
            "query": f"({query}) AND organism_id:9606",
            "format": "tsv",
            "fields": ",".join(fields),
            "size": batch_size,
        }

        print(
            f"Retrieving UniProt metadata: "
            f"{start + 1}-{start + len(batch)} / {len(accessions)}",
            flush=True,
        )

        response = get(f"{BASE}/uniprotkb/search", params=params)

        rows = list(csv.DictReader(
            io.StringIO(response.text),
            delimiter="\t"
        ))
        all_rows.extend(rows)

    return all_rows

def make_fasta(row):
    accession = row.get("Entry", "")
    entry_name = row.get("Entry Name", "")
    protein = row.get("Protein names", "")
    sequence = row.get("Sequence", "").replace(" ", "")

    header = f">{accession}|{entry_name} {protein}"
    wrapped = "\n".join(sequence[i:i + 60] for i in range(0, len(sequence), 60))
    return f"{header}\n{wrapped}\n"


def main():
    rbpdb_rows = read_rbpdb(INPUT)
    ensg_ids = sorted({x["ensembl_gene_id"] for x in rbpdb_rows if x["ensembl_gene_id"]})

    print(f"RBPDB rows read: {len(rbpdb_rows)}", flush=True)
    print(f"Unique Ensembl gene IDs: {len(ensg_ids)}", flush=True)

    mapping = map_ensembl_to_swissprot(ensg_ids)

    # Pick one reviewed accession per ENSG, preserving one-sequence-per-gene design.
    selected = {ensg: sorted(set(accs))[0] for ensg, accs in mapping.items()}
    accessions = sorted(set(selected.values()))

    print(f"ENSG IDs with reviewed UniProt mapping: {len(selected)}", flush=True)
    print(f"Selected unique Swiss-Prot accessions: {len(accessions)}", flush=True)
    print("Retrieving metadata in one batched request...", flush=True)

    metadata_rows = fetch_metadata(accessions)
    metadata_by_accession = {row["Entry"]: row for row in metadata_rows}

    output_rows = []
    fasta_records = []
    success = 0
    failed = 0

    for rbpdb_row in rbpdb_rows:
        ensg = rbpdb_row["ensembl_gene_id"]
        accession = selected.get(ensg)
        metadata = metadata_by_accession.get(accession) if accession else None

        if metadata:
            success += 1
            output_rows.append({
                **rbpdb_row,
                "retrieval_status": "success",
                "selected_uniprot_accession": accession,
                **{f"uniprot_{key.lower().replace(' ', '_')}": value
                   for key, value in metadata.items()},
            })
            fasta_records.append(make_fasta(metadata))
        else:
            failed += 1
            output_rows.append({
                **rbpdb_row,
                "retrieval_status": "no_reviewed_uniprot_mapping_or_metadata",
                "selected_uniprot_accession": accession or "",
            })

    csv_path = OUTDIR / "rbpdb_human_current_uniprot_swissprot.csv"
    fasta_path = OUTDIR / "rbpdb_human_current_uniprot_swissprot.fasta"

    columns = sorted({key for row in output_rows for key in row})
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(output_rows)

    with fasta_path.open("w", encoding="utf-8") as handle:
        handle.writelines(dict.fromkeys(fasta_records))

    print("\nFinished", flush=True)
    print(f"Successful retrievals: {success}", flush=True)
    print(f"Unsuccessful retrievals: {failed}", flush=True)
    print(f"CSV: {csv_path}", flush=True)
    print(f"FASTA: {fasta_path}", flush=True)


if __name__ == "__main__":
    main()