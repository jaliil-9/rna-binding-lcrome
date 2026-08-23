"""
Run self-hosted PlaToLoCo and export LCR intervals as CSV and FASTA files.

"""

import csv
import json
import sys
import time
from pathlib import Path

import requests


API = "http://127.0.0.1:5002/restapi"
POLL_SECONDS = 5

# PlaToLoCo response method -> output label.
METHOD_LABELS = {
    "SEG": "SEG_strict",
    "seg": "SEG_strict",
    "SEG_strict": "SEG_strict",
    "seg_strict": "SEG_strict",

    "CAST": "CAST",
    "cast": "CAST",

    "fLPS": "fLPS_strict",
    "FLPS": "fLPS_strict",
    "flps": "fLPS_strict",
    "fLPS_strict": "fLPS_strict",
    "flps_strict": "fLPS_strict",
}


def read_fasta(fasta_path):
    """Read FASTA exactly as text for submission to PlaToLoCo."""
    fasta_text = fasta_path.read_text().strip()
    if not fasta_text.startswith(">"):
        raise ValueError("Input must be a FASTA file whose first line starts with '>'.")
    return fasta_text + "\n"


def make_payload(fasta_text):
    """Create a PlaToLoCo request for SEG strict, CAST, and fLPS strict."""
    return {
        "name": "",
        "sequences": fasta_text,
        "methods": {
            "seg_default": False,
            "seg_intermediate": False,
            "seg_strict": True,
            "cast": True,
            "flps": False,
            "flps_strict": True,
            "simple": False,
            "gbsc": False,
        },
        "enrichment": {
            "pfam": False,
            "phobius": False,
            "aafrequency": False,
        },
        "params": {
            "seg": {"window": 12, "k1": 2.2, "k2": 2.5},
            "cast": {"threshold": 40, "matrix": 1},
            "flps": {
                "min_tract_len": 15,
                "max_tract_len": 500,
                "pval": 0.001,
                "regions": {"single": True, "multiple": True, "whole": False},
            },
            "simple": {},
            "gbsc": {},
        },
    }


def submit_job(payload):
    """Submit a job and return its PlaToLoCo token."""
    response = requests.put(f"{API}/query", json=payload, timeout=60)
    response.raise_for_status()
    return response.json()["token"]


def wait_for_job(token):
    """Poll PlaToLoCo until the job finishes or fails."""
    while True:
        response = requests.get(f"{API}/job/{token}", timeout=60)
        response.raise_for_status()
        status = response.json()["status"]
        print(f"Job status: {status}")

        if status == "FINISHED":
            return
        if status == "ERROR":
            raise RuntimeError("PlaToLoCo reported ERROR. Check: docker logs platoloco-api")

        time.sleep(POLL_SECONDS)


def get_rows(token):
    """Retrieve every protein result and return one row per LCR interval."""
    response = requests.get(f"{API}/proteins/{token}", timeout=60)
    response.raise_for_status()
    proteins = response.json()["proteins"]

    rows = []
    for index, protein in enumerate(proteins, start=1):
        protein_id = protein["id"]
        response = requests.get(f"{API}/proteins/{token}/{protein_id}", timeout=60)
        response.raise_for_status()
        details = response.json()

        accession = details.get("uniprot_id") or details["header"].split()[0].lstrip(">")
        sequence = details["sequence"]

        for result in details["data"]["wrapper"]:
            raw_method = result["method"]
            regions = result["regions"]

            print(
                f"{accession}: API method={raw_method}, "
                f"regions={len(regions)}"
            )

            output_method = METHOD_LABELS.get(raw_method)

            if output_method is None:
                raise ValueError(
                    f"Unexpected method label returned by PlaToLoCo: {raw_method}"
                )

            for region in regions:
                start = int(region["beg"])
                end = int(region["end"])
                rows.append(
                    {
                        "protein_id": accession,
                        "header": details["header"],
                        "method": output_method,
                        "start": start,
                        "end": end,
                        "length": end - start + 1,
                        "sequence": sequence[start - 1:end],
                        "description": region.get("description", ""),
                    }
                )

        print(f"Retrieved {index}/{len(proteins)} proteins")

    return rows


def write_outputs(rows, output_dir):
    """Write a CSV and FASTA file for SEG strict, CAST, and fLPS strict."""
    fields = [
        "protein_id", "header", "method", "start", "end",
        "length", "sequence", "description",
    ]

    for method in METHOD_LABELS.values():
        method_rows = [row for row in rows if row["method"] == method]

        with open(output_dir / f"{method}.csv", "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(method_rows)

        with open(output_dir / f"{method}.fasta", "w") as handle:
            for row in method_rows:
                handle.write(
                    f">{row['protein_id']} method={method} "
                    f"start={row['start']} end={row['end']} "
                    f"length={row['length']}\n{row['sequence']}\n"
                )


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python run_platoloco.py input.fasta output_directory")

    fasta_path = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = make_payload(read_fasta(fasta_path))
    (output_dir / "platoloco_request.json").write_text(json.dumps(payload, indent=2))

    token = submit_job(payload)
    print(f"Submitted job: {token}")
    wait_for_job(token)

    rows = get_rows(token)
    write_outputs(rows, output_dir)
    print(f"Detected {len(rows)} LCR intervals.")
    print(f"Outputs written to: {output_dir}")


if __name__ == "__main__":
    main()
