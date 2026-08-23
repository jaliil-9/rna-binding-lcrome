import argparse
import re
import time
import pandas as pd
import requests

GO_TERMS = [
    "GO:0003723", "GO:0044822", "GO:0003729", "GO:0003725",
    "GO:0003727", "GO:0019843", "GO:0000049", "GO:0017069",
    "GO:0030515", "GO:0106222", "GO:0061980", "GO:0003730",
    "GO:0048027", "GO:0000339", "GO:0002151", "GO:1990247",
    "GO:0008135",
]
KEYWORDS = [
    "RNA", "rRNA", "tRNA", "mRNA", "ribosomal", "ribonuclease",
    "ribonucleoprotein", "spliceosom", "RNA polymerase", "helicase",
    "RNA methyl",
]


def read_pfam_dat(path):
    records, name, accession, description = [], None, None, None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#=GF ID"):
                name = line.split(maxsplit=2)[2].strip()
            elif line.startswith("#=GF AC"):
                accession = line.split(maxsplit=2)[2].strip().split(".")[0]
            elif line.startswith("#=GF DE"):
                description = line.split(maxsplit=2)[2].strip()
            elif line.strip() == "//":
                if name and accession:
                    records.append((accession, name, description or ""))
                name, accession, description = None, None, None
    return pd.DataFrame(records, columns=["accession", "name", "description"])


def get_go_candidates(go_term):
    url = "https://www.ebi.ac.uk/interpro/api/entry/interpro/"
    params = {"go_term": go_term, "page_size": 200}
    rows = []

    while url:
        response = requests.get(
            url,
            params=params,
            headers={"Accept": "application/json"},
            timeout=60
        )

        if response.status_code == 204:
            break

        response.raise_for_status()
        data = response.json()

        for entry in data.get("results", []):
            pfam_members = entry["metadata"].get("member_databases", {}).get("pfam", {})
            for accession, name in pfam_members.items():
                rows.append((accession, name))

        url = data.get("next")
        params = None
        time.sleep(0.1)

    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mapping_xlsx", default="datasets/rbps_census/gerstberg/rna_binding_domains.xlsx", help="Updated historical-domain mapping Excel file")
    parser.add_argument("pfam38_dat", default="datasets/rbps_census/ensembl_v116/Pfam-A/Pfam-A38.2.hmm.dat", help="Pfam 38.2 Pfam-A.hmm.dat")
    parser.add_argument("-o", "--output", default="datasets/rbps_census/ensembl_v116/pfam38_new_rna_candidates.xlsx")
    args = parser.parse_args()

    historical = pd.read_excel(args.mapping_xlsx)
    historical_accessions = set(historical["accession"].dropna().astype(str))
    pfam = read_pfam_dat(args.pfam38_dat)

    go_rows = []
    for go_term in GO_TERMS:
        print(f"Querying {go_term}")
        for accession, name in get_go_candidates(go_term):
            go_rows.append((accession, name, "go"))
    go = pd.DataFrame(go_rows, columns=["accession", "name", "source"])

    pattern = re.compile("|".join(re.escape(k) for k in KEYWORDS), re.IGNORECASE)
    metadata = pfam[pfam["description"].str.contains(pattern, na=False)].copy()
    metadata = metadata[["accession", "name"]]
    metadata["source"] = "metadata"

    candidates = pd.concat([go, metadata], ignore_index=True)
    candidates = candidates.groupby(["accession", "name"], as_index=False)["source"].agg(
        lambda x: ";".join(sorted(set(x)))
    )
    candidates = candidates[~candidates["accession"].isin(historical_accessions)]
    candidates.sort_values(["source", "accession"]).to_excel(args.output, index=False)


if __name__ == "__main__":
    main()
