#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
import pandas as pd


def read_pfam_dat(path):
    records, name, accession = [], None, None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#=GF ID"):
                name = line.split(maxsplit=2)[2].strip()
            elif line.startswith("#=GF AC"):
                accession = line.split(maxsplit=2)[2].strip().split(".")[0]
            elif line.strip() == "//":
                if name and accession:
                    records.append((name, accession))
                name, accession = None, None
    return pd.DataFrame(records, columns=["pfam_name_38_2", "accession"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("domains_xlsx", help="Excel file containing historical Pfam names in column 1")
    parser.add_argument("pfam_dat", help="Pfam 38.2 Pfam-A.hmm.dat file")
    parser.add_argument("-o", "--output", default="pfam38_2_domain_mapping.xlsx")
    args = parser.parse_args()

    domains = pd.read_excel(args.domains_xlsx, usecols=[0])
    domains.columns = ["name"]
    domains["name"] = domains["name"].astype(str).str.strip()
    domains = domains[domains["name"].ne("") & domains["name"].ne("nan")].drop_duplicates()

    pfam = read_pfam_dat(args.pfam_dat)
    result = domains.merge(pfam, left_on="name", right_on="pfam_name_38_2", how="left")
    result["status"] = result["accession"].notna().map({True: "found", False: "not found"})
    result = result[["name", "accession", "status"]]
    result.to_excel(args.output, index=False)


if __name__ == "__main__":
    main()
