import argparse
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
    return pd.DataFrame(records, columns=["pfam_name", "accession"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mapping_xlsx", default="datasets/rbps_census/ensembl_v116/pfam38.2_rbd_candidates.xlsx", help="Output from map_pfam38_domains.py")
    parser.add_argument("pfam27_dat", default="datasets/rbps_census/ensembl_v116/Pfam-A/Pfam-A27.hmm.dat", help="Pfam 27 Pfam-A.hmm.dat")
    parser.add_argument("pfam38_dat", default="datasets/rbps_census/ensembl_v116/Pfam-A/Pfam-A38.2.hmm.dat", help="Pfam 38.2 Pfam-A.hmm.dat")
    parser.add_argument("-o", "--output", default="datasets/rbps_census/gerstberg/pfam_38_vs_27_domains_check.xlsx")
    args = parser.parse_args()

    mapping = pd.read_excel(args.mapping_xlsx)
    pfam27 = read_pfam_dat(args.pfam27_dat).rename(
        columns={"pfam_name": "pfam27_name", "accession": "pfam27_accession"}
    )
    pfam38 = read_pfam_dat(args.pfam38_dat).rename(
        columns={"pfam_name": "pfam38_name", "accession": "pfam38_accession"}
    )

    missing_idx = mapping.index[mapping["status"].eq("not found")]

    old = mapping.loc[missing_idx, ["name"]].merge(
        pfam27, left_on="name", right_on="pfam27_name", how="left"
    )
    old = old.merge(
        pfam38, left_on="pfam27_accession", right_on="pfam38_accession", how="left"
    )
    old.index = missing_idx

    mapping.loc[old.index, "accession"] = old["pfam27_accession"]

    renamed = old["pfam27_accession"].notna() & old["pfam38_accession"].notna()
    dead = old["pfam27_accession"].notna() & old["pfam38_accession"].isna()
    uncertain = old["pfam27_accession"].isna()

    mapping.loc[old.index[renamed], "status"] = "renamed"
    mapping.loc[old.index[dead], "status"] = "dead"
    mapping.loc[old.index[uncertain], "status"] = "uncertain"

    mapping.to_excel(args.output, index=False)


if __name__ == "__main__":
    main()
