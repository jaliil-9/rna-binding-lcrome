#!/usr/bin/env python3
"""Finalize the Phase 2.5 background proteome and RBP label set.

Actions:
1. Label patch: genes in the background that are RBPs per the master table but
   under a different (current) accession -> their background accessions are
   added to the RBP label set (label_source = gene_patch).
2. Append: master rows whose gene is absent from the background get their
   master-table sequence appended as a new FASTA record (label_source =
   master_appended). Uses the exact sequence the RBP-side LCRs were called on.
3. Emit background_labels.tsv: every background accession -> is_rbp + source.
4. Dump the still-unresolved rows (no_ensp_mapping) for UniProt ID mapping.

Inputs: resolve_absent TSV, reconciliation TSV, master xlsx, background FASTA.
"""
from __future__ import annotations

import argparse
import re

import pandas as pd

def read_fasta(path):
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line, []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolve", required=True)
    ap.add_argument("--rec", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", default="finalized")
    args = ap.parse_args()

    resolve = pd.read_csv(args.resolve, sep="\t")
    rec = pd.read_csv(args.rec, sep="\t")
    master = pd.read_excel(args.master, sheet_name="Combined")
    master.columns = [re.sub(r"[^a-z0-9]+", "", str(c).lower()) for c in master.columns]

    bg_records = list(read_fasta(args.fasta))
    bg_accs = []
    for h, _ in bg_records:
        rid = h[1:].split()[0]
        bg_accs.append(rid.split("|")[1] if rid.startswith(("sp|", "tr|")) else rid)

    # --- 1. label patch
    patch = resolve[resolve["bucket"] == "gene_in_background_unlabeled"].copy()
    patch = patch.assign(
        covering_accession=patch["covering_accession"].fillna("").str.split(";")
    ).explode("covering_accession")
    patch = patch[patch["covering_accession"] != ""]
    patch_accs = set(patch["covering_accession"])
    patch[["uniprot_accession", "gene_id", "covering_accession"]].to_csv(
        f"{args.out}_label_patch.tsv", sep="\t", index=False)

    # --- 2. append master sequences for genes absent from background
    missing = resolve[resolve["bucket"] == "gene_not_in_background"].copy()
    m_sub = master.dropna(subset=["uniprotaccession"]).drop_duplicates("uniprotaccession")
    missing = missing.merge(
        m_sub[["uniprotaccession", "uniprotentryname", "uniprotsequence"]],
        left_on="uniprot_accession", right_on="uniprotaccession", how="left")
    appended = []
    out_fasta = f"{args.out}_background.fa"
    with open(out_fasta, "w") as fh:
        for h, s in bg_records:
            fh.write(h + "\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
        for _, r in missing.iterrows():
            seq = str(r.get("uniprotsequence", "")).strip().upper()
            acc = str(r["uniprot_accession"]).strip()
            if not seq or seq == "NAN" or not acc or acc == "nan":
                continue
            name = str(r.get("uniprotentryname") or f"{acc}_HUMAN")
            fh.write(f">tr|{acc}|{name} master-appended RBP sequence\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
            appended.append(acc)

    # --- 3. label file
    matched_accs = set(rec.loc[rec["status"] == "accession_present", "uniprot_accession"])
    rows = []
    for acc in bg_accs:
        src = ("master_direct" if acc in matched_accs
               else "gene_patch" if acc in patch_accs else "")
        rows.append({"uniprot_accession": acc, "is_rbp": bool(src), "label_source": src})
    for acc in appended:
        rows.append({"uniprot_accession": acc, "is_rbp": True,
                     "label_source": "master_appended"})
    labels = pd.DataFrame(rows)
    labels.to_csv(f"{args.out}_background_labels.tsv", sep="\t", index=False)

    # --- 4. unresolved dump
    unres = resolve[resolve["bucket"] == "no_ensp_mapping"].copy()
    unres = unres.merge(
        m_sub[["uniprotaccession", "uniprotentryname", "uniprotreviewed",
               "uniprotlength", "uniprotsequence", "source"]],
        left_on="uniprot_accession", right_on="uniprotaccession", how="left")
    unres.to_csv(f"{args.out}_unresolved.tsv", sep="\t", index=False)
    has_acc = unres["uniprot_accession"].notna() & (unres["uniprot_accession"].astype(str) != "nan")
    has_seq = unres["uniprotsequence"].notna()

    n_rbp = int(labels["is_rbp"].sum())
    report = (
        f"background proteins (incl. appended) : {len(labels)}\n"
        f"  appended master sequences          : {len(appended)}\n"
        f"RBP-labeled                          : {n_rbp} "
        f"({n_rbp / len(labels):.1%} of background)\n"
        f"  via direct accession match         : {len(matched_accs)}\n"
        f"  via gene-level patch               : {len(patch_accs)}\n"
        f"  via appended sequences             : {len(appended)}\n"
        f"unresolved rows (need ID mapping)    : {len(unres)}\n"
        f"  of which have an accession         : {int(has_acc.sum())}\n"
        f"  of which have a sequence           : {int(has_seq.sum())}\n"
    )
    print("\n" + report)
    with open(f"{args.out}_report.txt", "w") as fh:
        fh.write(report)

if __name__ == "__main__":
    main()
