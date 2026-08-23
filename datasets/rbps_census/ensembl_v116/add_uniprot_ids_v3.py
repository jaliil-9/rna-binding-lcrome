#!/usr/bin/env python3
"""Annotate the canonical-proteome FASTA headers with UniProt accessions (v3).

v3 changes (after inspecting the r116 file):
- File layout confirmed: gene_stable_id, transcript_stable_id, protein_stable_id,
  xref, db_name, info_type, source_identity, xref_identity, linkage_type.
  Column detection is content-based (ENST column; accession column validated
  against the UniProt accession pattern; aborts if <80% valid).
- db_name-aware: when a transcript maps to multiple UniProt entries, the
  Uniprot/SWISSPROT row is preferred over Uniprot/SPTREMBL (master table is
  reviewed-centric). Counts reported either way.

Inputs
------
--fasta      ensembl_canonical_proteins.fa   (output of build_canonical_proteome.py)
--xref       Homo_sapiens.GRCh38.116.uniprot.tsv.gz
--selection  ensembl_canonical_selection.tsv (gene_id, transcript_id, rule)

Outputs
-------
<out>.fa               FASTA with ' uniprot:<ACC>' appended to each header
<out>_mapping.tsv      selection table + uniprot_accession column ("" if unmapped)
<out>_report.txt       QC counts + unmapped transcript list
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys

UNIPROT_ACC = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]"
    r"|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9])$"
)

def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def strip_version(ens_id: str) -> str:
    return str(ens_id).strip().split(".", 1)[0]

def load_xref(xref_path: str):
    """Return ({enst: acc}, stats). Swiss-Prot rows win over TrEMBL."""
    raw = []
    with open_maybe_gzip(xref_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        has_header = any(re.search(r"stable_id|xref|uniprot", c, re.I) for c in header)
        rows = [l.rstrip("\n").split("\t") for l in fh]
        if not has_header:
            rows = [header] + rows
            header = None
    n_cols = max(len(r) for r in rows[:50])
    enst_score = [0] * n_cols
    acc_score = [0] * n_cols
    for parts in rows[:50]:
        for i, v in enumerate(parts):
            v = v.strip()
            if re.fullmatch(r"ENST[0-9]+(?:\.\d+)?", v):
                enst_score[i] += 1
            if UNIPROT_ACC.fullmatch(v):
                acc_score[i] += 1
    tx_idx = max(range(n_cols), key=lambda i: enst_score[i])
    acc_idx = max(range(n_cols), key=lambda i: acc_score[i])
    db_idx = None
    if header:
        db_idx = next((i for i, c in enumerate(header) if re.search(r"db_name", c, re.I)), None)
    print(f"column sniffing: transcript={tx_idx} ({enst_score[tx_idx]}/50 ENST), "
          f"accession={acc_idx} ({acc_score[acc_idx]}/50 valid), "
          f"db_name={'absent' if db_idx is None else db_idx}")
    if acc_score[acc_idx] < 40:
        sys.exit(f"ABORT: column {acc_idx} does not validate as UniProt accessions "
                 f"({acc_score[acc_idx]}/50). Inspect the file and fix indices.")
    if tx_idx == acc_idx:
        sys.exit("ABORT: transcript and accession columns resolve to the same index.")

    tx_to_acc = {}
    tx_is_swissprot = {}
    n_rows = n_sp = multi = upgraded = 0
    for parts in rows:
        if len(parts) <= max(tx_idx, acc_idx):
            continue
        enst = strip_version(parts[tx_idx])
        acc = parts[acc_idx].strip()
        if not re.fullmatch(r"ENST[0-9]+", enst) or not UNIPROT_ACC.fullmatch(acc):
            continue
        n_rows += 1
        is_sp = db_idx is not None and len(parts) > db_idx and "SWISSPROT" in parts[db_idx].upper()
        n_sp += is_sp
        if enst not in tx_to_acc:
            tx_to_acc[enst] = acc
            tx_is_swissprot[enst] = is_sp
        else:
            multi += 1
            if is_sp and not tx_is_swissprot[enst]:
                tx_to_acc[enst] = acc          # upgrade TrEMBL -> Swiss-Prot
                tx_is_swissprot[enst] = True
                upgraded += 1
    print(f"xref: {n_rows} rows ({n_sp} Swiss-Prot) -> {len(tx_to_acc)} transcripts; "
          f"{multi} extra rows collapsed, {upgraded} upgraded to Swiss-Prot")
    return tx_to_acc, {"rows": n_rows, "swissprot_rows": n_sp, "multi": multi,
                       "upgraded": upgraded}

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--xref", required=True)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--out", default="ensembl_canonical_uniprot")
    args = ap.parse_args()

    tx_to_acc, stats = load_xref(args.xref)

    selected_enst = []
    n_mapped_sel = 0
    with open(args.selection) as fin, open(f"{args.out}_mapping.tsv", "w") as fout:
        fout.write(fin.readline().rstrip("\n") + "\tuniprot_accession\n")
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            enst = strip_version(parts[1])
            selected_enst.append(enst)
            acc = tx_to_acc.get(enst, "")
            n_mapped_sel += bool(acc)
            fout.write(line.rstrip("\n") + f"\t{acc}\n")

    tx_re = re.compile(r"transcript:(ENST[0-9]+)(?:\.\d+)?")
    n_in = n_annot = 0
    with open_maybe_gzip(args.fasta) as fin, open(f"{args.out}.fa", "w") as fout:
        for line in fin:
            if line.startswith(">"):
                n_in += 1
                m = tx_re.search(line)
                acc = tx_to_acc.get(m.group(1)) if m else None
                if acc:
                    n_annot += 1
                    line = line.rstrip("\n") + f" uniprot:{acc}\n"
            fout.write(line)

    unmapped = [t for t in selected_enst if t not in tx_to_acc]
    report = (
        f"xref rows parsed                     : {stats['rows']} "
        f"({stats['swissprot_rows']} Swiss-Prot)\n"
        f"transcripts with UniProt accession   : {len(tx_to_acc)}\n"
        f"extra rows collapsed                 : {stats['multi']} "
        f"({stats['upgraded']} upgraded TrEMBL->Swiss-Prot)\n"
        f"selected proteins                    : {len(selected_enst)}\n"
        f"selected proteins mapped to UniProt  : {n_mapped_sel} "
        f"({n_mapped_sel / max(1, len(selected_enst)):.1%})\n"
        f"FASTA records annotated              : {n_annot}/{n_in}\n"
    )
    print("\n" + report)
    with open(f"{args.out}_report.txt", "w") as fh:
        fh.write(report)
        if unmapped:
            fh.write("\nunmapped transcripts:\n" + "\n".join(sorted(unmapped)) + "\n")

if __name__ == "__main__":
    main()
