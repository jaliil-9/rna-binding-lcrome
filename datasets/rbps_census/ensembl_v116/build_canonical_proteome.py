#!/usr/bin/env python3
"""Build a one-protein-per-gene background proteome from Ensembl.

Selection rule per protein-coding gene:
    1. MANE Select transcript, if annotated;
    2. otherwise the Ensembl canonical transcript.

Genes with neither tag are reported but excluded (rare; typically
non-reference/patch genes). The pep.all FASTA is then filtered to the
selected transcripts.

Inputs
------
--gtf    Homo_sapiens.GRCh38.<release>.gtf.gz  (transcript lines carry the tags)
--fasta  Homo_sapiens.GRCh38.pep.all.fa.gz

Outputs
-------
<out>.canonical_selection.tsv   one row per selected transcript
<out>.canonical_proteins.fa     filtered protein FASTA (headers unchanged)
<out>.selection_report.txt      counts for QC / findings.md decision record
"""
from __future__ import annotations

import argparse
import gzip
import re
from collections import defaultdict

TAG_RE = {
    "gene_id": re.compile(r'gene_id "([^"]+)"'),
    "transcript_id": re.compile(r'transcript_id "([^"]+)"'),
    "gene_biotype": re.compile(r'(?:gene_biotype|gene_type) "([^"]+)"'),
    "transcript_biotype": re.compile(r'(?:transcript_biotype|transcript_type) "([^"]+)"'),
}

def strip_version(ens_id: str) -> str:
    return ens_id.split(".", 1)[0]


def parse_gtf(gtf_path: str):
    """Return {gene_id: {"mane": [enst...], "canonical": [enst...], "biotype": str}}"""
    genes = defaultdict(lambda: {"mane": [], "canonical": [], "biotype": None})
    n_transcript_lines = 0
    with open(gtf_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            n_transcript_lines += 1
            attrs = fields[8]
            gene_id = attrs_get(attrs, "gene_id")
            tx_id = attrs_get(attrs, "transcript_id")
            if not gene_id or not tx_id:
                continue
            g = genes[gene_id]
            g["biotype"] = attrs_get(attrs, "gene_biotype") or g["biotype"]
            tx_biotype = attrs_get(attrs, "transcript_biotype")
            is_pc = (tx_biotype == "protein_coding") or (
                tx_biotype is None and g["biotype"] == "protein_coding"
            )
            if not is_pc:
                continue
            if 'tag "MANE_Select"' in attrs:
                g["mane"].append(tx_id)
            if 'tag "Ensembl_canonical"' in attrs:
                g["canonical"].append(tx_id)
    print(f"GTF: parsed {n_transcript_lines} transcript lines, "
          f"{len(genes)} genes with protein-coding transcripts")
    return genes

def attrs_get(attrs: str, key: str) -> str | None:
    m = TAG_RE[key].search(attrs)
    return m.group(1) if m else None

def select_transcripts(genes):
    """One transcript per gene: MANE Select > Ensembl canonical."""
    selection = {}          # gene_id -> (enst, rule)
    unselected = []
    for gene_id, g in genes.items():
        if g["mane"]:
            selection[gene_id] = (g["mane"][0], "MANE_Select")
        elif g["canonical"]:
            selection[gene_id] = (g["canonical"][0], "Ensembl_canonical")
        else:
            unselected.append(gene_id)
    return selection, unselected

def filter_fasta(fasta_path: str, out_fasta: str, keep_enst: set[str]):
    """Keep FASTA records whose transcript:ENST (version-stripped) is selected."""
    tx_re = re.compile(r"transcript:(ENST[0-9]+)(?:\.\d+)?")
    gene_re = re.compile(r"gene:(ENSG[0-9]+)(?:\.\d+)?")
    n_in = n_out = 0
    kept_tx = set()
    with open(fasta_path, "r") as fin, open(out_fasta, "w") as fout:
        write = False
        for line in fin:
            if line.startswith(">"):
                n_in += 1
                tx = tx_re.search(line)
                tx_id = tx.group(1) if tx else None
                write = tx_id in keep_enst
                if write:
                    n_out += 1
                    kept_tx.add(tx_id)
            if write:
                fout.write(line)
    print(f"FASTA: {n_in} records in, {n_out} kept")
    return kept_tx

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", default="Homo_sapiens.GRCh38.116.gtf")
    ap.add_argument("--fasta", default="Homo_sapiens.GRCh38.pep.all.fa")
    ap.add_argument("--out", default="ensembl_canonical")
    args = ap.parse_args()

    genes = parse_gtf(args.gtf)
    selection, unselected = select_transcripts(genes)

    keep_enst = {strip_version(tx) for tx, _ in selection.values()}
    out_fasta = f"{args.out}_proteins.fa"
    kept_tx = filter_fasta(args.fasta, out_fasta, keep_enst)
    missing_in_fasta = keep_enst - kept_tx

    with open(f"{args.out}_selection.tsv", "w") as fh:
        fh.write("gene_id\ttranscript_id\trule\n")
        for gene_id, (tx, rule) in sorted(selection.items()):
            fh.write(f"{gene_id}\t{tx}\t{rule}\n")

    n_mane = sum(1 for _, r in selection.values() if r == "MANE_Select")
    n_can = sum(1 for _, r in selection.values() if r == "Ensembl_canonical")
    report = (
        f"protein-coding genes with >=1 pc transcript : {len(genes)}\n"
        f"selected (MANE Select)                       : {n_mane}\n"
        f"selected (Ensembl canonical)                 : {n_can}\n"
        f"genes without MANE/canonical tag (excluded)  : {len(unselected)}\n"
        f"selected transcripts missing from FASTA      : {len(missing_in_fasta)}\n"
    )
    print("\n" + report)
    with open(f"{args.out}_selection_report.txt", "w") as fh:
        fh.write(report)
        if missing_in_fasta:
            fh.write("missing transcripts:\n" + "\n".join(sorted(missing_in_fasta)) + "\n")

if __name__ == "__main__":
    main()
