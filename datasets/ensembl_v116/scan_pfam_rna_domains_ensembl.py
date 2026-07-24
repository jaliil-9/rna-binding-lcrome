#!/usr/bin/env python3
"""Scan selected RNA-related Pfam HMMs against Ensembl proteins and select one isoform/gene.

Default threshold reproduces the Gerstberger et al. 2014 E-value criterion:
full-sequence E-value < 0.01 AND independent-domain E-value < 0.01.
The input FASTA must contain all Ensembl protein isoforms, with headers carrying
'gene:ENSG...' fields, so isoforms can be collapsed afterwards.
"""
import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

DOMTBL_COLUMNS = [
    'target_name','target_accession','tlen','query_name','query_accession','qlen',
    'full_evalue','full_score','full_bias','domain_number','domain_of','c_evalue',
    'i_evalue','domain_score','domain_bias','hmm_from','hmm_to','ali_from','ali_to',
    'env_from','env_to','acc','description'
]

def run(cmd, **kwargs):
    print('+', ' '.join(map(str, cmd)))
    subprocess.run(list(map(str, cmd)), check=True, **kwargs)

def executable(hmmer_bin, name):
    p = Path(hmmer_bin).expanduser() / name if hmmer_bin else Path(name)
    if hmmer_bin and not p.exists():
        raise FileNotFoundError(f'Cannot find {p}')
    return str(p)

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def read_fasta(path):
    records, ident, desc, chunks = {}, None, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if ident:
                    records[ident] = {'sequence': ''.join(chunks), 'fasta_description': desc}
                desc = line[1:]
                ident = desc.split()[0]
                chunks = []
            else:
                chunks.append(line)
        if ident:
            records[ident] = {'sequence': ''.join(chunks), 'fasta_description': desc}
    return records

def ensembl_gene(description):
    m = re.search(r'(?:^|\s)gene:(ENSG\d+(?:\.\d+)?)', description)
    if not m:
        raise ValueError('No Ensembl gene:ENSG... field found in FASTA header: ' + description[:200])
    return m.group(1).split('.')[0]

def parse_domtbl(path, seqs, e_cutoff):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split(maxsplit=22)
            if len(fields) != 23:
                raise ValueError(f'Unexpected domtblout row ({len(fields)} fields): {line[:160]}')
            r = dict(zip(DOMTBL_COLUMNS, fields))
            for c in ['tlen','qlen','domain_number','domain_of','hmm_from','hmm_to','ali_from','ali_to','env_from','env_to']:
                r[c] = int(r[c])
            for c in ['full_evalue','full_score','full_bias','c_evalue','i_evalue','domain_score','domain_bias','acc']:
                r[c] = float(r[c])
            # target = Ensembl protein sequence; query = Pfam HMM in hmmscan output.
            if r['full_evalue'] >= e_cutoff or r['i_evalue'] >= e_cutoff:
                continue
            protein = r['query_name']
            if protein not in seqs:
                raise ValueError(f'{protein} from HMMER output not found in FASTA')
            r['human_sequence_id'] = protein
            r['ensembl_gene_id'] = ensembl_gene(seqs[protein]['fasta_description'])
            r['protein_length'] = len(seqs[protein]['sequence'])
            r['full_sequence'] = seqs[protein]['sequence']
            r['fasta_description'] = seqs[protein]['fasta_description']
            rows.append(r)
    return pd.DataFrame(rows)

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--hmmer-bin', required=True, help='Directory containing hmmscan, hmmfetch and hmmpress')
    ap.add_argument('--pfam-hmm', required=True, help='Pfam-A.hmm for the same Pfam release as the candidate list')
    ap.add_argument('--candidates-xlsx', required=True, help='Tiered candidate workbook from filter_rna_functional_pfam_domains.py')
    ap.add_argument('--candidates-sheet', default='retained_tiers')
    ap.add_argument('--ensembl-fasta', required=True, help='Ensembl pep.all FASTA: all isoforms, not pre-collapsed')
    ap.add_argument('--output', default='ensembl116_rna_pfam_hits_collapsed.xlsx')
    ap.add_argument('--workdir', default='hmmer_rna_pfam_work')
    ap.add_argument('--cpu', type=int, default=8)
    ap.add_argument('--evalue', type=float, default=0.01, help='Strict full-sequence and independent-domain E-value cutoff')
    ap.add_argument('--keep-workdir', action='store_true')
    args = ap.parse_args()

    work = Path(args.workdir); work.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_excel(args.candidates_xlsx, sheet_name=args.candidates_sheet)
    needed = {'accession','name','evidence_tier'}
    if not needed.issubset(candidates.columns):
        raise SystemExit(f'Candidate sheet must contain: {sorted(needed)}')
    candidates['accession'] = candidates['accession'].astype(str).str.replace(r'\.\d+$', '', regex=True)
    candidates = candidates.drop_duplicates('accession').copy()
    pfam_keys = work / 'selected_pfams.txt'

    # Map unversioned candidate accessions (e.g. PF00271) to the exact
    # versioned accession stored in this Pfam-A.hmm (e.g. PF00271.38).
    pfam_accession_map = {}
    with open(args.pfam_hmm, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("ACC"):
                versioned = line.split()[1]
                unversioned = versioned.split(".")[0]
                pfam_accession_map[unversioned] = versioned

    missing = sorted(set(candidates["accession"]) - set(pfam_accession_map))
    if missing:
        raise SystemExit(
            "Candidate accessions absent from this Pfam HMM file: "
            + ", ".join(missing)
        )

    pfam_keys.write_text(
        "\n".join(
            candidates["accession"].map(pfam_accession_map)
        ) + "\n"
    )

    hmmscan = executable(args.hmmer_bin, 'hmmscan')
    hmmfetch = executable(args.hmmer_bin, 'hmmfetch')
    hmmpress = executable(args.hmmer_bin, 'hmmpress')
    selected_hmm = work / 'selected_rna_pfams.hmm'
    # Pfam-A.hmm must have an .h3m index for hmmfetch. Create it once if absent.
    if not Path(str(args.pfam_hmm) + '.ssi').exists():
        run([hmmfetch, '--index', args.pfam_hmm])
    with open(selected_hmm, 'wb') as out:
        run([hmmfetch, '-f', args.pfam_hmm, pfam_keys], stdout=out)
    run([hmmpress, '-f', selected_hmm])

    domtbl = work / 'selected_rna_pfams.domtblout'
    tblout = work / 'selected_rna_pfams.tblout'
    log = work / 'hmmscan.log'
    if domtbl.exists() and domtbl.stat().st_size > 0:
        print(f"Reusing existing HMMER output: {domtbl}")
    else:
        with open(log, 'w') as out:
            run([hmmscan, '--cpu', args.cpu, '-E', args.evalue, '--domE', args.evalue,
                 '--domtblout', domtbl, '--tblout', tblout,
                 selected_hmm, args.ensembl_fasta], stdout=out)

    seqs = read_fasta(args.ensembl_fasta)
    hits = parse_domtbl(domtbl, seqs, args.evalue)
    if hits.empty:
        raise SystemExit('No hits passed the E-value threshold.')
    hits = hits.merge(candidates[['accession','name','evidence_tier']], how='left',
                      left_on='target_accession', right_on='accession')
    hits = hits.rename(columns={'query_name':'pfam_name','target_accession':'pfam_accession',
                                'name':'candidate_name'})
    # Each passing domain instance counts once. Select isoform with most domains; then longest.
    counts = hits.groupby('human_sequence_id').size().rename('number_of_domains').reset_index()
    proteins = pd.DataFrame([{'human_sequence_id': k, 'ensembl_gene_id': ensembl_gene(v['fasta_description']),
                              'protein_length': len(v['sequence']), 'full_sequence': v['sequence'],
                              'fasta_description': v['fasta_description']} for k,v in seqs.items()])
    scored = proteins.merge(counts, how='left', on='human_sequence_id')
    scored['number_of_domains'] = scored['number_of_domains'].fillna(0).astype(int)
    winners = (scored.sort_values(['ensembl_gene_id','number_of_domains','protein_length','human_sequence_id'],
                                  ascending=[True,False,False,True])
                     .drop_duplicates('ensembl_gene_id', keep='first'))
    winners['selection_rule'] = 'maximum passing Pfam-domain instances; tie: longest protein; final tie: lexical ENSP ID'
    selected_ids = set(winners.loc[winners.number_of_domains.gt(0), 'human_sequence_id'])
    selected_hits = hits[hits.human_sequence_id.isin(selected_ids)].copy()
    domain_summary = (selected_hits.groupby(['pfam_accession','pfam_name','candidate_name','evidence_tier'])
                      .agg(human_proteins=('human_sequence_id','nunique'), domain_instances=('human_sequence_id','size'),
                           best_independent_domain_evalue=('i_evalue','min'), best_domain_score=('domain_score','max'))
                      .reset_index().sort_values('human_proteins', ascending=False))
    metadata = pd.DataFrame({'field':['run_utc','pfam_hmm','pfam_hmm_sha256','candidates_xlsx','candidates_sha256','ensembl_fasta','ensembl_fasta_sha256','evalue_cutoff','hmmscan_command','selection_rule'],
     'value':[datetime.now(timezone.utc).isoformat(),str(args.pfam_hmm),sha256(args.pfam_hmm),str(args.candidates_xlsx),sha256(args.candidates_xlsx),str(args.ensembl_fasta),sha256(args.ensembl_fasta),str(args.evalue),f'{hmmscan} -E {args.evalue} --domE {args.evalue}', winners.selection_rule.iloc[0]]})
    with pd.ExcelWriter(args.output, engine='openpyxl') as xw:
        selected_hits.sort_values(['ensembl_gene_id','human_sequence_id','pfam_accession','ali_from']).to_excel(xw, sheet_name='selected_isoform_hits', index=False)
        hits.sort_values(['ensembl_gene_id','human_sequence_id','pfam_accession','ali_from']).to_excel(xw, sheet_name='all_isoform_hits', index=False)
        winners.sort_values('ensembl_gene_id').to_excel(xw, sheet_name='one_isoform_per_gene', index=False)
        domain_summary.to_excel(xw, sheet_name='domain_summary', index=False)
        metadata.to_excel(xw, sheet_name='run_metadata', index=False)
    print(f'Wrote {args.output}: {len(hits)} passing domain hits; {len(selected_hits)} retained after collapse; {len(selected_ids)} RBP-domain-positive genes.')
    print(f'Raw HMMER files: {domtbl}, {tblout}, {log}')
    if not args.keep_workdir:
        print('Work directory retained intentionally as an audit trail; delete it manually after validation.')

if __name__ == '__main__':
    main()
