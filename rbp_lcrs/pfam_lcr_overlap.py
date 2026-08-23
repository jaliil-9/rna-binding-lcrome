"""Annotate LCRs with Pfam domain overlaps using hmmscan."""

import argparse
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

COLS = [
    'pfam_name', 'pfam_accession', 'hmm_length', 'protein_id',
    'query_accession', 'protein_length', 'full_evalue', 'full_score',
    'full_bias', 'domain_number', 'domain_of', 'c_evalue', 'i_evalue',
    'domain_score', 'domain_bias', 'hmm_start', 'hmm_end', 'domain_start',
    'domain_end', 'env_start', 'env_end', 'accuracy', 'description',
]

AA = ''


def norm(x):
    return re.sub(r'[^a-z0-9]+', '', str(x).lower())


def pf(x):
    return re.sub(r'\.\d+$', '', str(x).strip().upper())


def pid(x):
    x = str(x).strip().split()[0]
    p = x.split('|')
    return (p[1] if len(p) > 1 and p[0].lower() in {'sp', 'tr'} else p[0]).upper()


def sheet(path, want):
    x = pd.ExcelFile(path)
    for s in x.sheet_names:
        if norm(s) == norm(want):
            return s
    raise ValueError(f'{want} not found; sheets: {x.sheet_names}')


def col(df, want):
    d = {norm(c): c for c in df.columns}
    if norm(want) not in d:
        raise ValueError(f'{want} missing; columns: {list(df.columns)}')
    return d[norm(want)]


def selected(tiers, status):
    a = pd.read_excel(tiers, sheet_name=sheet(tiers, 'retained_tiers'))
    b = pd.read_excel(status, sheet_name=sheet(status, 'Sheet1'))
    left = {pf(x) for x in a[col(a, 'accession')].dropna()}
    right = {
        pf(x)
        for x, s in zip(b[col(b, 'accession')], b[col(b, 'status')])
        if str(s).strip().lower() in {'found', 'renamed'}
    }
    keep = sorted(left | right)
    if not keep:
        raise ValueError('No Pfam accessions were retained from the input files.')
    return keep


def scan(fasta, hmm, dom, cpu):
    cmd = [
        'hmmscan', '--noali', '--cpu', str(cpu),
        '--domtblout', str(dom), str(hmm), str(fasta),
    ]
    print('Running:', ' '.join(cmd))
    with open(dom.with_suffix('.log'), 'w', encoding='utf8') as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)


def domains(dom, keep):
    rows = []
    with open(dom, encoding='utf8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            z = line.rstrip().split(maxsplit=22)
            if len(z) == 22:
                z.append('')
            if len(z) >= 23:
                rows.append(z)

    d = pd.DataFrame(rows, columns=COLS)
    if d.empty:
        return pd.DataFrame(columns=[
            'protein_id', 'pfam_accession', 'pfam_name', 'domain_start',
            'domain_end', 'i_evalue', 'domain_score', 'hmm_coverage',
        ])

    nums = [
        'hmm_length', 'i_evalue', 'domain_score', 'hmm_start', 'hmm_end',
        'domain_start', 'domain_end',
    ]
    d[nums] = d[nums].apply(pd.to_numeric, errors='coerce')
    d['protein_id'] = d.protein_id.map(pid)
    d['pfam_accession'] = d.pfam_accession.map(pf)
    d['hmm_coverage'] = (d.hmm_end - d.hmm_start + 1) / d.hmm_length
    d = d[
        d.pfam_accession.isin(keep)
        & (d.i_evalue <= 1e-3)
        & (d.hmm_coverage >= 0.30)
    ].copy()

    return (
        d[['protein_id', 'pfam_accession', 'pfam_name', 'domain_start',
           'domain_end', 'i_evalue', 'domain_score', 'hmm_coverage']]
        .drop_duplicates()
        .sort_values(['protein_id', 'domain_start'])
    )


def lcrs(path):
    s = sheet(path, 'all_results')
    d = pd.read_excel(path, sheet_name=s)
    c = {norm(x): x for x in d.columns}
    for x in ['protein_id', 'method', 'start', 'end']:
        if norm(x) not in c:
            raise ValueError(f'{x} missing in {s}: {list(d.columns)}')

    get = lambda x: d[c[norm(x)]] if norm(x) in c else ''

    o = pd.DataFrame({
        'source_method': get('source_method'),
        'protein_id': get('protein_id'),
        'header': get('header'),
        'method': get('method'),
        'lcr_start': get('start'),
        'lcr_end': get('end'),
        'sequence': get('sequence'),
        'description': get('description'),
    })
    o.protein_id = o.protein_id.map(pid) # type: ignore
    o.method = o.method.astype(str).str.strip()  # type: ignore
    o[['lcr_start', 'lcr_end']] = (
        o[['lcr_start', 'lcr_end']].apply(pd.to_numeric, errors='coerce')
    )
    o = o.dropna(subset=['protein_id', 'method', 'lcr_start', 'lcr_end']).copy()
    o[['lcr_start', 'lcr_end']] = o[['lcr_start', 'lcr_end']].astype(int)
    o['lcr_length'] = o.lcr_end - o.lcr_start + 1
    return o.drop_duplicates(['method', 'protein_id', 'lcr_start', 'lcr_end'])


def overlap(l, d):
    x = l.merge(d, on='protein_id', how='inner')
    x['overlap_start'] = x[['lcr_start', 'domain_start']].max(axis=1)
    x['overlap_end'] = x[['lcr_end', 'domain_end']].min(axis=1)
    x['overlap_length'] = (x.overlap_end - x.overlap_start + 1).clip(lower=0)
    x = x[x.overlap_length > 0].copy()
    x['lcr_domain_overlap_fraction'] = x.overlap_length / x.lcr_length
    x['overlap_class'] = pd.cut(
        x.lcr_domain_overlap_fraction,
        [0, 0.30, 0.70, 1],
        labels=['<=0.3', '0.3<x<0.7', '>=0.7'],
        include_lowest=True,
    )
    return x.sort_values(['method', 'protein_id', 'lcr_start', 'domain_start'])


def summary(all_lcr, ov):
    total = all_lcr.groupby(['method', 'protein_id']).size().rename('total_lcrs')
    hit = (
        ov.drop_duplicates(['method', 'protein_id', 'lcr_start', 'lcr_end'])
        .groupby(['method', 'protein_id'])
        .size()
        .rename('lcrs_overlapping_pfam')
    )
    out = pd.concat([total, hit], axis=1).fillna(0).reset_index()
    out[['total_lcrs', 'lcrs_overlapping_pfam']] = (
        out[['total_lcrs', 'lcrs_overlapping_pfam']].astype(int)
    )
    out['pct_lcrs_overlapping_pfam'] = (
        100 * out.lcrs_overlapping_pfam / out.total_lcrs
    ).round(2)
    return out


def safe(x):
    return re.sub(r'[\\/*?:\[\]]', '_', str(x))[:31]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--tiers-xlsx', required=True)
    p.add_argument('--status-xlsx', required=True)
    p.add_argument('--proteins-fasta', required=True)
    p.add_argument('--pfam-hmm', required=True)
    p.add_argument('--lcr-xlsx', required=True)
    p.add_argument('--out', default='pfam_lcr_overlap.xlsx')
    p.add_argument('--cpu', type=int, default=8)
    p.add_argument('--domtblout', default='pfam_scan.domtblout')
    p.add_argument('--skip-scan', action='store_true')
    a = p.parse_args()

    keep = selected(a.tiers_xlsx, a.status_xlsx)
    print(f'Retained Pfam accessions: {len(keep)}')

    dom = Path(a.domtblout)
    if not a.skip_scan:
        scan(a.proteins_fasta, a.pfam_hmm, dom, a.cpu)
    if not dom.exists():
        raise FileNotFoundError(f'Domain table not found: {dom}')

    d = domains(dom, set(keep))
    l = lcrs(a.lcr_xlsx)
    o = overlap(l, d)
    s = summary(l, o)

    with pd.ExcelWriter(a.out, engine='openpyxl') as w:
        pd.DataFrame({'pfam_accession': keep}).to_excel(
            w, sheet_name='Retained_Pfam', index=False
        )
        d.to_excel(w, sheet_name='Pfam_hits', index=False)
        o.to_excel(w, sheet_name='LCR_Pfam_overlaps', index=False)
        s.to_excel(w, sheet_name='Summary', index=False)
        for m, z in o.groupby('method', sort=True):
            z.to_excel(w, sheet_name=safe(m), index=False)

    print(f'Wrote {a.out}; retained hits={len(d):,}; '
          f'LCR-domain overlaps={len(o):,}')


if __name__ == '__main__':
    main()