# RNA-binding protein superclass assignment — Stage 1 report

**Project:** Mapping RNA binders / functional classification of LCR-associated RBPs  
**Stage:** 1 — Gerstberger S3 transfer plus UniProt annotation fallback  
**Run date:** 2026-07-29  
**Input proteins:** 2,098  
**Output:** `rbp_rna_classification.xlsx`

## Purpose

This stage assigns every protein in the combined human RBP candidate set to one primary RNA-target superclass. The goal is a transparent, protein-level starting dataset for downstream LCR analyses, not a claim that every candidate is a directly experimentally validated RNA binder.

The classification follows the target framework in Gerstberger, Hafner and Tuschl (2014), which manually curated 1,542 human RBPs/RNP components and assigned a predominant RNA-target category from literature evidence. The original study treats ribosomal proteins separately from pre-rRNA factors and retains an `unknown` class for RBPs whose natural RNA targets were unresolved.

## Input data

### Candidate protein set

The input workbook was:

```text
combined_rbp_pfam38_rbpdb_modomics_uniprot.xlsx
```

It contained 2,098 protein-level rows assembled from the project’s combined RBP evidence sources. The classifier used the following metadata fields when Gerstberger S3 did not provide a match:

| Input field | Use in fallback annotation |
|---|---|
| `uniprot_accession` | Primary protein identifier and output identifier |
| `uniprot_entry_name` | Gene-symbol fallback for S3 matching |
| `uniprot_protein_name` | Protein-name keyword evidence |
| `go_biological_process` | Biological-process evidence |
| `go_molecular_function` | Molecular-function evidence |
| `subcellular_location` | Supporting UniProt context |
| `uniprot_domains` | Supporting UniProt domain/context evidence |
| `ensembl_protein_id` | Additional possible S3 protein-ID match |

### Curated reference

The reference workbook was:

```text
41576_2014_BFnrg3813_MOESM25_ESM.xls
```

The classifier reads its `RBP table` sheet. Four S3 fields were used:

| S3 field | Use |
|---|---|
| `protein id` | Direct curated match key |
| `gene name` | Fallback curated match key |
| `consensus RNA target` | Assigned primary RNA-target superclass |
| `putative RNA target` | Secondary class when distinct from the consensus target |
| `supporting evidence (# pubmed ID)` | Preserved in the evidence text |

No new manual target assignments were added at this stage.

## Controlled vocabulary

Every row receives one of the following primary labels:

| Label | Interpretation |
|---|---|
| `ribosomal protein` | Structural ribosomal proteins; intentionally distinct from pre-rRNA factors |
| `mRNA` | mRNA/pre-mRNA processing, splicing, export, stability, localization, translation-related pathways, or mRNA binding |
| `tRNA` | tRNA biogenesis, charging, modification, transport, splicing, or binding |
| `pre-rRNA` | rRNA processing, ribosome biogenesis, and pre-rRNA-associated pathways |
| `snRNA` | snRNA/snRNP pathways |
| `snoRNA` | snoRNA/snoRNP/scaRNA pathways |
| `ncRNA` | Other non-coding RNA classes, including miRNA, piRNA, lncRNA, 7SK, 7SL, Y RNA, vault RNA, telomerase RNA, RNase P and RNase MRP |
| `diverse` | Broad/non-specific RNA turnover or binding to diverse RNA classes, including the RNA exosome |
| `unknown` | RBP status or RNA-related evidence exists, but no natural target class is assigned |

The S3 source label `ribosome` was explicitly normalized to `ribosomal protein`. This correction recovered the expected curated ribosomal-protein group.

## Assignment procedure

### 1. Curated S3 transfer

For each input row, the script first attempts to find a Gerstberger S3 record using:

1. `uniprot_accession` against S3 `protein id`
2. `ensembl_protein_id` against S3 `protein id`, where available
3. Gene symbol derived from `uniprot_entry_name` against S3 `gene name`

On a match:

- `consensus RNA target` becomes `rna_primary_class`.
- A distinct valid `putative RNA target` becomes `rna_secondary_class`.
- Supporting S3 PubMed/evidence text is copied into `evidence_matched`.
- The row receives `confidence = very high` and `evidence_source = literature based (Gerstberger S3)`.

`very high` therefore denotes a direct transfer from the manually curated 2014 S3 table. It does **not** mean that the assigned RNA target is always experimentally resolved: an S3 record with consensus `unknown` is still a high-provenance curated unknown.

### 2. UniProt fallback

Rows without an S3 match are assigned using target-specific keywords in protein names, GO biological-process and molecular-function annotations, subcellular-location comments, and domain annotations.

The keyword rules deliberately do **not** use generic RNA-binding features alone (for example RRM, KH, helicase, zinc finger, G-patch). Those features support RNA-related candidacy but do not uniquely identify the RNA substrate class.

Examples of target-specific fallback evidence include:

- `mRNA`: messenger RNA, pre-mRNA, splicing, spliceosome, hnRNP, polyadenylation, mRNA export/localization/decay.
- `tRNA`: transfer RNA, aminoacyl-tRNA, anticodon, tRNA processing/modification/splicing.
- `pre-rRNA`: pre-rRNA, rRNA processing/maturation/modification, ribosome biogenesis, 18S/28S/5.8S/5S rRNA.
- `snRNA`: snRNA, snRNP, small nuclear ribonucleoprotein, U1/U2/U4/U5/U6.
- `snoRNA`: snoRNA, snoRNP, scaRNA, box C/D, box H/ACA.
- `ncRNA`: miRNA, piRNA, PIWI, lncRNA, 7SK, 7SL, Y RNA, vault RNA, telomerase RNA, RNase P/MRP.
- `diverse`: RNA exosome, general RNA turnover, RNA surveillance/degradation, RNA:DNA hybrid.

If multiple classes match, the primary class is selected by the number of independent metadata fields supporting it, followed by the number of matching terms. Other matched classes are retained as semicolon-separated secondary labels.

### 3. Confidence rules

| Confidence | Rule |
|---|---|
| `very high` | Direct Gerstberger S3 match |
| `high` | At least two distinct UniProt metadata fields support the primary class; or a target-specific protein name alone supports it |
| `medium` | The class is supported only by GO, location, and/or domain metadata |
| `low` | Weak protein-name-only assignment not meeting the target-specific naming criterion |
| `unknown` | No class-specific evidence found |

## Results

### Primary superclass counts

| Primary class | Proteins | Share of 2,098 |
|---|---:|---:|
| mRNA | 835 | 39.8% |
| unknown | 432 | 20.6% |
| tRNA | 187 | 8.9% |
| ribosomal protein | 181 | 8.6% |
| pre-rRNA | 131 | 6.2% |
| ncRNA | 131 | 6.2% |
| snRNA | 108 | 5.1% |
| diverse | 51 | 2.4% |
| snoRNA | 42 | 2.0% |
| **Total** | **2,098** | **100.0%** |

### Confidence counts

| Confidence | Proteins | Share of 2,098 |
|---|---:|---:|
| very high | 1,632 | 77.8% |
| high | 89 | 4.2% |
| medium | 44 | 2.1% |
| low | 10 | 0.5% |
| unknown | 323 | 15.4% |
| **Total** | **2,098** | **100.0%** |

### Key validation observation

After mapping the S3 label `ribosome` to the project label `ribosomal protein`, the run recovered **169 curated ribosomal proteins** as `very high` confidence. This closely matches the 169 ribosomal proteins reported in the Gerstberger census. The final `ribosomal protein` total is 181 because the remaining 12 proteins are additional UniProt-derived assignments.

## Interpretation

The output is suitable as a first-pass functional map of the combined candidate set:

- Most rows are anchored to a directly matched curated S3 record.
- `mRNA` is the largest class, as expected for an RBP resource enriched in mRNA-associated proteins.
- `unknown` should be retained, not forced into another class. It has two distinct meanings: curated S3 proteins with unresolved natural targets, and non-S3 candidates with no sufficiently specific UniProt target annotation.
- The classes beyond the 2014 census are expected to contain proteins contributed by the project’s other sources and current UniProt annotation.

## Limitations

1. **Historical reference.** Gerstberger S3 reflects a 2014 curation and its identifiers, evidence, and target calls may not capture later studies.
2. **One primary class is a simplification.** Many proteins bind multiple RNA classes. The consensus target is a predominant pathway label, not a complete binding profile.
3. **Fallback rules are annotation-driven.** UniProt/GO terms can be incomplete, broad, inferred electronically, or reflect pathway membership rather than direct RNA contact.
4. **`diverse` needs careful interpretation.** Generic nuclease or RNA-turnover wording can obscure substrate specificity; this class should be manually reviewed before biological conclusions are drawn.
5. **Candidate-set composition matters.** These counts describe the 2,098-protein combined input, not the whole human proteome or a new exhaustive RBP census.

