# RNA-binding Pfam pilot: execution report

## Scope
Human Ensembl GRCh38 protein translations were scanned with selected Pfam 38.2 HMMs to build a domain-supported candidate protein set for later curation and LCR/SEG analysis.

## Pipeline
Candidate models -> HMMER hmmfetch/hmmpress -> hmmscan against Homo_sapiens.GRCh38.pep.all.fa -> parse domain calls -> one selected isoform per ENSG -> comparison with 2014 census.

## Run statistics
| Metric | Result |
|---|---:|
| Scanned Pfam HMMs | 172 |
| E-value thresholds | full E <= 0.01; independent domain E <= 0.01 |
| Passing isoform-level domain calls | 15,570 |
| Hit-positive isoforms | 10,969 |
| Selected representative protein records | 965 |
| Distinct selected gene symbols | 859 |
| Retained domain instances | 1,275 |
| Hit-producing HMMs | 140 |
| Models with no selected-isoform hit | 32 |

## 2014 comparison
| Measure | Count |
|---|---:|
| 2014 Gerstberger et al. census entries | 1,542 |
| Current distinct protein symbols | 859 |
| Shared symbols | 299 |
| Current-only symbols | 560 |
| 2014-only symbols | 1,243 |
| Nonredundant union | 2,102 |

These are symbol-level, deduplicated counts; current output has 965 ENSP records because several records map to repeated symbols.

## Key technical corrections
1. Candidate combination initially stopped before writing because validation raised SystemExit. A non-blocking writer was provided.
2. hmmfetch failed because whitespace-containing HMM names were used as keys (e.g., Influenza...). The script was changed to map unversioned PF accessions to exact versioned accessions in Pfam-A.hmm.
3. hmmscan completed but parser treated the profile target as the protein query. The orientation was corrected: hmmscan target = Pfam HMM; query = human protein. Existing domtblout was reused instead of rescanning.
4. Excel warned that a sequence exceeded its 32,767-character per-cell limit. FASTA, not Excel sequence cells, remains the authoritative sequence input.
5. The output hit tables have an unresolved candidate-metadata join defect: candidate accession/name/tier columns are blank. Raw domtblout is preserved and supports rebuilding correct joins.

## Current status
The result is a broad, unfiltered domain-supported candidate screen, not a final RBP census. The 54 GO-derived entries include viral, bacterial and indirect RNA-function families; PF15337 (Vasculin-like 1) is clearly DNA/promoter-associated and should be excluded from an RNA-domain list. The 299 shared 2014 proteins are validation/benchmark overlap, not evidence that historical domains were reused. The 560 current-only symbols are candidates requiring review.

## Deliverables
- ensembl116_rna_pfam_hits_collapsed.xlsx: raw/collapsed HMMER workbook
- hmmer_rna_pfam_work/selected_rna_pfams.domtblout: audit-grade raw HMMER domains
- review_and_rebuild_pfam_hits.py: model-review/rebuild helper
- current_pfam_vs_2014_census_overlap.csv: earlier mass-spec-sheet overlap
- current_vs_2014_RBP_census_new_symbols.csv: 560 current-only symbols
- current_vs_2014_RBP_census_overlap_symbols.csv: 299 shared census symbols

## Next actions
1. Review the 172 models using a human-relevant RNA-function allow-list; exclude virus/phage/toxin/CRISPR and non-RNA models.
2. Rebuild hit annotations from raw domtblout with correct Pfam accession joining.
3. Re-collapse reviewed hits to one declared representative protein per gene using Ensembl mapping.
4. Join UniProt/Swiss-Prot and literature evidence, with separate target-RNA labels and direct-versus-indirect evidence.
5. Freeze final protein FASTA, then run SEG and quantify LCR/domain overlap.
