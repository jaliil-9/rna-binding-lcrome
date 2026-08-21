# Background Proteome Construction and RBP Label Reconciliation

**Status:** background construction complete; ready for six-method LCR detection run.
**Date:** August 17, 2026.

---

## 1. Objective

This document records the construction of the shared background proteome and the RBP label layer for runs 1–2, including the QC that motivated each design decision.

## 2. Design decision: one protein per gene

The full Ensembl translation set (~300k records, `pep.all`) was rejected as a background because:

- **Statistical unit integrity.** Isoforms of one gene share most of their sequence; the same LCR would be detected multiple times as near-identical intervals, introducing pseudo-replication at the protein level — the level designated as the robustness criterion in Phase 2.4.1. RBPs are enriched for long, multi-isoform genes, so the bias would be directional (inflating apparent RBP LCR incidence).
- **Label compatibility.** RBP labels (master table, RBP2GO, RBPbase) are UniProt-centric, one canonical sequence per protein.
- **Field convention.** Reference sets in the comparable literature are one-protein-per-gene (e.g., the RBP2GO refinement study's whole-proteome reference, 20,741 proteins).

## 3. Route evaluation: Ensembl canonical vs. UniProt reference proteome

### 3.1 Ensembl r116 canonical set (constructed, then rejected)

A MANE Select ∪ Ensembl Canonical set was built from the r116 GTF (`chr_patch_hapl_scaff`) and `pep.all.fa`:

| QC metric | Value |
|---|---|
| Genes with ≥1 protein-coding transcript | 78,941 |
| Selected via MANE Select | 19,229 |
| Selected via Ensembl Canonical (fallback) | 544 |
| Excluded (no representative tag; patch/alt-loci) | 59,168 |
| Selected transcripts missing from FASTA | 0 |
| Transcripts mapped to UniProt accessions (r116 xref TSV) | 19,695 / 19,773 (99.6%) |

Reconciliation against the master RBP table (ENSP-identity and sequence-identity based) showed this set is **not sequence-compatible with the Phase 2.4.1 RBP arm**:

| Reconciliation bucket | Rows | Share |
|---|---|---|
| Exact translation present (master ENSP in set) | 1,061 | 55.1% |
| — of which sequence-identical to `uniprot_sequence` | 1,028 | — |
| Gene present, different isoform | 703 | 36.5% |
| Gene absent | 161 | 8.4% |

Ensembl's MANE-driven canonical choice diverges from the UniProt canonical isoform for ~45% of the RBP set. Mixing sequence universes would make RBP-side LCR coordinates non-transferable for those genes.

**Decision:** rejected. The UniProt reviewed reference proteome was adopted instead — native accession join, sequences identical to the RBP arm by construction, and direct comparability with the RBP2GO reference sets.

### 3.2 UniProt reference proteome (UP000005640, reviewed) — adopted

Download: UniProtKB REST stream, `query=(proteome:UP000005640) AND (reviewed:true)`, FASTA format → 20,416 records.

Reconciliation against the master table (1,924 distinct non-null accessions; 144 rows lack an accession):

| Metric | Value |
|---|---|
| Accession present in background | 1,463 (76.0%) |
| — sequence-identical to `uniprot_sequence` | 1,463 / 1,463 (100%) |
| Absent accessions | 462 |

The 462 absent accessions (605 master rows) decomposed as: 461 unreviewed (TrEMBL) accessions of Pfam27/Pfam38-census vintage, plus 144 rows with no accession at all (collapsing to one null bucket). Sources: pfam27 (400 rows), pfam38 (200), pfam27;pfam38 (3), rbpdb (2).

## 4. Gene-level resolution of the absent rows

Each absent row was resolved to its gene via the master `ensembl_protein_id` and the Ensembl r116 xref (ENSP→ENSG→current UniProt accessions):

| Bucket | Rows | Interpretation / action |
|---|---|---|
| `gene_already_rbp_labeled` | 115 | Legacy duplicates; gene already labeled via another master row. No action. |
| `gene_in_background_unlabeled` | 201 (199 distinct genes) | **Label leak**: gene present in background under a current reviewed accession but would have been labeled non-RBP. Fixed by adding the covering accessions to the RBP label set (`gene_patch`). |
| `gene_not_in_background` | 12 | Gene absent from the reviewed proteome. 11 sequences appended to the background directly from the master table (`master_appended`; the exact sequences used for RBP-side LCR detection). 1 row had no usable sequence — excluded, documented. |
| `no_ensp_mapping` | 277 | Not resolvable locally (see §6). |

The 199-gene patch matters: without it, ~10% of the RBP set would have been false negatives inside the negative set.

## 5. Final deliverables

| File | Content |
|---|---|
| `finalized_background.fa` | 20,427 records: 20,416 UniProt reviewed canonical proteins + 11 appended master sequences (`>tr|<acc>|... master-appended`) |
| `finalized_background_labels.tsv` | Per-protein `is_rbp` + `label_source` ∈ {`master_direct`, `gene_patch`, `master_appended`} |
| `finalized_label_patch.tsv` | The 199 patched genes: master accession → covering background accession |
| `finalized_unresolved.tsv` | The 277 unresolved rows for optional follow-up |

**Label summary:** 1,669 RBP-labeled proteins (8.2% of the background): 1,463 direct accession matches + 199 patch accessions + 11 appended (the label file is authoritative; overlaps between sources are deduplicated).

## 6. Known limitations and open items

- **277 unresolved master rows.** 134 have both an accession and a sequence and are recoverable via the UniProt ID-mapping service (handles retired/merged accessions); ~143 rows have neither accession, current ENSP, nor sequence and are likely unrecoverable without returning to the original Pfam27/Pfam38 source lists. Decision deferred; not blocking. Maximum potential impact bounded at ~7% of master rows.
- **1 gene** of the 12 background-absent genes had no usable master sequence; excluded.
- **Ensembl r116 divergence** (§3.1) is retained as the documented justification for the UniProt route; the Ensembl-derived files are kept for provenance but are not analysis inputs.
- **Rung 3** (consensus non-RBP census vs. RBP2GO/RBPbase) requires an additional label join against those catalogs; not yet performed.

## 7. Instructions for the detection run

1. Run all six method configurations (CAST, SEG, SEG_intermediate, fLPS, LCRFinder, AlcoR) on `finalized_background.fa` with parameters identical to the RBP-side run.
2. FASTA headers are UniProt-style (`>sp|<ACC>|<NAME> ...`). When building the combined LCR table, parse the accession as the middle `|`-delimited field and construct `protein_id` accession-first (`<ACC>|...`), so the existing `extract_accession()` join logic works unchanged.
3. Run `pc_analysis_annotation.py` once on the **combined RBP ∪ background LCR table** — the cation–π signature threshold is a dataset-relative quantile (`cooc_top_quantile = 0.75`), so both arms must be annotated in a single run to share one threshold and one signature taxonomy.
4. Labels join via `finalized_background_labels.tsv` (`uniprot_accession` → `is_rbp`).

## 8. Scripts produced (this phase)

| Script | Role |
|---|---|
| `build_canonical_proteome.py` | Ensembl route: MANE/canonical selection + FASTA filtering |
| `add_uniprot_ids_v3.py` | Ensembl route: UniProt xref annotation (content-validated column detection; Swiss-Prot preference) |
| `length_concordance_qc.py` | Ensembl route: accession-level length concordance vs. master |
| `reconcile_master_background_v2.py` | ENSP/accession/sequence-level reconciliation (both routes) |
| `diagnose_absent.py` | Classification of absent accessions (reviewed status, source) |
| `resolve_absent_genes.py` | Gene-level bucketing of absent rows |
| `finalize_background.py` | Label patch + sequence append + final label file |

## 9. Next steps

1. Six-method LCR detection on `finalized_background.fa`.
2. Combined feature extraction + behavior annotation (RBP ∪ background).
3. `phase2_5_compare.py`: RBP vs. background comparison (incidence, architecture, signature specificity/occupancy, charge landscape) with length-confounder handling (stratification, matched subset, logistic adjustment).
4. Optional: UniProt ID-mapping pass over `finalized_unresolved.tsv` (134 recoverable rows).
5. Rung 3: extend labels with RBP2GO 2.0 and RBPbase memberships.
