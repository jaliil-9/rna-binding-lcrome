# Phase 2.4.1 — Quantitative and Diversity Analysis of LCRs: Context, Design, and Findings

**Status:** analysis complete; all six methods validated at LCR and protein levels; continuous-metric and machine-checked concordance layers complete; findings finalized.
**Date:** August 2026.

---

## 1. Project context

This work is part of the PhD research program on **low-complexity regions (LCRs) in human RNA-binding proteins (RBPs)**, within the broader benchmarking effort across LCR detection tools: **CAST, SEG, fLPS, LCRFinder, and AlcoR**, with SEG additionally run in an intermediate-stringency variant (`SEG_intermediate`). Earlier phases produced per-method LCR calls, a rule-based physicochemical behavior annotation layer (signatures such as `basic_enriched_region`, `arginine_rich_rna_contact_region`, `sr_rs_repeat_region`), RNA-target superclass assignments, and domain-position classifications.

## 2. Phase 2.4.1 objective

A first **quantitative and diversity analysis** of LCR annotations, examining how LCR characteristics vary by **RNA-target superclass** (mRNA, tRNA, pre-rRNA, snRNA, snoRNA, ncRNA, ribosomal protein, diverse, unknown), stratified by detection method, across five analysis levels:

1. Domain-position enrichment relative to RNA classes
2. Physicochemical-property enrichment relative to RNA classes
3. Sequence-level properties (length, amino-acid enrichment)
4. LCR multiplicity per protein (1 vs. multiple), per method
5. Overall residue coverage and physicochemical-property fractions

## 3. Input data

| Input | Role | Key columns |
|---|---|---|
| `lcr_annotations.xlsx` (Sheet1, 20,954 rows) | LCR calls with annotations | `protein_id` (pipe-delimited composite), `source_method`, `start`, `end`, `length`, `sequence`, `rna_target_superclass`, `domain_position_class`, `primary_physicochemical_annotation` |
| `lcr_features.xlsx` (3 sheets: composition, distribution, co_occurrence; 20,954 rows each) | Per-LCR measurements | `frac_A`…`frac_Y`, property fractions, `fcr`, `ncpr`, `frac_GS`, RG/SR/GS repeat columns, run/gap distribution statistics, `cooc_positive_aromatic` |
| Master RBP table (`Combined` sheet) | Protein lengths only | `uniprot_accession`, `uniprot_length` |

## 4. Measured metrics and units of measurement

### 4.1 Identifiers and coordinates (per LCR)

| Metric | Unit | Description |
|---|---|---|
| `uniprot_accession` | — (identifier) | UniProt accession parsed from the composite protein ID (text before the first `\|`) |
| `method` | — (categorical) | Detection tool; canonical labels: CAST, SEG, SEG_intermediate, FLPS, LCRFinder, AlcoR |
| `start`, `end` | amino-acid position (1-based) | LCR interval boundaries |
| `length` | amino acids (aa) | LCR length |
| `rna_target_superclass` | — (categorical) | RNA-target superclass of the parent protein (9 classes) |
| `domain_position_class` | — (categorical) | `domain_intrinsic`, `domain_edge`, `domain_adjacent`, `interdomain_linker`, `distal_terminal`, `unclassified_no_pfam` |
| `primary_physicochemical_annotation` | — (categorical) | Primary physicochemical behavior signature |

### 4.2 Composition metrics (per LCR)

All fraction metrics are **unitless ratios in [0, 1]**.

| Metric | Unit | Description |
|---|---|---|
| `frac_A` … `frac_Y` | fraction (0–1) | Per-residue amino-acid fractions (20 columns) |
| `frac_polar` | fraction (0–1) | Polar residues (R, N, D, Q, E, H, K, S, T, W, Y) |
| `frac_hydrophobic` | fraction (0–1) | Hydrophobic residues (A, C, G, H, I, K, L, M, F, T, W, V, Y) |
| `frac_strong_hydro` | fraction (0–1) | Strong hydrophobes (V, I, L, M, F, W, Y) |
| `frac_aromatic` | fraction (0–1) | Aromatic residues (F, Y, W, H) |
| `frac_disorder` | fraction (0–1) | Disorder-promoting residues, TOP-IDP set (A, R, G, Q, S, P, E, K) |
| `frac_positive` (f+) | fraction (0–1) | Positively charged at pH 7.4 (R, K) |
| `frac_negative` (f−) | fraction (0–1) | Negatively charged at pH 7.4 (D, E) |
| `fcr` | fraction (0–1) | Fraction of charged residues, f+ + f− (Das & Pappu 2013) |
| `ncpr` | unitless (−1 to +1) | Net charge per residue, f+ − f− |
| `frac_GS` | fraction (0–1) | Combined glycine + serine fraction |

### 4.3 Repeat metrics (per LCR)

| Metric | Unit | Description |
|---|---|---|
| `n_RG_motifs`, `n_SR_motifs`, `n_GS_motifs` | count (integer) | Non-overlapping motif matches (RG-like = `RG{1,2}\|RPR`; RS-like = `RS\|SR`; GS-like = `GS\|SG`) |
| `*_repeat_cluster`, `*_repeat_region` | — (boolean) | ≥3 motifs with spacing ≤10 aa (Gerstberger et al. 2014) |

### 4.4 Distribution metrics (per LCR; per property track: polar, hydrophobic, aromatic, disorder, charged)

| Metric pattern | Unit | Description |
|---|---|---|
| `<track>_n_runs` | count (integer) | Number of contiguous property runs |
| `<track>_mean_run_length`, `<track>_max_run_length` | aa | Run-length statistics |
| `<track>_mean_gap` | aa | Mean gap between runs |
| `<track>_cv_gap` | unitless | Coefficient of variation of gap lengths |
| `<track>_label` | — (categorical) | `compact` (≤3 runs, max run > 35% of L), `dispersed`, `insufficient` |

### 4.5 Co-occurrence metric (per LCR)

| Metric | Unit | Description |
|---|---|---|
| `cooc_positive_aromatic` | fraction of windows (0–1) | Fraction of windows of size `max(5, L // 4)` containing ≥1 positive and ≥1 aromatic residue (cation–π potential) |

### 4.6 Derived analysis metrics (computed in Phase 2.4.1)

| Metric | Unit | Level | Description |
|---|---|---|---|
| `coverage_per_lcr` | fraction of protein length (0–1) | LCR | `length / uniprot_length` |
| `n_lcr` | count (integer); categorized 1 / 2 / 3+ | Protein | LCRs per protein per method |
| `lcr_residues` | aa | Protein | Σ LCR lengths per protein per method |
| `coverage` | fraction of protein length (0–1) | Protein | `lcr_residues / uniprot_length` |
| `property_fraction_<signature>` | fraction of protein length (0–1) | Protein | Residues in LCRs with primary signature *s* / `uniprot_length` |
| `log2_obs_exp` | log2 fold ratio (unitless) | Contingency-table cell | 0 = expected, +1 = 2× enriched, −1 = 2× depleted |
| `proportion_within_rna_class` | fraction (0–1) | Contingency-table cell | Column-normalized share within an RNA class |
| Summary statistics | source-metric unit | Group | n, median, Q1, Q3, IQR (= Q3 − Q1), mean, SD |
| `p_value` / `p_fdr` | probability (0–1) | Test | Chi-square / Kruskal–Wallis; BH-FDR within analysis block |

## 5. Key design decisions

- **Positional analysis**: five true position classes on Pfam-classified proteins; `unclassified_no_pfam` analyzed as a separate binary question.
- **Physicochemical analysis**: primary annotations only; both sequence-count and residue-count bases reported.
- `description` column excluded; measurements from the metrics workbook only.
- **No filtering**: AlcoR retained; no LCRFinder length threshold; overlaps fixed upstream.
- **Pseudo-replication**: every analysis at LCR level and protein level; protein-level significance is the robustness criterion.
- **Method as stratum**: all analyses per method; cross-method comparison via the consolidation layer.

## 6. Statistical contract

| Question type | Measure | Test |
|---|---|---|
| Categorical × RNA class | log2(observed/expected) per cell | Chi-square |
| Continuous × RNA class | Median + IQR and mean ± SD | Kruskal–Wallis |

BH-FDR per analysis block, α = 0.05; log2 uniform; effect sizes lead. Visualization: histograms and stacked bars only; fraction histograms grouped into multi-panel figures.

## 7. Scope limitations

- Analyses are **compositional over LCR-positive proteins** (no LCR-negative background; the master table lacks RNA-class labels).
- Categorical table p-values are **omnibus (table-level)**; cell claims rest on log2 effect sizes.
- Coverage assumes non-overlapping intervals within method (verified upstream).
- Class imbalance: mRNA dominates; `diverse`, `ribosomal protein`, and `unknown` are small — per-class n quoted on interpretation.

## 8. Data-integration validation

A staged join diagnostic (accession → +method → +start → +end → ±1 tolerance → interval overlap, per method) identified and corrected three integration defects:

1. **Method-label mismatch**: annotations `FLPS`/`SEG` vs. measurements `fLPS_strict`/`SEG_strict`; resolved by a canonical alias map on both sides.
2. **Accession extraction** from pipe-delimited composite IDs (text before the first `|`); initial parsing had silently broken the UniProt-length join.
3. **Multi-sheet measurements** (composition, distribution, co_occurrence): all three sheets now merged on interval keys.

Post-fix join QC: **100% measurement match and 100% UniProt-length match for all six methods**; exact interval equality 100% per method, confirming identical coordinate sets.

## 9. Outputs and organization

- Per-method results: `results/methods/<method>/{tables,figures}`.
- Cross-method consolidation: `results/summary/{tables,figures}` (master stacked tables, headline numbers, class × method pivots, signature-enrichment matrix, comparison heatmaps/bars).
- Deep-read extraction: `results/summary/deep_read/{tables,figures}` (charge-landscape summaries, significant coverage cells, concordance check).
- `RESULTS_GUIDE.md` and `README_deep_read.md`: reading order and question→output maps.

---

## 10. Findings

### 10.1 Detection methods fall into three archetypes (headline numbers)

| Method | n LCRs | n proteins | Median LCR length | Median coverage | % multi-LCR | Mean LCRs/protein |
|---|---|---|---|---|---|---|
| AlcoR | 2,804 | 968 | 25 aa | 0.295 | 61.1% | 2.90 |
| CAST | 1,880 | 1,055 | 76 aa | 0.294 | 44.6% | 1.78 |
| FLPS | 8,099 | 1,697 | 17 aa | 0.166 | 80.8% | 4.77 |
| LCRFinder | 6,371 | 1,008 | 6 aa | 0.047 | 86.2% | 6.32 |
| SEG_intermediate | 1,307 | 622 | 26 aa | 0.075 | 48.2% | 2.10 |
| SEG | 493 | 296 | 21 aa | 0.044 | 35.5% | 1.67 |

- **Long-window** (CAST), **short-motif** (LCRFinder), and **intermediate** (FLPS, AlcoR, SEG_intermediate, SEG) archetypes differ by an order of magnitude in call length and coverage.
- Cross-method comparisons are therefore meaningful only on **protein-level** quantities; LCR-level metrics compare different objects across methods.
- Method effects dominate RNA-class effects in scale; class structure is read *within* each method.

### 10.2 Validated cross-method results (protein level)

**Consensus board — final (6/6 methods assessed):**

| Signal | Support | Status |
|---|---|---|
| Ribosomal proteins: single basic/Arg-rich LCR at a domain boundary; no linker/terminal/no-Pfam LCRs | 5/5 powered methods (SEG underpowered, directional) | **Headline result** |
| pre-rRNA/snoRNA: basic + acidic LCRs | 6/6 | **Validated** |
| snRNA (and mRNA): SR/RS repeats; snRNA polyampholyte + linker tendency | 6/6 | **Validated** |
| mRNA: GS-rich, RG/RGG, aromatic-sticker, disorder-spacer composition; multi-LCR in short-call methods | 5/6 (composition); multi-LCR method-dependent | **Validated** |
| ncRNA: ~40% of proteins lack Pfam annotation | 6/6 (+0.94…+1.70) | **Validated — most consistent result** |
| ncRNA: aromatic-aggregation-prone LCRs | Strong cells in CAST and FLPS at both levels | **Validated** |
| tRNA: few, domain-intrinsic LCRs | 5/6 (CAST null = window-length artifact) | **Validated** |
| diverse: hydrophobic/aggregation-prone LCRs | 3/5 powered methods | **Validated, secondary thread** |
| tRNA: cation-π enrichment | contradicted at protein level | **Rejected (short-window artifact)** |
| Arg-specific (vs. generic basic) enrichment in nucleolar classes | 3/4 composition methods | Supported; claim charge-driven basic/Arg-rich |

**Headline result — the ribosomal-protein LCR architecture.** Across AlcoR, FLPS, CAST, LCRFinder, and SEG_intermediate, ribosomal proteins (and ribosome-biogenesis proteins more broadly) are domain-dense proteins carrying **exactly one dominant basic/Arg-rich LCR positioned at a domain boundary**: single-LCR enrichment in every method (up to 95.7% single in CAST), `domain_edge`/`domain_intrinsic` enriched (CAST edge +1.54, 69.6% of proteins; SEG_intermediate intrinsic +3.24), with near-total exclusion from interdomain linkers (0–0.8% of proteins), distal termini (2–12%), and no-Pfam proteins (≤2%).

**Nucleolar axis (pre-rRNA, snoRNA).** Both classes carry basic + acidic LCRs in all six methods. Acidic arm: `acidic_region` presence in 92.9% of SEG-detected pre-rRNA proteins (+1.89), 67.6% snoRNA in FLPS (+0.67). Basic arm: `basic_enriched_region` enriched in pre-rRNA and snoRNA in all methods (e.g., CAST +1.66/+1.71). The Arg-specific refinement holds in AlcoR/FLPS/CAST but not LCRFinder (6 aa peaks rarely satisfy the combined Arg rule) — the claim is charge-driven basic/Arg-rich.

**Splicing axis (snRNA, mRNA).** `sr_rs_repeat_region` enriched in snRNA in all six methods at protein level (+0.49…+1.26; 15–32% of proteins), with `mixed_charge_polyampholyte` support (FLPS +0.56, CAST +0.62) and an `interdomain_linker` tendency in snRNA (CAST +1.00, SEG_intermediate +0.77, FLPS +0.30), consistent with RS domains between RRMs in SR proteins.

**mRNA class.** Compositionally distinct: `gs_rich_neutral_region`, `rg_rgg_repeat_region`, `aromatic_sticker_region`, `disorder_rich_spacer_region` enriched across methods (disorder-spacer in 76% of LCRFinder mRNA proteins), basic signatures depleted. Multi-LCR architecture appears only in short-call methods (FLPS 72% 3+, LCRFinder 80% 3+) — method-dependent, reflecting call granularity.

**ncRNA class.** The single most consistent cross-method result: ~30–44% of ncRNA-class proteins lack Pfam annotation, enriched in all six methods (+0.94, +0.98, +1.11, +1.01, +1.70, +1.56). `aromatic_aggregation_prone_region` enrichment is machine-confirmed strong at both levels in CAST (+2.34/+2.35) and FLPS (+1.50/+1.69). ncRNA-binding proteins emerge as a partially uncharacterized, LCR-driven class — a priority for follow-up.

**tRNA class.** Few (mostly single) LCRs per protein (4/5 powered methods), `domain_intrinsic` placement (AlcoR +0.72, FLPS +0.57, LCRFinder +0.82; CAST null attributable to 76 aa windows exceeding domain size). tRNA enzymes are domain-dense with small embedded LCRs.

**diverse class.** Hydrophobic, aggregation-prone LCRs enriched in three powered methods at protein level (FLPS +1.10, LCRFinder +1.00, SEG_intermediate +3.53), with elevated `unmapped_physicochemical_properties` (SEG_intermediate +2.94) — atypical LCRs partially outside the current taxonomy; flagged for taxonomy review.

**unknown class.** Heterogeneous top signatures (Gln-rich, hydrophobic, Pro-rich), elevated multi-LCR and distal-terminal LCRs; retained as the control/QC narrative, as designed.

### 10.3 Continuous charge landscape validates the categorical results

The threshold-based signature calls are backed by distribution-level differences in charge metrics (per class × method, KW-tested):

- **Ribosomal basic charge is distribution-level, not a threshold artifact**: median NCPR positive in **all six methods** (+0.082…+0.267; the only class positive everywhere); median `frac_positive` the class maximum in 5/6 methods (CAST 0.33, FLPS 0.29, LCRFinder 0.37, SEG_intermediate 0.34, AlcoR 0.17).
- **The nucleolar classes are bipolar**: pre-rRNA combines high `frac_positive` (CAST 0.24, second to ribosomal) with the highest `frac_negative` of any class in 4–5/6 methods, and near-zero median NCPR in long-window methods (CAST +0.05, FLPS +0.06) with maximal FCR (class maximum in 5/6, 0.28–0.81). Near-zero NCPR + maximal FCR = coexisting basic and acidic tracts, not neutrality.
- **SEG variants resolve the acidic pole**: SEG pre-rRNA median NCPR = −0.74 (FCR 0.81); SEG_intermediate = −0.50. The earlier signature-level difference (SEG's top pre-rRNA call = `acidic_region` vs. CAST/FLPS = `basic_enriched_region`) is therefore a **stringency effect selecting which pole of the same bipolar landscape is detected**, not a method disagreement.
- **mRNA neutrality confirmed continuously**: median NCPR ≈ 0 in every method; FCR in the lower band — mRNA LCRs are genuinely the least charged population.
- **snRNA positive pole in SEG** (median NCPR +0.35): consistent with Arg within RS runs; SEG's strict windows isolate exactly those runs.
- **AlcoR is the only method with no class structure in FCR** (KW FDR 0.18) — expected, as it detects coiled-coil propensity rather than composition bias; its contribution is positional/structural.
- FLPS/LCRFinder show large within-class IQRs (bimodality from call granularity: many zero-charge micro-calls plus dense ones) — another reason protein-level summaries carry the interpretation.

### 10.4 Coverage and signature-occupancy statistics (KW-tested)

- **Coverage differences across RNA classes are significant in every method** (protein-level KW FDR: FLPS 6.4e-27, LCRFinder 9.9e-17, CAST 6.6e-06, SEG_intermediate 8.0e-05, SEG 7.2e-04, AlcoR 0.019). Ribosomal protein coverage is the class maximum in 4/6 methods (CAST tops at mRNA, FLPS at snoRNA); tRNA is the lowest/near-lowest in CAST, FLPS, and LCRFinder.
- **The ribosomal LCR occupies the largest share of its parent protein in all six methods at LCR level** (`coverage_per_lcr` class maximum: CAST median 26% of protein length, AlcoR 13%, SEG_intermediate 13%, FLPS 10%, SEG 8.7%, LCRFinder 2.4%).
- **Signature occupancy (property fractions among carrier proteins)**:
  - Arg-rich contact regions cover a median **28% of protein length in CAST and AlcoR** ribosomal proteins (vs. 2–5% elsewhere; AlcoR FDR 0.006, FLPS 0.129, FDR 1.3e-15); basic regions cover 13–28% in CAST/FLPS/SEG_intermediate ribosomal proteins. Combined with the single-LCR architecture and boundary placement: **one dominant basic LCR occupying a substantial fraction of the protein at a domain boundary** — an RNA-contact arm, not a small patch.
  - Acidic-region occupancy is the class maximum in pre-rRNA (FLPS 0.076, LCRFinder 0.026) and elevated in snoRNA.
  - SR/RS repeat regions occupy 12–14% of protein length in mRNA/snRNA carriers (SEG_intermediate, FDR 0.0075) — a substantial structural element.
  - Disorder-spacer *presence* is depleted in ribosomal proteins, but its *fraction among carriers* is the class maximum (CAST 0.12, SEG_intermediate 0.20): presence and magnitude are distinct statements; report both.
  - Rows with n = 1–4 proteins (e.g., FLPS gs_rich ribosomal n=3) are descriptive only and excluded from formal claims.

### 10.5 Machine-checked concordance validation

`concordance_check.csv` places LCR-level log2 (sequence-count and residue-count bases) beside protein-level log2 for every method × feature × class cell, with computed flags: `concordant` (same direction both levels), `discordant` (opposite non-zero directions), `strong_cell` (concordant **and** ≥2-fold at both levels, |log2| ≥ 1).

- **All headline cells are machine-confirmed**: ribosomal Arg-rich/basic (CAST strong: +2.67/+1.95 and +2.28/+1.59), domain-edge and domain-intrinsic placement (strong in CAST, SEG_intermediate, LCRFinder, SEG), distal-terminal and linker/no-Pfam exclusions (strong in 5/5 powered methods), nucleolar basic (CAST pre-rRNA/snoRNA strong), snRNA SR/RS (AlcoR strong +1.81/+1.45), ncRNA no-Pfam (strong in CAST and SEG_intermediate; concordant elsewhere).
- **Depletions are as consistent as enrichments**: basic/Arg-rich depleted in mRNA and snRNA; aromatic-sticker, GS-rich, and disorder-spacer depleted across the nucleolar classes; SR/RS and RG/RGG depleted in ncRNA, pre-rRNA, and tRNA. Classes systematically avoid each other's LCR types.
- **Discordant cells are predominantly near-zero sign flips** (sampling noise). One genuine pattern: FLPS ribosomal `basic_enriched_region` is LCR-level enriched (+1.05) but protein-level neutral (−0.06) — within-protein clustering of FLPS calls inflating LCR-level counts; a worked example of the pseudo-replication the two-level rule was designed to catch, cited in methods.
- **Calibration note**: the `strong_cell` bar is strict — several validated results are "concordant but not strong" because protein-level effects are genuinely smaller (e.g., AlcoR domain_edge × ribosomal +0.86/+0.47). Strong cells = headline claims; concordant cells = support layer.

### 10.6 Rejected hypotheses and method artifacts

- **tRNA cation-π enrichment**: large LCR-level log2 values (LCRFinder +3.44, SEG +4.36) trace to single short-window observations; protein-level tables flat in CAST and LCRFinder. **Rejected as a short-window artifact.**
- **CAST positional nulls**: 76 aa median windows cannot fit inside typical domains, suppressing `domain_intrinsic` calls and inflating boundary ambiguity (edge vs. intrinsic assignments depend on window length). Positional comparisons across methods must account for call length.
- **LCRFinder single-residue signatures** (glycine/serine-rich) are near-tautological at 6 aa median length; its motif-level signatures (SR/RS, RG/RGG) are its trustworthy contributions.
- **SEG extreme log2 values** derive from tiny windows and small counts; n quoted alongside throughout.
- **SEG strict's role** is the strictness control: patterns confirmed in SEG are high-confidence; its nulls are underpowered, not contradictory.
- **AlcoR** measures coiled-coil propensity rather than composition bias; its hydrophobic/domain-intrinsic signals read accordingly.

## 11. Outputs inventory (for documentation)

- `headline_numbers.csv` — per-method key figures (§10.1).
- `master_lcr_categorical.csv` / `master_protein_categorical.csv` — all enrichment tables, both levels.
- `master_lcr_continuous.csv` / `master_protein_continuous.csv` — all metric summaries with median/IQR/mean/SD/KW p/FDR.
- `median_coverage_by_class_method.csv`, `median_length_by_class_method.csv` — class × method pivots.
- `top_signature_per_class.csv` — strongest enriched signature per class per method (LCR level; table-level p).
- `charge_landscape_summary.csv` + `charge_landscape_median_<metric>.csv` — charge metrics per class × method with KW statistics (§10.3).
- `significant_coverage_cells.csv` — KW-tested coverage and property-fraction cells (§10.4).
- `concordance_check.csv` — machine-checked LCR↔protein concordance with strong-cell flags (§10.5).
- Figures: coverage/length/signature heatmaps (summary), NCPR/frac-positive median heatmaps (deep_read), per-method stacked bars, multi-panel fraction histograms.

## 12. Next steps

1. Detailed documentation/manuscript drafting from this record.
2. Optional per-cell Fisher-exact pass if cell-level significance claims are required for specific headline cells.
3. Follow-up analysis of the ncRNA (no-Pfam, aromatic-aggregation) and diverse (hydrophobic-LCR) classes — the two most novel threads.
4. Taxonomy review for hydrophobic/aggregation-prone LCRs given the `unmapped` enrichment in the diverse class.
