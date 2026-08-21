# RBP LCRs against a proteome-wide background: Context, Design, and Findings

**Status:** analysis complete (rungs 1–2); digest finalized.
**Date:** August 20, 2026.

---

## 1. Objective and design

Quantitative analysis was compositional over LCR-positive RBPs (no LCR-negative background; findings.md §7). Background analysis re-measures the same LCR universe against a **proteome-wide background**, in three planned rungs:

1. Proteome-wide (descriptive only — the proteome contains the RBPs, so no inferential test is possible)
2. **RBP vs. non-RBP background** (primary inferential contrast; all reported tests)
3. Consensus non-RBP census (master ∩ RBP2GO ∩ RBPbase) — *not yet performed*

Two analysis modes share one engine:

- **Mode A (group):** all RBPs vs. background — incidence, architecture, composition/charge, signature spectrum.
- **Mode B (class):** each RNA-target superclass vs. the non-RBP background, same blocks, plus per-class incidence (enabled by full-protein RNA-class labels, `rbp_rna_classification.xlsx`).
- **Validation grid:** quantitative analysis internal log2(obs/exp) joined against the external log2 odds ratio per class × signature × method.

## 2. Background construction (summary; full record in ensembl_background_construction.md)

- One protein per gene; Ensembl r116 canonical route evaluated and **rejected** (36.5% of RBP genes would be represented by a non-UniProt-canonical isoform; 8.4% absent).
- **Adopted:** UniProt reviewed reference proteome (UP000005640) = 20,416 proteins + 11 appended master sequences → **20,427 proteins**, 1,669 RBP-labeled (1,463 direct + 199 gene-level patch + 11 appended).
- All six method configurations (CAST, SEG, SEG_intermediate, FLPS, LCRFinder, AlcoR) run on this FASTA with RBP-run parameters: **184,966 LCRs**; `pc_analysis_annotation.py` run once over the full set → single signature taxonomy, single cation-π quantile.
- Join QC: 100% length match; 1,669/1,669 RBPs class-labelled; 249 classified proteins unmapped (documented unresolved set).

## 3. Statistical contract

Fisher exact + Haldane-corrected log2 OR (binary); Mann–Whitney + Cliff's delta (continuous); BH-FDR per block × method; α = 0.05; effect sizes lead; n < 5 → `descriptive_only`. Length confounding addressed by reference-length-quintile stratification of incidence. Scripts: `phase2_5_compare.py`, `phase2_5_digest.py`; outputs under `results/methods/<method>/tables/mode{A,B}_*.csv`, `results/summary/tables/validation_grid.csv`, `results/digest/`.

## 4. Findings

### 4.1 Headline incidence (Mode A)

| Method | RBP carriers | Background carriers |
|---|---|---|
| CAST | 55.7% (930/1669) | 43.4% |
| FLPS | 89.3% | 87.5% (saturated) |
| LCRFinder | 53.1% | 43.7% |
| SEG | 15.8% | 8.5% |
| SEG_intermediate | 32.6% | 22.4% |
| AlcoR | 43.7% | 43.4% (null) |

**LCRs are modestly RBP-enriched, not RBP-exclusive.** Every method with dynamic range shows higher RBP carrier rates; FLPS saturates (89% of all proteins), AlcoR is flat (coiled-coil propensity, not composition bias). The RBP character is in *which* LCRs, not *whether* LCRs.

### 4.2 Signature specificity of the RBPome (Mode A; median log2 OR across methods)

| Signature | Support | Median log2 OR | Interpretation |
|---|---|---|---|
| rg_rgg_repeat_region | 6/6 enriched | +3.19 | ~10× enriched; RBP-marking (3.2% vs 0.3% carriers) |
| sr_rs_repeat_region | 5/6 | +2.78 | RBP-marking (4.3% vs 0.3%) |
| basic_enriched_region | 6/6 | +1.42 | Charge regions mark the RBPome globally |
| arginine_rich_rna_contact_region | 6/6 | +0.99 | " |
| acidic_region | 6/6 | +0.91 | Enriched but background-common (13.1% vs 8.1%) — not specific |
| mixed_charge_polyampholyte | 5/5 | +0.77 | Enriched |
| polar_linker | 3/5 | +0.78 | Modest, method-mixed |
| gs_rich_neutral_region | 5/6 | +0.46 | Modest |
| disorder_rich_spacer_region | 3/6 | +0.33 | Modest |
| aromatic_sticker_region | 2/6 (+1 depl.) | +0.47 | Heterogeneous |
| cation_pi_rich_neighborhood | 0/4 (all positive) | +0.53 | Direction only; quantile-based signature weak at carrier level |
| Single-AA (S, G, Q) | 0 enriched | ≈ 0 | Proteome-generic |
| **hydrophobic_region** | **4/4 depleted** | **−1.37** | **RBPs avoid hydrophobic LCRs** (3.9% vs 6.3%) |
| unmapped_physicochemical_properties | null | −0.06 | Taxonomy coverage comparable across arms |

### 4.3 Class incidence vs background (B1; per-class carrier rates)

| Class (n) | Enriched/powered | Median log2 OR | Verdict |
|---|---|---|---|
| mRNA (658) | 5/5 | +1.33 | **Enriched — most LCR-rich class; length-independent** (sig. in length bins 2–5 in 4 methods) |
| snoRNA (39) | 4/5 | +1.59 | Enriched; strongest fold, small n |
| snRNA (93) | 4/5 | +0.90 | Enriched |
| pre-rRNA (113) | 3/6 | +0.84 | Enriched (powered methods) |
| ncRNA (118) | 1/6 | +0.61 | Weakly enriched |
| unknown (284) | 3/6 | +0.43 | Mildly enriched — consistent with uncharacterized/LCR-driven |
| ribosomal protein (157) | 0/5; 4/5 depleted | −1.26 | **Depleted vs pooled background, but largely length-confounded**: within the shortest bin the depletion attenuates and loses significance in most methods (CAST bin ≤211: 26.7% vs 20.0%, ns; FLPS +0.35 ns; LCRFinder −0.96, FDR 0.08). Honest formulation: low carrier rate, high occupancy among carriers (cf. §4.6) |
| **tRNA (162)** | **0/5; 3/5 depleted** | **−0.91** | **True, length-independent depletion** — persists in length bins across methods (CAST −1.69**/−1.18**/−1.87** in bins 3–5; LCRFinder −1.20*/−1.43**/−2.10**; FLPS −2.03** bin 3; SEG_int bins 3–4 sig.). tRNA-binding enzymes exclude LCRs at matched protein lengths |
| diverse (45) | 0/4 | −0.15 | Null |

### 4.4 External validation board (Phase 2.4.1 internal log2 × external log2 OR)

Externally validated cells (class × signature, support across 6 methods):

- **6/6:** mRNA × rg_rgg_repeat_region (external median +4.16); unknown × acidic_region (+1.18)
- **5/6:** mRNA × sr_rs (+3.88), polar_linker (+1.59), gs_rich (+1.19), disorder_rich_spacer (+0.99); snRNA × sr_rs (+4.14); pre-rRNA × acidic (+1.99); snoRNA × acidic (+2.52)
- **4/6:** mRNA × aromatic_sticker (+1.25); pre-rRNA × basic (+2.31); snRNA × basic (+2.00); unknown × basic (+1.56)

**Systematic observation:** external log2 ORs are larger than internal log2(obs/exp) for the same cells — the within-RBP marginal baseline of Phase 2.4.1 was conservative; class specificity is stronger against the proteome.

**New external signals (internally invisible):** mRNA × acidic_region (+1.09, 4/6 new) and mRNA × basic_enriched_region (+1.08, 4/6 new). These were internally *depleted* (below the RBP average) yet lie *above the proteome background*. **Within-RBP depletion ≠ absence of enrichment vs. background** — both statements required in the manuscript framing.

**RBP-generic reclassifications:** cation-π neighborhoods, hydrophobic regions (diverse), disorder-spacers (ncRNA), and parts of the basic/Arg-rich signal are RBPome-wide features rather than class-specific ones. The map now has two layers: what marks the RBPome globally (§4.2) vs. what distinguishes classes (§4.4).

**Ribosomal-protein compositional claims:** no RP cell reaches ≥4-method external validation; basic/Arg-rich carrier-rate enrichment is largely RBP-generic or underpowered. The RP claim now rests on (i) validated depletions (disorder_rich_spacer 3/6; acidic directional), (ii) the Phase 2.4.1 positional/single-LCR architecture (out of Phase 2.5 scope by design), and (iii) occupancy among carriers (§4.6, pending). Reword, not retract.

### 4.5 Class charge landscape vs background (B2)

*Pending readout of `digest_charge_landscape.csv` — external test of the §10.3 charge claims (ribosomal NCPR-positive, nucleolar bipolar, mRNA neutral).*

### 4.6 Occupancy among carriers (key signatures)

*Pending readout of `digest_occupancy_key_signatures.csv` — magnitude vs. presence for Arg-rich (ribosomal), acidic (nucleolar), aromatic-aggregation (ncRNA).*

## 5. Scope limitations (Phase 2.5-specific)

- `external_underpowered` = absence of evidence (n < 5 or no FDR), not rejection; rare signatures × small classes land there by construction (ncRNA aromatic-aggregation: 1 validated + 4 underpowered — not externally confirmed *at carrier-rate level*).
- The validation grid tests carrier rate; occupancy/magnitude claims live in §4.6.
- FLPS incidence is saturated and AlcoR composition-free; both are near-null instruments for incidence.
- 249 class-labelled proteins unmapped to the background (documented); 199 gene_patch + 11 appended proteins flagged for a sensitivity rerun (not yet executed).
- Domain-position and RBD stratification deliberately out of scope (see thread decision, 2026-08-18); the ribosomal positional headline is therefore not externally testable here.

## 6. Next steps

1. Read `digest_charge_landscape.csv` and `digest_occupancy_key_signatures.csv`; complete §4.5–4.6.
2. Sensitivity rerun excluding gene_patch/master_appended proteins (`--sensitivity` flag, to implement).
3. ncRNA aromatic-aggregation and tRNA-depletion threads: decide occupancy-level follow-up.
4. Rung 3: consensus non-RBP census vs. RBP2GO 2.0 + RBPbase; RBD stratification becomes available there.
