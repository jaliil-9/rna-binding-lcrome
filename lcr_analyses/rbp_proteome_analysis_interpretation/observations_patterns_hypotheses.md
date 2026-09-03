# Observations, Patterns, Hypotheses

**Scope.** Consolidated discussion record of the RNA-binder LCR pilot: external RBP-vs-background comparison (Mode A), per-class external analysis (Mode B), the three-layer annotation system (RNA target → domain position → physicochemical signature), unsupervised clustering, feature-contribution analysis, and cross-method agreement.

**Data basis.** `project_outcome_report.md` (six detection methods: CAST, fLPS, LCRFinder, SEG, SEG_intermediate, AlcoR; internal + external tables; PAM global k=4 and carriers k=5); `lcr_annotations.xlsx` (20,954 RBP-arm LCRs with primary/secondary signatures, positions, sequences); `protein_clusters.csv` × 5 (carrier-level PAM k=2–10, hierarchical, HDBSCAN labels for CAST, fLPS, LCRFinder, SEG, SEG_intermediate). Background: 20,427 UniProt reviewed proteins, 1,669 RBP-labeled (8.2%).

**Statistical conventions.** Binary: Fisher exact + Haldane-corrected log2OR. Continuous: Mann–Whitney/Kruskal–Wallis + Cliff's δ. BH-FDR per block × method. `descriptive_only` (n < 5) rows are not inferential. RNA class was never a clustering feature.

---

## 1. Mode A — what distinguishes RBP-associated LCRs

### Observations
- **Incidence.** RBPs carry LCRs more often than non-RBPs in every composition-sensitive method (SEG +1.01, SEG-int +0.75, CAST +0.72, LCRFinder +0.55 log2OR; all FDR < 1e-12), robust to length stratification. fLPS is saturated (~88% everywhere); AlcoR is null (it measures coiled coils, not composition).
- **Architecture (carriers).** RBP carriers have more LCRs, higher coverage, longer LCRs (δ +0.04 to +0.19 depending on method; AlcoR flat).
- **Composition.** The most concordant result in the project: RBP LCRs are more polar (δ +0.11–0.28), more disorder-prone, more charge-rich (FCR δ +0.15–0.24), modestly cationic (NCPR δ +0.05–0.15), hydrophobe-depleted — all six methods, same direction.
- **Signature carrier rates (A3).** Universal enrichment: RG/RGG (6/6, log2OR +2.9 to +4.7), SR/RS (5/6), basic, Arg-rich, acidic (6/6), GS-rich and polyampholyte (5/6). Universal depletion: hydrophobic (6/6; CAST effectively excludes it: 2 RBP carriers vs 290 background). Single-AA enrichment signatures (Gln/Gly/Pro/Ser) are flat — RBP LCR identity is patterning/charge, not raw composition excess.
- **Occupancy ≪ carriage.** Occupancy among carriers is weak (δ ≈ 0.1–0.17, mostly null) while carrier odds ratios reach 2–20×. The RBPome's signature landscape is a prevalence phenomenon.

### Patterns
- The RBP–non-RBP distinction is not "has LCR" (43% of background proteins carry CAST LCRs) but: higher carrier probability + more elaborate architecture + a shifted physicochemical regime.
- Effect sizes scale with caller stringency (strict SEG largest ORs); base rates scale with leniency (fLPS). Cross-method OR comparison is qualitative only.
- Cation–π signature is fragile (dataset-relative 75th-percentile cutoff; 4-carrier enrichments in LCRFinder/SEG are descriptive-only). Aromatic sticker is the one genuine cross-method contradiction (fLPS −0.35* depletion vs LCRFinder +0.70* enrichment) — unresolved, likely a length-context effect on the 0.08 aromatic-fraction rule.

### Hypotheses
- **H1.** LCR propensity is permissive; signature content is discriminative.
- **H2.** The polar/cationic/disordered triad reflects electrostatic RNA compatibility + solubility constraints in RNA-rich environments; hydrophobe depletion is the negative-design signature.
- **H3.** RNA interaction is recruited at the protein level (carrier rate), not dosed within proteins (occupancy). Corollary: class-level claims should rest on carrier-rate logic; occupancy is descriptive.
  - *Why occupancy stays flat despite more/longer LCRs:* occupancy is per-signature and length-normalized; extra LCR mass spreads across 18 signature buckets (one primary per LCR), RBP protein length inflates the denominator, and conditioning on carriers deletes the extensive margin where the action is.
- **H4.** The A2 composition shift may be partly mediated by the enriched-signature subset (acidic/basic/Arg LCRs are charge-rich by construction) — test by recomputing A2 within signature-homogeneous strata.
- **Exception to H3:** ribosomal carriers — Arg occupancy δ +0.67, basic δ +0.65 (3–4/5 methods significant; medians 2.4–4.6× background carriers). The ribosomal basic/Arg island is dosage-loaded, not just recruited.

---

## 2. Mode B — per-class external profile

### B1 incidence (class vs background)
LCR-enriched classes: snoRNA (+1.4 to +2.1 in 4/6), pre-rRNA (3/6), mRNA (5/6), snRNA (4/6), ncRNA (2/6, weakest), unknown (3/6, modest). LCR-depleted: ribosomal proteins (6/6 significant, −0.75 to −2.30) and tRNA (4/6). diverse: null (n=45, underpowered). Depletions are absolute (below the 43% background baseline), not relative-RBPome effects.

### B2 composition — charge balance is the class axis
All carrying classes share the polar/disorder/charge-rich envelope; they separate on charge balance:
- Extreme cationic polyelectrolyte: ribosomal (NCPR δ +0.44, driven by positive fraction +0.46).
- Charge-maximal, net-neutral (strong-polyampholyte profile): pre-rRNA (FCR +0.57, both charge fractions up) and snoRNA (FCR +0.48), with aromatic depletion.
- Moderately cationic: mRNA (NCPR +0.11), snRNA (+0.16).
- tRNA's rare LCRs are charge-rich but balanced (FCR +0.32, NCPR ~0).
- diverse: flat (no significant composition effect).

### B3 signatures and internal×external verdicts
- **mRNA:** broad spectrum — RG/RGG +4.2 (6/6), SR/RS +3.9 (5/6) plus polar linker, aromatic sticker, GS, acidic, basic, polyampholyte, disorder-spacer; only hydrophobic depleted. Verdicts: SR/RS externally validated 5/6, RG/RGG 4/6, GS 3/6, aromatic sticker 2/6, disorder-spacer 2/6 — the strongest class-defining record.
- **snoRNA:** RG/RGG +4.8 (5/6; SEG +7.29 = largest class×signature OR in the project), basic +3.0 (3/6), acidic +2.5 (5/6). Verdicts: acidic/basic EV 2/6 (tentative); RG/RGG external-only.
- **snRNA:** SR/RS +4.1 (5/6), RG/RGG +4.0 (5/6), basic +2.0 (4/6). No ≥2-method validated verdicts — above background but not distinctive within the RBPome.
- **pre-rRNA:** basic +2.3 (4/6), acidic +2.0 (5/6), Arg +1.5 (4/6). Verdicts: acidic EV 4/6, basic EV 3/6, Arg EV 2/6 — the most internally-anchored charge class.
- **unknown:** diluted generic palette (acidic 6/6, basic 5/6, Arg 5/6, modest magnitudes). Nothing depleted.
- **ncRNA:** acidic/basic enriched externally but every verdict is new_external_signal (internally null) → RBP-generic charge, no class distinctiveness. Internal RG/RGG depletion (CAST −2.47, FDR 0.044). fLPS aromatic-aggregation-prone EV is anecdote-level (≈3 carriers).
- **ribosomal:** small basic/Arg island (2–3/6) amid broad depletion (aromatic-sticker −2.4, GS −2.4, cation–π −2.2, disorder-spacer −1.8, acidic −1.6, hydrophobic −1.9). Only Arg-rich EV 2/6.
- **tRNA:** nothing up; validated disorder-spacer depletion (2/6).
- **diverse:** flat; underpowered, not evidence of absence.

### Class-size caveat (verdict-layer bias)
The internal reference is ~1/3 mRNA proteins, so repeat-rich small classes (snRNA n=93) show near-zero internal point estimates for RG/RGG (≈0 vs external +3.3 to +5.5) — reference composition, not just power. For SR/RS, snRNA's internal estimates are positive (+0.7 to +1.8) but underpowered. Fix: read verdicts with point estimates; run a direct snRNA-vs-mRNA contrast or a pooled repeat-axis (mRNA ∪ snRNA) vs rest test. The verdict system structurally couples class size to interpretability (diverse is nearly all external_underpowered).

---

## 3. The three modules

### Module 1 — repeat axis (mRNA, snRNA)
RG/RGG and SR/RS repeats define a shared spliceosomal/mRNA-processing interaction axis; mRNA is its statistical center of mass. snRNA carries the same repeats at the same rates and forms no separate cluster anywhere → functional subset of the axis (H-R2 holds). Where repeats are dosed, it is modest: mRNA SR/RS occupancy ~2× background carriers (δ +0.30, 2/6). snRNA's SR/RS tilt relative to RG/RGG is the open, power-limited question.
- **H-R1.** RG/RGG (multivalent G4-RNA contact) + SR/RS (phospho-gated switches) = one interaction axis for mRNA/spliceosome biology.
- **H-R2.** snRNA proteins are a subset of that axis (shared machinery; not a separate regime).

### Module 2 — charge patch (pre-rRNA, snoRNA)
Acidic + basic co-enrichment with the dataset's most extreme charge composition (FCR-maximal, NCPR≈0, polar-maximal, aromatic-depleted) = strong-polyampholyte profile. pre-rRNA is the internally anchored class (EV 4/6 acidic, 3/6 basic).
- **H-C1 (simplified).** Charge-dense, net-balanced LCRs are reversible droplet glue for assembling the ribosome-building machinery (nucleolus physics; NPM1 paradigm).
- **Cis/trans resolution (from `lcr_annotations.xlsx`).** True cis-polyampholyte LCRs are rare (3–7% of charge-family LCRs). The duality is carried by charge-dual proteins: ~50% of pre-rRNA/snoRNA multi-LCR carriers hold separate acidic + basic/Arg LCRs in one chain (NPM1-style), with the rest splitting into acidic-only and basic-only proteins. In mRNA/snRNA (and snoRNA, LCRFinder) charge-duality is significantly enriched above independence even after controlling for LCR count; in pre-rRNA it sits at independence (near-saturated marginals suffice). Strict SEG-int sees pre-rRNA as acidic-only (−5.19 anti-pairing) — the basic side is caller-dependent.
- **H-C2.** snoRNA's RG/RGG recruitment (SEG +7.29) points to G4 recognition in pre-rRNA/snoRNA guide biology — an extension of the repeat axis into the nucleolus.

### Module 3 — LCR avoidance (ribosomal, tRNA)
- Incidence depleted in all informative methods (ribosomal 6/6; tRNA 4/6); down to ~2% carriage under strict SEG.
- Length stratification: tRNA avoidance is length-robust (depleted within mid/long bins across CAST/fLPS/LCRFinder/SEG-int). Ribosomal depletion is unambiguous unstratified but noisy within bins (small, short-skewed class) — magnitude to be re-stated after a length-matched subset comparison.
- The ribosomal island: residual carriers keep basic/Arg patches (caller-split rates; occupancy 2.4–4.6× — the H3 exception), compositionally the most cationic LCRs in the dataset, positioned domain-edge/domain-intrinsic.
- **H-A1.** Ribosomal LCR biology is structural (Arg/Lys tails contacting rRNA cores), explaining the inverted presence-vs-dosage rule.
- **H-A2.** tRNA binders are purely RBD-centric; validated disorder-spacer depletion suggests active exclusion of disordered spacers around folded tRNA-binding domains.

### Non-module classes
- **unknown:** heterogeneous mixture ≈ diluted RBPome average; internally enriched for distal-terminal, 3+ multiplicity, cation–π. Expect clustering to fragment it.
- **ncRNA:** mild incidence enrichment; generic charge signatures (all new_external_signal); no cluster of its own; internal RG/RGG depletion. Not a fifth regime.
- **diverse:** underpowered (n=45).

---

## 4. Clustering

### Setup reminder
Per-method, per-analysis (global = all RBPs; carriers = LCR carriers only); PAM + hierarchical + HDBSCAN; k=2–6 initially; user-side manual inspection of all k and all three algorithms led to PAM as primary, k={4,5}, cross-algorithm ARI = 1. The report consolidates PAM global k=4 and carriers k=5. RNA class is metadata only.

### Observations
- Global k=4 in every method isolates a zero-signature non-carrier mega-cluster; ribosomal/tRNA classes collapse into it (71–98%). Global clustering ≈ Mode A incidence; the biology resolves at carriers k=5.
- Carriers k=5 recovers a recurring archetype vocabulary across the five compositional callers (AlcoR's clusters are uniformly aromatic-sticker and uninterpretable in this taxonomy):
  1. acidic, net-negative (all 5 methods)
  2. basic/Arg, net-positive (all 5; holds 52–70% of ribosomal carriers)
  3. neutral disorder-spacer (all 5; SEG splits it in two)
  4. sticker/hydrophobic (CAST, fLPS, LCRFinder)
  5. charge-dual, NCPR≈0 (fLPS, LCRFinder only)
  6. repeat clusters RG/RGG and SR/RS (SEG, SEG-int only; lenient callers dissolve repeats into spacer/charge clusters)
  7. cation–π (CAST only)
- Module confrontation: avoidance confirmed (non-carrier collapse + ribosomal basic island); charge patch confirmed with the cis/trans structure intact (pre-rRNA → acidic/basic clusters, snoRNA → charge-dual/acidic); repeat axis partially confirmed (strict-caller resolution only; snRNA merges with mRNA there — H-R2 supported).
- Class→cluster mapping is sharp for small extreme classes (pre-rRNA 93% into one SEG cluster; ribosomal 70%; tRNA 67%), diffuse for mRNA/unknown.

### Feature contribution (rebuilt feature panel vs reported k=5 labels)
Feature panel reconstructed exactly from `lcr_annotations.xlsx` (sequence-derived composition metrics; validated: per-cluster FCR/NCPR medians match the report to 2 decimals). Three measures per method (Kruskal ε², mutual information, ExtraTrees permutation importance; classifier recovers labels at 88–93% balanced accuracy):
- **Univariate:** FCR (ε² 0.64) and NCPR (ε² 0.52) are the top separators → the charge axis is the primary driver. 
- **Conditional:** continuous charge metrics are ~redundant (permutation importance ≈ 0.01) given the charge-signature binaries; the practical drivers are `sig_acidic` (0.116), `sig_disorder_spacer` (0.077), positions (distal-terminal 0.069, domain-adjacent 0.055), `sig_basic` (0.042), `sig_aromatic_sticker` (0.039), `sig_Arg` (0.039).
- **Block ordering:** signatures ≫ positions > composition > architecture ≈ 0 — retroactively validates the length-exclusion design.
- Parsimony option: drop FCR/NCPR/frac± (absorbed by binaries) or vice versa; the one thing continuous charge uniquely adds is the charge-dual gradient (NCPR≈0, FCR-high).
- Caveat: rebuilt PAM (documented Gower recipe) reproduces reported boundaries at ARI 0.32–0.49 (script-level weighting not in the workflow doc), but all cluster archetypes and medians match — conclusions concern the reported partition.

### Hypotheses
- **H-cl1.** ~6 recurrent archetypes are the candidate functional classes for the pilot's classification goal.
- **H-cl2.** Charge sign (NCPR) is the primary organizing axis of carrier space; charge density (FCR) and repeat/spacer content are secondary.
- **H-cl3.** tRNA-carrier LCRs concentrate in sticker/hydrophobic clusters (37% in both LCRFinder and fLPS) — a genuinely different LCR usage; inspect those proteins directly.

---

## 5. Cross-method agreement

- Raw pairwise ARI on shared carriers: 0.08–0.17 across detector families (SEG vs SEG-int 0.37 — same lineage). Exact partitions do not transfer; each caller defines a different LCR population.
- Archetype-level consensus across the 1,213 proteins covered by ≥2 callers: mean agreement 0.61; full 5/5 consensus for 19.5% of proteins.
- Stability tiers: acidic (pairwise recurrence 0.55; 33% full consensus) and basic/Arg (0.40; 34%) are the robust core; neutral-spacer is solid (0.49) but blurs with sticker/hydrophobic (mutual confusion 0.44); charge-dual is resolved only by lenient callers (12% recurrence, 1% full); repeat clusters exist only in strict callers (0% full consensus); cation–π is a CAST-only artifact.
- Confusion directions are chemically sensible (repeat-SR → basic/Arg or charge-dual where no repeat cluster exists; charge-dual → acidic where no dual cluster exists).
- Class tie-back: ribosomal → basic/Arg (most consistent class: 0.81 mean, 57% full); pre-rRNA → acidic 42% + basic 22% + dual 21%; snoRNA → acidic 43% + dual 25%; tRNA/diverse → sticker/hydrophobic; mRNA → spacer plurality (42%) but spread; snRNA/ncRNA/unknown spread without dominant type.

### Consequence for the pilot deliverable
The functional classification should be consensus archetypes with per-protein agreement as a confidence score; the ~20% full-consensus proteins form the high-confidence core. Report: acidic and basic/Arg as stable classes; sticker/hydrophobic with the spacer-boundary caveat; repeats as strict-caller-resolved classes; charge-dual with the cis/trans analysis attached.

---

## 6. Consolidated hypothesis board

| ID | Statement | Status |
|---|---|---|
| H1 | LCR propensity permissive, signature content discriminative | Supported (Mode A) |
| H2 | Polar/cationic/disordered triad = RNA compatibility + solubility design | Supported (Mode A, all methods) |
| H3 | Recruitment at protein level (carriage) ≫ dosage (occupancy) | Supported except ribosomal island |
| H4 | Composition shift mediated by enriched-signature subset | Open (stratify A2 by signature) |
| H-R1 | RG/RGG + SR/RS = shared mRNA/spliceosome interaction axis | Supported |
| H-R2 | snRNA is a subset of the repeat axis, not a separate regime | Supported (RG/RGG); SR/RS tilt open (power) |
| H-C1 | Charge-patch LCRs = droplet glue for ribosome-biogenesis machinery | Refined: mostly protein-level duality (trans-in-cis), cis-LCRs rare |
| H-C2 | snoRNA RG/RGG = G4 recognition in the nucleolus | Open (external-only signal) |
| H-A1 | Ribosomal island = structural rRNA-contact tails (dosage-loaded) | Supported; length-matched restatement pending |
| H-A2 | tRNA = purely RBD-centric, active spacer exclusion | Supported |
| H-cl1 | ~6 consensus archetypes = candidate functional classes | Supported, with stability tiers |
| H-cl2 | NCPR = primary axis of carrier space | Supported (top univariate; conditionally redundant with binaries) |
| H-cl3 | tRNA-carrier LCRs = sticker/hydrophobic mode | Open; inspect members |

## 7. Caveats and limitations
- FLPS saturates incidence; AlcoR measures coiled coils (null instrument for composition); strict SEG has tiny carrier sets (n=296 proteins). Cross-method effect sizes are qualitative.
- Cation–π and aromatic-sticker signatures need threshold/length-context review before claims lean on them.
- Ribosomal depletion magnitude needs a length-matched subset; verdict grid structurally disfavors small classes.
- Occupancy/architecture are carrier-conditional quantities (selection differs between arms).
- Archetype mapping thresholds were hand-set (documented and reproducible); agreement is computed on shared-carrier proteins (overweights strong-LCR proteins).
- Clustering boundaries are method-dependent (ARI 0.08–0.37); only archetypes transfer.

## 8. Next steps
1. Direct snRNA-vs-mRNA SR/RS contrast (or pooled repeat-axis test) to settle the H-R2 residue.
2. Length-matched subset restatement of ribosomal/tRNA avoidance.
3. Stratified A2 (composition within signature-homogeneous LCR sets) for H4.
4. Consensus-archetype assignment table with per-protein confidence (plurality fraction), incl. high-confidence core export.
5. Member inspection: tRNA sticker/hydrophobic carriers (H-cl3); snoRNA RG/RGG carriers vs known G4-binders (H-C2); medoid recurrence (e.g., Q96DP5 basic/Arg medoid in both LCRFinder and fLPS).
6. Threshold review for cation–π (dataset-relative cutoff) and aromatic-sticker (length context) before publication claims.
7. Optional: lncRNA-aware ncRNA subclassification to test whether the ncRNA heterogeneity hides an LCR mode in a subclass.
