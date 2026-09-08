# Physicochemical Behavior Taxonomy for RBP Low-Complexity Regions

**Scope.** This document specifies the rule-based physicochemical annotation system for low-complexity regions (LCRs) identified in human RNA-binding proteins (RBPs). Each LCR receives one or more behavior signatures derived from sequence-encoded physicochemical properties. Signatures are **hypothesis-oriented descriptors of expected physical behavior**, not direct functional assignments. They are designed to be combined downstream with (i) RNA-target superclass and (ii) domain-position class, yielding a three-layer functional description: *RNA target → domain location → physicochemical behavior*.

**Implementation.** `pc_analysis_annotation.py`, with residue properties from `aa-physicochemical-properties.csv`.

---

## 1. Measuring logic: two levels

The system operates on two levels. Residue-level assignments feed sequence-level metrics; signatures are rules over sequence-level metrics.

### 1.1 Residue level — property encoding

Each of the 20 amino acids is encoded with five overlapping (non-exclusive) physicochemical attributes, based on Taylor/Jalview-style property sets, with charge defined at pH 7.4:

| Attribute | Values | Membership |
|---|---|---|
| **Charge (pH 7.4)** | Pos / Neg / Neu | Pos: R, K · Neg: D, E · Neu: all others (H treated as neutral) |
| **Polar** | Yes / No | R, N, D, Q, E, H, K, S, T, W, Y |
| **Hydrophobic** | Yes / No | A, C, G, H, I, K, L, M, F, T, W, V, Y |
| **Aromatic** | Yes / No | F, Y, W, H |
| **Disorder-promoting** | Yes / No | A, R, G, Q, S, P, E, K |

The disorder-promoting set follows the TOP-IDP disorder/order classification (Campen et al. 2008, *Protein Pept Lett* 15:956, [doi:10.2174/092986608785849164](https://doi.org/10.2174/092986608785849164)): disorder-promoting residues A, R, G, Q, S, P, E, K; order-promoting W, C, F, I, Y, V, L, N; neutral D, H, M, T. The use of such binary residue groupings for sequence analysis is established practice:

> "One can generalize the calculation of patterning parameters to any arbitrary binary sequence-patterning parameter... Examples include hydrophobic patterning, whereby all residues are grouped into hydrophobic or non-hydrophobic sets, disorder- or order-promoting residues, or neutral polar residues or all other residues."
> — Holehouse et al. (2017), *Biophys J* 112:16–21, [doi:10.1016/j.bpj.2016.11.3200](https://doi.org/10.1016/j.bpj.2016.11.3200)

### 1.2 Sequence level — three metric families

Given a residue-encoded LCR sequence of length *L*:

**Composition** (what is present)
- Per-residue fractions `frac_X` for all 20 amino acids
- Property fractions: `f_polar`, `f_hydrophobic`, `f_aromatic`, `f_disorder`
- `frac_strong_hydro`: combined fraction of strong hydrophobes (V, I, L, M, F, W, Y)
- `frac_GS`: combined Gly + Ser fraction
- Charge metrics: `f+` (fraction positive), `f−` (fraction negative),
  **FCR** = f+ + f−, **NCPR** = f+ − f− (Das & Pappu 2013; Holehouse et al. 2017)

**Distribution** (how it is arranged)
For each property track (polar, hydrophobic, aromatic, disorder, charged), run/gap statistics are computed: number of runs, mean/max run length, mean gap, gap CV. Each track is labeled:
- **compact** — few runs with at least one long block (≤3 runs and max run > 35% of *L*)
- **dispersed** — otherwise
- **insufficient** — property absent

**Co-occurrence** (what co-localizes locally)
- Sliding-window score for the **positive × aromatic** pair (cation–π potential): fraction of windows of size `max(5, L // 4)` containing at least one positively charged and at least one aromatic residue. Scores are classified relative to the dataset-wide 75th percentile.

**Repeat detection**
- Non-overlapping motif matches (regex), followed by a cluster criterion: a repeat region requires **≥3 motifs with consecutive start-to-start spacing ≤10 residues**, following the Gerstberger et al. (2014) definition of RG/RGG regions as *"at least three RG/RGG repeats spaced 10 amino acids or fewer apart"* (*Nat Rev Genet* 15:829–845, [doi:10.1038/nrg3813](https://doi.org/10.1038/nrg3813)).
- Motif sets: **RG-like** = `RG{1,2}|RPR`; **RS-like** = `RS|SR`; **GS-like** = `GS|SG`.

---

## 2. Global thresholds

All thresholds are provisional, dataset-calibrated values pending a dedicated calibration phase; literature-derived values are marked †.

| Parameter | Value | Used by |
|---|---|---|
| `f_charge_min` | 0.25 | basic / acidic |
| `ncpr_min` | 0.10 | basic / acidic |
| `frac_R_min` | 0.12 | Arg-rich basic |
| `fcr_min` | 0.30 * | polyampholyte (Das & Pappu 2013) |
| `ncpr_neutral_max` | 0.05 | polyampholyte |
| `f_polar_min` | 0.50 | polar linker |
| `fcr_polar_linker_max` | 0.30 | polar linker |
| `f_hydro_low` | 0.20 | polar linker |
| `frac_strong_hydro_min` | 0.30 | hydrophobic region |
| `f_arom_min` | 0.08 | aromatic sticker / cluster |
| `f_disorder_min` | 0.60 | disorder-rich |
| `frac_S_min` / `frac_Q_min` / `frac_G_min` / `frac_P_min` | 0.30 / 0.20 / 0.25 / 0.15 | single-AA enrichment |
| `f_GS_min` | 0.40 | GS-rich |
| `frac_GS_each_min` | 0.10 | GS-rich (both G and S required) |
| `min_repeats` / `max_repeat_spacing` | 3 * / 10 * | repeat clusters (Gerstberger et al. 2014) |
| co-occurrence cutoff | dataset 75th percentile | cation–π |

Each LCR may receive multiple signatures; the first matched rule (in the order below) is primary, the rest secondary. Triggered rules are reported with their metric values in the `annotation_evidence` field.

---

## 3. Signature catalogue

### 3.1 Repeat-driven signatures

#### 3.1.1 `rg_rgg_repeat_region`
**Rule.** RG-like repeat cluster: ≥3 `RG/RGG/RPR` motifs spaced ≤10 aa apart.
**Metrics.** `n_RG_motifs`, `RG_repeat_cluster`; thresholds `min_repeats = 3`, `max_repeat_spacing = 10` (Gerstberger et al. 2014).
**Interpretation.** Multivalent RNA-binding and phase-separation driver; preferential recognition of G-quadruplex-containing RNA.

> "RGG/RG motifs are RNA binding segments found in many proteins that can partition into membraneless organelles. They occur in the context of low-complexity disordered regions and often in multiple copies."
> — Chong, Vernon & Forman-Kay (2018), *J Mol Biol* 430:4650–4665, [doi:10.1016/j.jmb.2018.06.014](https://doi.org/10.1016/j.jmb.2018.06.014)

> "Although short RGG/RG-containing regions can sometimes form high-affinity interactions with RNA structures, multiple RGG/RG repeats are generally required for high-affinity binding, suggestive of the dynamic, multivalent interactions that are thought to underlie phase separation in formation of cellular membraneless organelles."
> — Chong, Vernon & Forman-Kay (2018), ibid.

> "The arginine/glycine-rich RGG domain, found in over 1000 human RBPs, has been identified as a key driver of phase separation, and has been shown to be implicated in the recognition and preferential binding towards mRNAs containing G-quadruplex (G4) structure."
> — Vidal Ceballos et al. (2025), *Sci Rep* 15:9295, [doi:10.1038/s41598-025-88499-y](https://doi.org/10.1038/s41598-025-88499-y)

#### 3.1.2 `sr_rs_repeat_region`
**Rule.** RS-like repeat cluster: ≥3 `RS/SR` motifs spaced ≤10 aa apart.
**Metrics.** `n_SR_motifs`, `SR_repeat_cluster`; same cluster thresholds.
**Interpretation.** Phosphorylation-gated RNA/protein interaction switch; hallmark of SR-family splicing factors.

> "SR proteins are characterized by the presence of a C-terminal domain enriched with the Arginine (R) and Serine (S) amino acid sequences (RS domain) and an N-terminal RNA recognition domain (RRM domain)."
> — Jeong (2017), *Mol Cells* 40:1–9, [doi:10.14348/molcells.2017.2319](https://doi.org/10.14348/molcells.2017.2319)

> "Phosphorylation switches the RS domain of the serine/arginine-rich splicing factor 1 from a fully disordered state to a partially rigidified arch-like structure."
> — Xiang et al. (2013), *Structure* 21:2162–2174, [doi:10.1016/j.str.2013.09.014](https://doi.org/10.1016/j.str.2013.09.014)

#### 3.1.3 `gs_rich_neutral_region`
**Rule.** GS-like repeat cluster (≥3 `GS/SG` motifs, spacing ≤10), **or** `frac_GS ≥ 0.40` with both `frac_G ≥ 0.10` and `frac_S ≥ 0.10`.
**Metrics.** `frac_GS`, `frac_G`, `frac_S`, `GS_repeat_cluster`; thresholds `f_GS_min = 0.40`, `frac_GS_each_min = 0.10`.
**Interpretation.** Chemically neutral spacer behaving approximately as an ideal (Gaussian) chain; flexible connector between interaction modules.

> "Recent work has implicated the polar residues glycine (G) and serine (S) as chemically neutral spacer residues, in agreement with observations that (GS) repeat sequences behave as physical instantiations of ideal (Gaussian) chains."
> — Holehouse et al. (2021), *Biochemistry* 60:3566–3581, [doi:10.1021/acs.biochem.1c00465](https://doi.org/10.1021/acs.biochem.1c00465)

### 3.2 Charge-driven signatures

Charge logic follows the Das–Pappu framework, in which FCR and NCPR discriminate conformational classes:

> "Polyampholytes are either strong (FCR ≥ 0.3) or weak (FCR < 0.3) and can be neutral (NCPR ∼ 0) or have a net charge."
> — Das & Pappu (2013), *PNAS* 110:13392–13397, [doi:10.1073/pnas.1304749110](https://doi.org/10.1073/pnas.1304749110)

> "Polyelectrolytic IDPs with NCPR below a threshold value of 0.25 adopt compact globular ensembles, whereas sequences that lie above this threshold adopt well-solvated expanded coils and even stiff rod-like conformations."
> — Holehouse et al. (2017), *Biophys J* 112:16–21, [doi:10.1016/j.bpj.2016.11.3200](https://doi.org/10.1016/j.bpj.2016.11.3200)

#### 3.2.1 `basic_enriched_region`
**Rule.** `f+ ≥ 0.25` and `NCPR ≥ 0.10`.
**Metrics.** `frac_positive`, `ncpr`; thresholds `f_charge_min = 0.25`, `ncpr_min = 0.10`.
**Interpretation.** Net-positive polyelectrolyte adopting expanded, solvated conformations; electrostatic interaction potential with the negatively charged RNA backbone.

#### 3.2.2 `arginine_rich_rna_contact_region`
**Rule.** Basic-region conditions **and** `frac_R ≥ 0.12`.
**Metrics.** as above plus `frac_R`; threshold `frac_R_min = 0.12`.
**Interpretation.** Specific RNA-contact module: arginine mediates base-specific hydrogen bonding and π-stacking in addition to backbone electrostatics, distinguishing it from lysine-driven basic regions.

> "Arginine can interact with nucleotide bases via hydrogen bonding and π-stacking; thus, nucleotide conformers that provide access to the bases provide enhanced opportunities for RGG interactions."
> — Chong, Vernon & Forman-Kay (2018), *J Mol Biol* 430:4650–4665, [doi:10.1016/j.jmb.2018.06.014](https://doi.org/10.1016/j.jmb.2018.06.014)

#### 3.2.3 `acidic_region`
**Rule.** `f− ≥ 0.25` and `NCPR ≤ −0.10`.
**Metrics.** `frac_negative`, `ncpr`; thresholds `f_charge_min = 0.25`, `ncpr_min = 0.10`.
**Interpretation.** Net-negative functional module; local clusters of acidic charge are documented drivers of complex coacervation and regulatory interactions.

> "Pak et al. used localCIDER to identify local clusters of high negative charge in the disordered region of the Nephrin intracellular domain that drives phase separation via complex coacervation."
> — Holehouse et al. (2017), *Biophys J* 112:16–21, [doi:10.1016/j.bpj.2016.11.3200](https://doi.org/10.1016/j.bpj.2016.11.3200)

#### 3.2.4 `mixed_charge_polyampholyte`
**Rule.** `FCR ≥ 0.30` and `|NCPR| ≤ 0.05` and charge distribution **dispersed**.
**Metrics.** `fcr`, `ncpr`, charged-track distribution; thresholds `fcr_min = 0.30` (Das & Pappu 2013), `ncpr_neutral_max = 0.05`.
**Interpretation.** Strong polyampholyte with well-mixed opposite charges: salt-sensitive, conformationally tunable chain; the dispersed-pattern requirement acts as a low-κ proxy.

> "The fraction of charged residues discriminates between weak and strong polyampholytes. Using atomistic simulations, we show that weak polyampholytes form globules, whereas the conformational preferences of strong polyampholytes are determined by a combination of fraction of charged residues values and the linear sequence distributions of oppositely charged residues."
> — Das & Pappu (2013), *PNAS* 110:13392–13397, [doi:10.1073/pnas.1304749110](https://doi.org/10.1073/pnas.1304749110)

### 3.3 Polarity / hydropathy signatures

#### 3.3.1 `polar_linker`
**Rule.** `f_polar ≥ 0.50` and `f_hydrophobic ≤ 0.20` and `FCR ≤ 0.30` and polar distribution **dispersed**.
**Metrics.** `f_polar`, `f_hydrophobic`, `fcr`, polar-track distribution; thresholds `f_polar_min = 0.50`, `f_hydro_low = 0.20`, `fcr_polar_linker_max = 0.30`.
**Interpretation.** Uncharged, solvent-exposed polar spacer with low intrinsic cohesion; the FCR cap restricts the signature to the uncharged polar spacers of the stickers-and-spacers framework.

> "Charged and polar residues interspersed between stickers do not display strong interaction patterns, and instead act as spacers that mediate the contacts among stickers."
> — Martin et al. (2020), *Science* 367:694–699, [doi:10.1126/science.aaw8653](https://doi.org/10.1126/science.aaw8653)

#### 3.3.2 `hydrophobic_region`
**Rule.** `frac_strong_hydro ≥ 0.30` (V, I, L, M, F, W, Y) and hydrophobic distribution **compact**.
**Metrics.** `frac_strong_hydro`, hydrophobic-track distribution; threshold `frac_strong_hydro_min = 0.30`.
**Interpretation.** Cohesive, aggregation-prone segment; atypical within LCRs, which are generally depleted of hydrophobic content. Hits are interpreted as candidate local core/aggregation modules.

> "IDPs fail to fold autonomously, their sequences are deficient in hydrophobic groups and enriched in polar and charged residues."
> — Das & Pappu (2013), *PNAS* 110:13392–13397, [doi:10.1073/pnas.1304749110](https://doi.org/10.1073/pnas.1304749110)

### 3.4 Aromatic / interaction-sticker signatures

#### 3.4.1 `aromatic_sticker_region`
**Rule.** `f_aromatic ≥ 0.08` and aromatic distribution **dispersed**.
**Metrics.** `f_aromatic`, aromatic-track distribution; threshold `f_arom_min = 0.08`.
**Interpretation.** Uniformly patterned aromatic stickers supporting multivalent π-mediated interactions; promotes liquid–liquid phase separation while inhibiting aggregation.

> "We also show that uniform patterning of aromatic residues is a sequence feature that promotes LLPS while inhibiting aggregation."
> — Martin et al. (2020), *Science* 367:694–699, [doi:10.1126/science.aaw8653](https://doi.org/10.1126/science.aaw8653)

> "Evolutionary analysis has argued that evenly distributed hydrophobic and/or aromatic residues facilitate liquid-like condensates and prevent aggregation."
> — Holehouse et al. (2021), *Biochemistry* 60:3566–3581, [doi:10.1021/acs.biochem.1c00465](https://doi.org/10.1021/acs.biochem.1c00465)

#### 3.4.2 `aromatic_aggregation_prone_region`
**Rule.** `f_aromatic ≥ 0.08` and aromatic distribution **compact**.
**Metrics.** as above.
**Interpretation.** Clustered aromatic patterning; favors amorphous aggregation over liquid-like condensates.

> "While Aro^Perfect formed spherical droplets..., Aro^Patchy formed large amorphous aggregates... the uniform distribution of aromatic residues along LCDs favors solubility and LLPS over aggregation."
> — Martin et al. (2020), *Science* 367:694–699, [doi:10.1126/science.aaw8653](https://doi.org/10.1126/science.aaw8653)

#### 3.4.3 `cation_pi_rich_neighborhood`
**Rule.** Positive × aromatic co-occurrence score in the **top quartile** of the dataset (75th-percentile cutoff).
**Metrics.** `cooc_positive_aromatic` (window `max(5, L // 4)`); cutoff computed per dataset.
**Interpretation.** Frequent local positive–aromatic contacts (effectively Arg–Tyr/Phe in LCRs); stabilizes structures and enables binding to ligands, drugs, and nucleic acids; a core interaction mode in RBP condensates.

> "Previous studies identified arginine and tyrosine as stickers in FUS and other FET family proteins. An analytical model shows that the saturation concentrations of these proteins are inversely proportional to the product of the numbers of arginine and tyrosine residues in a given sequence."
> — Martin et al. (2020), *Science* 367:694–699, [doi:10.1126/science.aaw8653](https://doi.org/10.1126/science.aaw8653)

### 3.5 Disorder signature

#### 3.5.1 `disorder_rich_spacer_region`
**Rule.** `f_disorder ≥ 0.60` over the disorder-promoting set (A, R, G, Q, S, P, E, K).
**Metrics.** `f_disorder`; threshold `f_disorder_min = 0.60`.
**Interpretation.** Entropy-rich, conformationally heterogeneous spacer; the residue set follows the TOP-IDP disorder-propensity classification (Campen et al. 2008, [doi:10.2174/092986608785849164](https://doi.org/10.2174/092986608785849164)), applied as a binary residue grouping per localCIDER practice (Holehouse et al. 2017, quote in §1.1).

### 3.6 Single-residue enrichment signatures

**Rules.** `serine_rich_region`: `frac_S ≥ 0.30` · `glutamine_rich_region`: `frac_Q ≥ 0.20` · `glycine_rich_region`: `frac_G ≥ 0.25` · `proline_rich_disordered_region`: `frac_P ≥ 0.15`.
**Metrics.** per-residue fractions; thresholds as listed (provisional, dataset-calibrated).
**Interpretation.** Composition-dominated regions characteristic of LCRs: Ser-rich (phosphorylation substrate potential), Gln-rich (prion-like character), Gly-rich (flexibility, RG-context), Pro-rich (backbone rigidification within disorder, SH3/WW-motif potential). These are empirical composition descriptors; their thresholds are calibrated on the present dataset rather than taken from literature constants.

### 3.7 Fallback

#### 3.7.1 `unmapped_physicochemical_properties`
**Rule.** Assigned when no signature fires.
**Interpretation.** No detectable compositional, distributional, co-occurrence, or repeat bias at current thresholds; reported as unmapped rather than forced into a behavior class.

---

## 5. Threshold status

Thresholds marked * derive directly from the cited literature (`fcr_min`, repeat-cluster parameters). All others are provisional values calibrated qualitatively on the current RBP LCR dataset and will be revisited in a dedicated calibration phase. The annotation schema (signatures, metric families, two-level logic) is independent of specific threshold choices.
