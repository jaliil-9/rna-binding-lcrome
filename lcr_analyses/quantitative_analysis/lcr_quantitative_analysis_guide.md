# Results Guide — (LCR diversity across RNA classes)

**Companion to:** `lcr_quantitative_analysis.py` outputs → `results/methods/<method>/tables/<method>_results.xlsx` + `figures/`.

## 1. What this analysis asks

How do LCR characteristics vary across **RNA-target superclasses** (mRNA, tRNA, pre-rRNA, snRNA, snoRNA, ncRNA, ribosomal protein, diverse, unknown), within the RBP set, per detection method?

**Reference design:** every contrast is **one-vs-rest** — each class is compared against all other classes pooled. There is no external (non-RBP) background in this phase.

## 2. Statistical contract 

| Question type | p-value | Effect size |
|---|---|---|
| Categorical (carrier rates, 2×2 per cell) | Fisher exact | log2 odds ratio (Haldane-corrected) |
| Continuous (fractions, coverage, charge) | Kruskal–Wallis omnibus per metric | Cliff's δ per class vs. rest |

- BH-FDR within analysis block × method; α = 0.05.
- `descriptive_only` = any cell/group n < 5 → read, never claim.
- **Protein level carries the claims; LCR level is descriptive** (pseudo-replication rule).
- Effect sizes lead; p-values gate.

## 3. The workbook, sheet by sheet

### `categorical_lcr` — LCR-level categorical enrichments (descriptive)

One row = one feature value × one RNA class. Covers: domain-position classes, Pfam-classified vs. not, physicochemical signatures (sequence-count and residue-count bases).

| Column | Meaning |
|---|---|
| `analysis` | Sub-analysis: position / no_pfam / signature, and the count basis |
| `feature_value` | The category tested (e.g. `domain_edge`, `sr_rs_repeat_region`) |
| `rna_class` | The class arm of the one-vs-rest contrast |
| `n_pos_group` / `n_tot_group` | LCRs with the feature / all LCRs, in this class |
| `n_pos_ref` / `n_tot_ref` | Same, all other classes pooled (the reference) |
| `rate_group` / `rate_ref` | The two frequencies — the raw interpretable numbers |
| `log2_or` | Effect size. 0 = parity; +1 = 2× odds vs. other classes; −1 = half |
| `p_value` | Fisher exact, this 2×2 cell |
| `p_fdr` | BH-adjusted within block × method |
| `descriptive_only` | n < 5 anywhere in the 2×2 → do not claim |

### `categorical_protein` — protein-level categorical enrichments (inferential)

Identical columns; unit = protein. Contains: LCR multiplicity (1 / 2 / 3+), position-presence, signature-presence (does the protein carry ≥1 LCR of that type). **This sheet is where categorical claims live.**

### `continuous_lcr` / `continuous_protein` — metric distributions per class

One row = one metric × one class. LCR level: length, coverage_per_lcr, all composition/property fractions, fcr/ncpr. Protein level: coverage, per-signature property fractions.

| Column | Meaning |
|---|---|
| `rna_class`, `metric` | Row identity |
| `n` | Observations in this class |
| `median` | Robust central tendency (fractions are heavily skewed) |
| `mean`, `sd` | Parametric summary, for literature comparability |
| `cliffs_delta_vs_rest` | Effect size vs. all other classes. Sign = direction; \|δ\| ≈ 0.11 / 0.28 / 0.43 ≈ small / medium / large |
| `p_value` | KW omnibus for the metric (same value repeated on each of the metric's class rows — it is a per-metric, not per-cell, test) |
| `p_fdr` | BH-adjusted within block |

### `qc` — audit trail

One-row overview (n LCRs, n proteins, measurement/length join rates) + proteins per RNA class. Every class-size caveat in the findings traces here.

## 4. Figures (unchanged)

`figures/` per method: stacked bars (position/signature proportions per class), multiplicity and coverage histograms, multi-panel fraction histograms, coverage-vs-n_LCR scatter. Histograms only, per plot policy.

## 5. Reading examples

- `categorical_protein`, row `signature_presence × sr_rs_repeat_region × snRNA`, log2_or = +1.0, p_fdr < 0.01 → snRNA proteins have 2× the odds of carrying an SR/RS-repeat LCR compared with all other RBP classes; significant after FDR; claimable.
- `continuous_lcr`, row `ncpr × ribosomal protein`, δ = +0.4 → RP LCRs are shifted positive in net charge vs. other classes, medium-large effect; the KW p on that metric tells you charge varies across classes overall.
- Any row with `descriptive_only = True` → quote n, make no claim.
