# Results Guide — (RBP LCRs vs. proteome background)

**Companion to:** `comparison.py` outputs → `results/methods/<method>/tables/<method>_results.xlsx` + `results/summary/tables/`.

## 1. What this analysis asks

Are the LCR features found within the RBP set **enriched relative to the full proteome** — globally (Mode A) and per RNA-target class (Mode B)? The reference is the **non-RBP background** (20,427-protein UniProt reviewed proteome, 1,669 RBP-labeled; see `phase2_5_background_construction.md`).

Two modes:
- **Mode A (group):** all RBPs vs. background — what marks the RBPome overall.
- **Mode B (class):** each RNA class vs. background — what marks each class.

## 2. Statistical contract

| Question type | p-value | Effect size |
|---|---|---|
| Binary (carrier rates) | Fisher exact | log2 odds ratio (Haldane-corrected) |
| Continuous | Mann–Whitney | Cliff's δ (group vs. background) |

MW is KW with two groups — the same test. BH-FDR per block × method; α = 0.05; effect sizes lead; n < 5 → `descriptive_only`. Rung 1 (proteome-wide) is descriptive only: the proteome contains the RBPs, so it cannot be tested against.

## 3. The workbook, sheet by sheet

### `A_incidence` — does the RBPome carry LCRs at all? (binary)

| Column | Meaning |
|---|---|
| `block` | `A1_incidence` (overall) or `A1_incidence_by_length` |
| `group` | The arm tested (`True` = RBPs) |
| `length_bin` | Background-length quintile; blank = unstratified. The confounder control — RBPs are longer on average |
| `n_pos_group` / `n_tot_group` | Carrier / total proteins in the arm |
| `n_pos_ref` / `n_tot_ref` | Same in the background |
| `rate_group` / `rate_ref` | Carrier rates — the headline percentages |
| `log2_or`, `p_value`, `p_fdr`, `descriptive_only` | Fisher effect size, p, FDR, n-floor |

### `A_architecture` — how much LCR do carriers carry? (continuous, among carriers)

Rows: `coverage`, `n_lcr`, `mean_lcr_length`. Columns per arm: `n`, `median`, `mean`, `sd`; then `cliffs_delta`, `p_value` (MW), `p_fdr`, `descriptive_only`.

### `A_composition` — what are RBP LCRs made of? (continuous, LCR level)

One row per composition metric (20 aa fractions + 10 property metrics + length), RBP LCRs vs. background LCRs. Same continuous columns. Contains the charge landscape (rows `ncpr`, `fcr`, `frac_positive`, `frac_negative`).

### `A_signatures` — which LCR types are RBP-specific? (mixed)

Per signature, two row sets distinguished by `measure`:
- `carrier_rate` (binary columns): is the signature enriched among RBP proteins?
- `occupancy_among_carriers` (continuous columns): among carriers, what fraction of the protein does the signature occupy? Presence and magnitude are separate statements.

### `B_incidence` / `B_composition` / `B_signatures` — per RNA class

Identical structures; `group` = RNA-target superclass, reference = non-RBP background. B1: which classes are LCR-driven at all (enabled by full-protein class labels). B2: the external charge/composition landscape per class. B3: class signature specificity vs. background.

## 4. Summary level (`results/summary/tables/`)

### `headline_numbers.csv`
Per method: n LCRs, carrier counts and totals per arm — the denominator record for every rate.

### `validation_grid.csv` — Phase 2.4.1 × Phase 2.5

One row per signature × class × method. Both sides carry the **same metric (log2 OR)** with different references:

| Column | Meaning |
|---|---|
| `internal_log2_or` / `internal_fdr` | Phase 2.4.1: class vs. other RBP classes |
| `external_log2_or` / `external_fdr` | Phase 2.5: class vs. background |
| `n_carriers` / `n_class` | External-side power |
| `verdict` | See below |

**Verdicts:** `externally_validated` (enriched in both) · `rbp_generic` (internal enrichment that is an RBPome-wide feature, not class-specific) · `internal_only_artifact` (internal enrichment contradicted externally) · `new_external_signal` (internally null/depleted but above background) · `validated_depletion` (depleted in both) · `concordant_null` · `external_underpowered` (n < 5 — absence of evidence, not rejection).

### `validation_consensus.csv`
Verdict tallies per class × signature across the six methods — the cross-method board.

## 5. Reading examples

- `A_signatures`, row `rg_rgg_repeat_region × carrier_rate`, log2_or = +3.2 → an RBP has ~9× the odds of carrying an RG/RGG-repeat LCR vs. a background protein; the signature is RBP-marking.
- `B_incidence`, row `tRNA`, unstratified −0.91 and length bins similarly negative → tRNA-binders are genuinely LCR-avoidant, not just short.
- Validation grid: internal +0.5 / external +4.2 → class specificity stronger against the proteome than against other RBPs — internal obs/exp was conservative.

## 6. Known scope limits

- FLPS incidence is saturated (~87% of all proteins); AlcoR measures coiled-coil propensity, not composition — both are near-null instruments for incidence.
- Domain-position and RBD stratification are out of scope for rungs 1–2 (return at rung 3).
- Composition/architecture are computed among carriers; non-carriers enter only through incidence.
