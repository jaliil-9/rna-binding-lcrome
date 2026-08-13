# Phase 2.4.1 — Results Guide

## Reading order
1. `summary/tables/headline_numbers.csv` — per-method key figures; check methods detect comparable amounts before interpreting biology.
2. `summary/figures/coverage_heatmap.png` and `length_heatmap.png` — do methods grossly agree on scale?
3. `summary/figures/signature_enrichment_heatmap.png` — cross-method view of physicochemical enrichment by RNA class.
4. Per-method folders (`methods/<method>/`) — drill-down archive; open only when a summary signal needs detail.

## Which output answers which question
| Analysis question | Summary layer | Per-method source |
|---|---|---|
| 1. Domain position x RNA class | `master_lcr_categorical.csv` (analysis=`position_lcr`), `master_protein_categorical.csv` (`position_presence`) | `lcr_position_stacked.png` |
| 2. Physicochemical signature x RNA class | `signature_enrichment_key.csv`, `top_signature_per_class.csv`, `signature_enrichment_heatmap.png` | `lcr_signature_*_stacked.png` |
| 3. Sequence-level properties | `master_lcr_continuous.csv`, `length_heatmap.png` | `lcr_length_histogram.png`, `lcr_fraction_histograms.png` |
| 4. LCR multiplicity per protein | `headline_numbers.csv` (`pct_multi_lcr`, `mean_n_lcr`), `n_lcr_by_method.png` | `protein_multiplicity_stacked.png`, `protein_n_lcr_histogram.png` |
| 5. Coverage and property fractions | `coverage_heatmap.png`, `master_protein_continuous.csv` | `protein_coverage_histogram.png`, `protein_signature_fraction_histograms.png` |

## Conventions
- Categorical enrichment: log2(observed/expected); 0 = expected, +1 = 2x enriched, -1 = 2x depleted.
- Continuous: median, IQR, mean, SD, Kruskal-Wallis p (BH-FDR within level and method).
- Robustness rule: a signal is robust if it holds at protein level; LCR level is supporting evidence.
