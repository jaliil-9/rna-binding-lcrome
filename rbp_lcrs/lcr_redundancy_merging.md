# Within-method merging of LCR calls

**Phase:** 2.4 — Map LCRs' 3-layer annotation
**Date:** 2026-08-12
**Script:** `merge_lcr_within_method.py`
**Input:** `lcr_methods_combined.xlsx` (sheet `all_results`, raw calls — retained as the frozen reference)
**Output:** `lcr_methods_merged.xlsx` (`all_results` = merged table, `summary` = stats, one sheet per method)

## Rationale

Each LCR caller produces redundant calls within a single protein through three
mechanisms: boundary jitter (same region, ±1–2 aa), containment (one call inside
another), and exact duplicates. These inflate per-method statistics (counts,
composition distributions), so calls are merged **within each
(protein_id × method) group only** — never across methods, never across
proteins. Cross-method structure is preserved deliberately: methods are analyzed
separately, and a "detected by k/6 methods" consensus axis is planned as a
robustness layer.

## Rules

| Tier | Condition | Action |
|------|-----------|--------|
| T0 | Identical (start, end) | Keep one copy |
| T1 | Call contained inside another | Keep the outer interval |
| T2 | Overlap ≥ 0.7 of the **shorter** interval | Merge into union span |

- Overlap metric is overlap-of-shorter, not IoU, so a short call is never
  absorbed by a large neighbor on a small shared patch.
- Sweep is single-pass, intervals sorted by (start, end).
- Tier 3 (fragment bridging across small gaps) is **off** in this version.

## Provenance columns (on every merged row)

- `n_members` — number of raw calls absorbed (1 = untouched)
- `member_intervals` — coordinate spans, e.g. `5-24;5-25`
- `merge_tier` — rule(s) applied, e.g. `T1T2`, `none` for singletons
- `raw_row_ids` — 0-based row indices into the input `all_results` sheet

## Results

| method | merged_calls | multi_member | raw_calls_absorbed | pct_redundant |
|---|---|---|---|---|
| AlcoR | 2804 | 0 | 2804 | 0.0 |
| CAST | 1880 | 415 | 2509 | 25.1 |
| LCRFinder | 6371 | 0 | 6371 | 0.0 |
| SEG_intermediate | 1307 | 0 | 1307 | 0.0 |
| SEG_strict | 493 | 2 | 495 | 0.4 |
| fLPS_strict | 8099 | 2009 | 11490 | 29.5 |

**Total: 20,954 merged calls from 24,670 raw calls (~15% redundant overall).**

## Interpretation

Merge rates match each caller's known mechanism:

- **CAST (25.1%)** and **fLPS_strict (29.5%)** — iterative window scanning and
  boundary refinement produce shifted/contained re-calls; merging removes these
  detector artifacts.
- **SEG_strict (0.4%)** and **SEG_intermediate (0.0%)** — SEG emits
  non-overlapping segments by construction; near-zero redundancy is expected and
  confirms the merge is not over-aggressive.
- **AlcoR (0.0%)** and **LCRFinder (0.0%)** — no within-method overlap, but
  their redundancy is *fragmentation* (one tract split into short calls across
  small gaps), which is Tier 3 territory and intentionally not applied yet.

## Caveats / next steps

1. `sequence` is blank in merged rows — union spans must be re-sliced from
   parent protein sequences.
2. All composition/annotation columns (`frac_*`, `fcr`, `ncpr`, labels,
   physicochemical annotations) must be **recomputed** on merged spans —
   never averaged from members.
3. Domain-context columns (nearest Pfam, distance, flanking) likewise need
   recomputation against the new coordinates.
4. Evaluate Tier 3 (gap-bridging, same-class only) separately for
   LCRFinder/AlcoR; calibrate gap threshold from the observed gap distribution.
5. Merging precedes all downstream overlays (pLDDT, phyloP, PTMs) — those
   aggregate per-residue data over these coordinates.
