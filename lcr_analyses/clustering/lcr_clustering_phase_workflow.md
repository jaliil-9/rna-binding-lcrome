# LCR clustering phase: workflow and design decisions

## Purpose

This phase performs unsupervised clustering of RNA-binding proteins (RBPs) from their low-complexity-region (LCR) phenotypes.

The goal is exploratory: identify whether proteins form interpretable groups from LCR architecture, domain position, physicochemical signature, and—in LCR carriers—LCR composition. RNA-target superclass labels are retained as metadata only. They are **not** input features and are not used to choose clustering parameters in this phase.

The phase is designed to be simple, method-specific, reproducible, and directly consistent with the preceding internal RBPome and external proteome analyses.

---

## Scope and decisions

### Per-method analysis only

Each LCR detection method is processed independently:

- AlcoR
- CAST
- FLPS
- LCRFinder
- SEG
- SEG_intermediate

There is no cross-method feature consensus, pooled segmentation, or merged clustering dataset at this stage.

The comparison between methods is deferred until all methods have independently produced their cluster assignments and summaries.

### Two complementary analyses

Each method produces two separate clustering analyses.

| Analysis | Protein population | Biological question |
|---|---:|---|
| **Analysis A: global LCR phenotype** | All RBPs, including proteins without a detected LCR | Do broad LCR architectures divide the RBPome? |
| **Analysis B: carrier phenotype** | Only proteins with at least one LCR from the current method | Do LCR-carrying proteins separate into detailed positional, signature, and composition phenotypes? |

The analyses are intentionally separate. Non-carriers have meaningful values for LCR count and coverage, but do not have a biologically defined mean LCR length or LCR composition. Filling those quantities with zero would create an artificial compositional phenotype.

### Protein is the clustering unit

One row in every clustering matrix represents one protein.

Raw annotation and measurement tables contain one row per detected LCR. These LCR-level rows are aggregated before clustering so a protein with multiple LCRs does not contribute multiple observations. This follows the existing analysis principle that protein-level features support the main biological interpretation, while LCR-level measurements are descriptive context.

### No new biological feature system

The clustering uses the architecture, position, signature, and physicochemical metrics already generated for the internal and external analyses. The only additional computational operation is aggregation of existing LCR-level physicochemical values to one protein-level value in Analysis B.

### Simple length handling

Raw protein length is not included as a clustering feature.

Length is handled through existing feature choices:

- LCR coverage is normalized by protein length.
- Analysis A uses a coarse LCR multiplicity category instead of raw LCR count.
- Protein length is retained in output files for later quality control and interpretation.

No regression residualization, length-derived composite metric, or length-matched clustering is used in this initial workflow.

---

## Input files

### Core inputs

| Input | Role |
|---|---|
| `lcr_annotations.xlsx` | LCR coordinates, LCR length, detection method, RNA class, domain-position class, and primary physicochemical signature |
| `lcr_features.xlsx` | Per-LCR physicochemical composition metrics |
| `combined_rbp_pfam38_rbpdb_modomics_uniprot.xlsx` | Protein universe and UniProt protein lengths |
| `rbp_rna_classification.xlsx` | RNA-target superclass metadata for output and later interpretation |

### Required annotation fields

The clustering scripts standardize the following annotation fields:

```text
proteinid
sourcemethod
start
end
length
rnatargetsuperclass
domainpositionclass
primaryphysicochemicalannotation
```

### Required LCR property metrics

Analysis B uses the existing physicochemical property panel:

```text
frac_polar
frac_hydrophobic
frac_strong_hydro
frac_aromatic
frac_disorder
frac_positive
frac_negative
fcr
ncpr
frac_GS
```

---

## Analysis A: global LCR phenotype

### Population

All RBPs in the master RBP universe are included for each detection method.

A protein without a detected LCR for the current method remains in the matrix. It receives:

```text
n_lcr = 0
lcr_residues = 0
coverage = 0
multiplicity_cat = "0"
all domain-position presence values = 0
all signature-presence values = 0
```

### Feature panel

Analysis A uses only protein-level LCR architecture and categorical-presence features.

| Feature block | Variables | Representation |
|---|---|---|
| Multiplicity | `multiplicity_cat` | Categorical: `0`, `1`, `2`, `3+` |
| LCR amount | `coverage` | Continuous |
| Domain context | domain-intrinsic, domain-edge, domain-adjacent, interdomain-linker, distal-terminal | Binary presence per protein |
| Signature identity | full primary-signature catalogue | Binary presence per protein |

### Excluded features

The following are intentionally excluded from Analysis A:

- `mean_lcr_length`, because it is undefined for non-carriers;
- LCR composition metrics, because they are undefined for non-carriers;
- raw protein length;
- raw amino-acid fractions;
- custom LCR-complexity, mutation-to-repeat, residualized, or PCA-derived features.

### Interpretation target

Analysis A is expected to reveal broad LCR architectures, for example:

- no-LCR / low-LCR proteins;
- single-LCR proteins;
- multi-LCR proteins;
- proteins with domain-associated LCRs;
- proteins with terminal or interdomain LCRs;
- proteins carrying repeat-, acidic-, basic-, or arginine-rich LCR signatures.

These are phenotype descriptions, not predefined RNA-target classes.

---

## Analysis B: LCR-carrier phenotype

### Population

Only proteins with `n_lcr > 0` for the current detection method are included.

This makes mean LCR length and LCR composition defined for every included protein.

### Feature panel

| Feature block | Variables | Representation |
|---|---|---|
| Architecture | `n_lcr`, `coverage`, `mean_lcr_length` | Continuous |
| Domain context | domain-intrinsic, domain-edge, domain-adjacent, interdomain-linker, distal-terminal | Binary presence per protein |
| Signature identity | full primary-signature catalogue | Binary presence per protein |
| LCR composition | 10 existing property metrics | Continuous, LCR-residue-weighted protein average |

### Protein-level LCR composition

For a protein with multiple LCRs, each LCR-level physicochemical metric is aggregated with LCR length as its weight:

\[
\mathrm{protein\ property} =
\frac{\sum_{\ell=1}^{n_{\mathrm{LCR}}}
\left(\mathrm{LCR\ length}_{\ell} \times \mathrm{property}_{\ell}\right)}
{\sum_{\ell=1}^{n_{\mathrm{LCR}}}\mathrm{LCR\ length}_{\ell}}.
\]

This preserves the existing composition measurements while giving one value per protein. It prevents a protein containing several short LCRs from being treated as several independent clustering observations.

### Signature catalogue

The initial scripts use the following existing primary-signature values:

```text
rg_rgg_repeat_region
sr_rs_repeat_region
gs_rich_neutral_region
basic_enriched_region
arginine_rich_rna_contact_region
acidic_region
mixed_charge_polyampholyte
polar_linker
hydrophobic_region
aromatic_sticker_region
aromatic_aggregation_prone_region
cation_pi_rich_neighborhood
disorder_rich_spacer_region
serine_rich_region
glutamine_rich_region
glycine_rich_region
proline_rich_disordered_region
unmapped_physicochemical_properties
```

A signature that is absent for a given method is kept as an all-zero column during feature construction and should then be removed before clustering if it has no variance.

---

## Mixed-data distance

PAM and hierarchical clustering use the same simple Gower-style mixed distance.

### Feature treatment

| Feature type | Treatment |
|---|---|
| Binary position/signature presence | Hamming mismatch distance |
| Categorical multiplicity category in Analysis A | Exact-match / mismatch distance |
| Continuous features | Robust scaling followed by Manhattan distance |

Continuous features are robust-scaled with `RobustScaler`, which centres by the median and scales using the interquartile range. This reduces the influence of skewed distributions without adding a new biological transformation.

### Feature-block weighting

The current implementation gives the broad feature types comparable influence rather than allowing large signature or composition blocks to dominate solely because they contain more columns.

For Analysis A, the blocks are:

- binary position/signature presence;
- categorical multiplicity;
- continuous coverage.

For Analysis B, the blocks are:

- binary position/signature presence;
- continuous architecture/composition variables.

The exact current weighting scheme is defined in each script’s `gower_mixed_distance()` function and should remain fixed across all methods for comparability.

---

## Clustering algorithms

Each method and each analysis uses three complementary clusterers.

| Algorithm | Input | Role |
|---|---|---|
| PAM / k-medoids | Precomputed Gower-style distance matrix | Primary partitioning method; yields representative medoid proteins |
| Average-linkage hierarchical clustering | Same precomputed Gower-style distance matrix | Shows related or nested partitions at the same candidate values of `k` |
| HDBSCAN | Robust-scaled numeric/encoded feature matrix | Detects compact density-supported groups and can label mixed proteins as noise |

### PAM initialization

The global script uses the installed `kmedoids` package, whose valid initialization values include `build`, `random`, and `first`. It therefore uses:

```python
init="build"
```

The carrier script currently imports `KMedoids` from `sklearn_extra`; that package supports `k-medoids++`. The package-specific initialization setting must match the imported implementation.

### Hierarchical input requirement

Hierarchical clustering must receive the condensed form of the precomputed distance matrix:

```python
from scipy.spatial.distance import squareform

condensed_distance = squareform(distance, checks=False)
Z = linkage(condensed_distance, method="average")
```

Passing a square distance matrix directly to `linkage()` is incorrect because it may be interpreted as an observation-by-feature matrix rather than as pairwise distances.

### HDBSCAN output

HDBSCAN does not use a fixed number of clusters.

Its outputs are:

```text
cluster_hdbscan
hdbscan_confidence
```

`cluster_hdbscan = -1` denotes an observation classified as noise or unresolved. This is a valid outcome and must not be forced into another cluster.

The confidence output should use the one-dimensional native vector:

```python
out_df["hdbscan_confidence"] = clusterer.probabilities_
```

Do not call:

```python
membership.max(axis=1)
```

when `membership` already equals `clusterer.probabilities_`; it is one-dimensional and causes an axis error.

---

## Candidate cluster counts

PAM and hierarchical clustering are run for the same fixed candidate range:

```python
K_VALUES = [2, 3, 4, 5, 6]
```

No single value of `k` is assumed to be biologically correct.

| Candidate `k` | Initial role |
|---:|---|
| 2 | Broad split, such as LCR-poor versus LCR-rich phenotype |
| 3 | First coarse multi-phenotype solution |
| 4 | Moderate resolution without excessive fragmentation |
| 5 | Tests whether a broad phenotype contains a meaningful subtype |
| 6 | Upper initial bound before likely over-fragmentation |

Values 7–9 are intentionally excluded from the first version. They can be considered later only if the results at `k = 6` still contain large, internally heterogeneous groups.

HDBSCAN is run once with its current `min_cluster_size` and `min_samples` settings. A parameter sweep is deferred to the later stability and validation phase.

---

## Per-method execution workflow

For each detection method:

1. Filter annotation and feature rows to proteins in the RBP universe.
2. Construct either the all-RBP protein table (Analysis A) or carrier-only protein table (Analysis B).
3. Add RNA-target superclass as metadata only.
4. Build the method-specific feature matrix.
5. Remove zero-variance features, if present.
6. Compute one Gower-style distance matrix for PAM and hierarchical clustering.
7. Run PAM for every value in `K_VALUES`.
8. Run hierarchical clustering for every value in `K_VALUES` using the same distance matrix.
9. Run HDBSCAN.
10. Save protein-level assignments and algorithm-specific cluster summaries.

No labels, enrichment tests, external-background comparisons, or cross-method comparisons are used during these ten steps.

---

## Output structure

Each method and analysis should have its own directory.

```text
lcr_analyses/clustering/
├── global/
│   ├── AlcoR/
│   ├── CAST/
│   ├── FLPS/
│   ├── LCRFinder/
│   ├── SEG/
│   └── SEG_intermediate/
└── carriers/
    ├── AlcoR/
    ├── CAST/
    ├── FLPS/
    ├── LCRFinder/
    ├── SEG/
    └── SEG_intermediate/
```

### Protein assignment table

Each method directory contains one `protein_clusters.csv` file with one protein per row.

Expected columns:

```text
uniprot_accession
uniprot_length
rnaprimaryclass
method
cluster_pam_k2
cluster_pam_k3
cluster_pam_k4
cluster_pam_k5
cluster_pam_k6
cluster_hierarchical_k2
cluster_hierarchical_k3
cluster_hierarchical_k4
cluster_hierarchical_k5
cluster_hierarchical_k6
cluster_hdbscan
hdbscan_confidence
```

RNA class and protein length are output metadata. They are not clustering features.

### Cluster-summary tables

For each PAM and hierarchical run, save a separate summary:

```text
cluster_summary_pam_k2.csv
cluster_summary_pam_k3.csv
cluster_summary_pam_k4.csv
cluster_summary_pam_k5.csv
cluster_summary_pam_k6.csv

cluster_summary_hierarchical_k2.csv
cluster_summary_hierarchical_k3.csv
cluster_summary_hierarchical_k4.csv
cluster_summary_hierarchical_k5.csv
cluster_summary_hierarchical_k6.csv
```

Summaries contain cluster size and feature profiles:

- architecture medians;
- position-presence proportions;
- signature-presence proportions;
- carrier-analysis composition medians;
- RNA-class counts as descriptive metadata.

### Manifest

The top-level `method_manifest.csv` records one row per:

```text
method × analysis × algorithm × k
```

For HDBSCAN, record method, analysis, protein count, selected cluster count, and noise count/fraction when this summary is later added.

---

## Deferred work

The following are deliberately outside the present workflow-building phase:

- Choosing a final preferred value of `k`.
- Silhouette or other internal partition diagnostics.
- HDBSCAN parameter sweeps.
- Bootstrap/subsampling stability analysis.
- Cluster agreement across PAM, hierarchical, and HDBSCAN.
- Cross-method cluster comparison.
- RNA-class enrichment/depletion tests.
- Protein-length distribution checks across final clusters.
- External non-RBP proteome comparison.
- Biological naming of final clusters.
- Visualization choices, including heatmaps, UMAP, or PCA.

These steps should be implemented only after the method-specific clustering outputs are generated and inspected.

---

## Practical interpretation guardrails

1. A cluster is an unsupervised LCR phenotype, not automatically an RNA-target class.
2. RNA-class counts in cluster summaries are descriptive metadata only at this stage.
3. A dense HDBSCAN core can be useful, but a large HDBSCAN noise fraction is also an informative result.
4. PAM and hierarchical clustering are intentionally run across several `k` values; no one output is privileged before later evaluation.
5. Do not compare absolute coverage, LCR count, or cluster number across callers before the separate comparative phase. Detection methods define different LCR populations.
6. Do not treat an all-zero signature feature as evidence of biological absence; it may be absent because the current detection/annotation method did not produce it.
7. Keep all protein identifiers and raw output metadata so cluster membership can later be traced back to LCR coordinates, domains, and sequences.

---
