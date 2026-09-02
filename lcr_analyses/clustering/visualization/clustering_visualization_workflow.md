# Phase 3.2.2 — Clustering visualization and results

**Project:** Mapping RNA binders
**Status:** carrier-only PAM (k=5) visualization completed for all methods; global results reserved.
**Updated:** 2026-09-02

## Purpose

This phase visualizes the method-specific unsupervised clustering of RNA-binding proteins (RBPs) by their low-complexity-region (LCR) phenotypes.

The visualization workflow has two sequential biological questions:

1. **What does each cluster represent biologically?**
2. **How are RNA-target superclasses distributed across the discovered LCR phenotypes?**

The goal is to make clustering results interpretable without conflating an unsupervised LCR phenotype with an RNA-target functional class.

## Scope and principles

- Clustering remains **method-specific**. AlcoR, CAST, fLPS, LCRFinder, SEG, and SEG_intermediate are visualized independently.
- The protein remains the biological and statistical unit.
- PAM is the primary partitioning solution for representative-protein visualization because it provides medoids that are real proteins.
- RNA-target superclass is metadata for interpretation. It is not a clustering feature and must not be used to choose a clustering solution.
- All cluster membership, class-distribution, and enrichment comparisons are restricted to the same method × analysis population.
- Figures should be clean, minimal, and auditable. Numerical values shown in figures must be saved in companion tables.

## Selected clustering solutions

The initial visualization workflow uses one selected PAM solution per analysis:

| Analysis | Population | Selected PAM solution | Main purpose |
|---|---:|---:|---|
| Analysis A: global LCR phenotype | All RBPs, including non-carriers | PAM k = 4 | Broad RBPome LCR architecture |
| Analysis B: carrier phenotype | Proteins carrying at least one LCR for the current method | PAM k = 5 | Detailed LCR positional, signature, and composition phenotypes |

These selected resolutions are visualization targets. They do not prevent later stability, agreement, or alternative-resolution analyses.

## Question 1: biological meaning of clusters

### Interpretation target

A cluster represents a recurring **LCR phenotype**: a group of proteins with similar LCR architecture, positional context, primary physicochemical signature, and, in the carrier analysis, LCR composition.

A cluster is not automatically:

- an RNA-target superclass;
- a protein family;
- a pathway;
- a cellular compartment; or
- a mechanistic functional module.

Cluster labels should initially be phenotype descriptions, for example:

- `high-coverage cation-pi-rich LCR proteins`;
- `domain-adjacent acidic LCR proteins`;
- `terminal aromatic-sticker LCR proteins`; or
- `basic arginine-rich LCR proteins`.

### Cluster cards

The primary visualization is one compact **cluster card** per PAM cluster.

Each card combines a cluster-level profile with one representative real protein: the PAM medoid.

#### Cluster-level profile

Report the following values for the entire cluster:

| Quantity | Summary | Rationale |
|---|---|---|
| Cluster size | Protein count, `n` | Records the denominator and prevents overinterpretation of small clusters |
| LCR coverage | Median protein-level coverage | Describes the amount of protein sequence occupied by detected LCRs |
| Main domain-position feature | Highest protein-level presence fraction and count | Identifies the most common LCR positional context |
| Second domain-position feature | Second-highest protein-level presence fraction and count | Preserves a major secondary positional pattern |
| Main primary signature | Highest protein-level presence fraction and count | Identifies the dominant LCR signature phenotype |
| Second primary signature | Second-highest protein-level presence fraction and count | Preserves a major secondary signature pattern |
| NCPR, Analysis B only | Median length-weighted protein NCPR | Describes net LCR charge direction and magnitude |
| FCR, Analysis B only | Median length-weighted protein FCR | Describes total charged-residue content |

Position and signature values are computed at protein level. A protein contributes at most once to a given position or signature presence category, even if it has multiple LCRs of that type.

For Analysis B, protein-level composition is aggregated from LCR-level measurements with LCR length as the weight:

\[
\mathrm{protein\ property} =
\frac{\sum_{\ell=1}^{n_{\mathrm{LCR}}}
\left(\mathrm{LCR\ length}_{\ell} \times \mathrm{property}_{\ell}\right)}
{\sum_{\ell=1}^{n_{\mathrm{LCR}}}\mathrm{LCR\ length}_{\ell}}.
\]

#### Profile-bar encoding

- Coverage, position presence, signature presence, and FCR are shown as left-to-right horizontal bars with their exact numerical values alongside.
- Position and signature bars show both percentage and protein count, for example `55.1% (130/236)`.
- NCPR uses a rightward bar whose length represents `abs(NCPR)`.
- NCPR colour encodes sign:
  - positive NCPR: muted red;
  - negative NCPR: muted yellow;
  - zero: grey.
- The signed numerical NCPR value remains visible, for example `+0.184` or `−0.143`.

This avoids a visually confusing backward bar while retaining both charge magnitude and direction.

#### PAM medoid

The PAM medoid is an actual protein within a cluster that minimizes summed Gower distance to all other members of that cluster:

\[
\mathrm{medoid}(c) =
\underset{i \in c}{\operatorname{argmin}}
\sum_{j \in c}D_{ij}.
\]

The medoid is an illustrative representative, not a claim that every cluster member has identical protein architecture or molecular function.

Medoid calculation must use the same feature construction and Gower-style mixed distance used by the original PAM clustering solution:

- binary domain-position and signature presence;
- categorical multiplicity in Analysis A;
- continuous architecture and composition features;
- the same scaling and block-weighting rules;
- zero-variance feature removal.

#### Protein-coordinate architecture panel

The lower section of each card shows the medoid protein as a linear coordinate map.

| Visual element | Meaning |
|---|---|
| Main horizontal backbone | Full UniProt protein length |
| Domain track | Pfam domain blocks with name and start–end coordinates |
| LCR track | All detected LCR intervals on the same protein-coordinate line |
| LCR colour | Primary physicochemical signature |
| LCR label | Signature and start–end coordinates |
| Domain label | Pfam name and start–end coordinates |

All LCRs for a medoid are drawn on the same backbone coordinate line. This preserves their actual relative locations and avoids creating an artificial separate row for proteins with several LCRs.

A medoid without retained Pfam coordinates remains valid. Its card shows the protein backbone and LCR intervals but no domain blocks. This means no retained Pfam annotation was available in the plotting input; it does not establish that the protein lacks a folded domain.

### Cluster-card inputs and outputs

Required inputs:

```text
protein_clusters.csv
lcr_annotations.xlsx
lcr_features.xlsx, composition sheet
pfam_lcr_overlap.xlsx, Pfam_hits sheet
```

Core inputs include:

```text
protein clusters:
uniprot_accession, uniprot_length, cluster_pam_k4 / cluster_pam_k5

LCR annotations:
protein_id, source_method, start, end, length,
domain_position_class, primary_physicochemical_annotation

LCR composition:
protein_id, source_method, start, end, length,
physicochemical property metrics including fcr and ncpr

Pfam coordinates:
protein_id, pfam_name, domain_start, domain_end
```

Protein identifiers are standardized to canonical UniProt accessions before joining. Composite IDs are reduced to the text before the first `|` character.

Outputs per method × analysis:

```text
cluster_card_<method>_<analysis>_<pam_solution>_cluster_<id>.png
cluster_card_summary_<method>_<analysis>_<pam_solution>.csv
reconstructed_features_<method>_<analysis>_<pam_solution>.csv
```

The summary table records every quantity displayed in the card, including medoid identity, cluster size, profile values, and selected top position/signature features.

## Question 2: RNA-class distribution across clusters

### Rationale

RNA-target classes are highly unbalanced, particularly because mRNA and unknown classes are large. A plot showing the RNA-class composition of each cluster, \(P(\mathrm{class}\mid\mathrm{cluster})\), is therefore visually dominated by the largest classes and is not the first descriptive view.

The first question should instead be:

> Within a given RNA-target class, how are proteins distributed across the discovered LCR phenotypes?

### First descriptive figure

Use a vertical 100% stacked bar chart of:

\[
P(\mathrm{cluster}\mid\mathrm{RNA\ class}) =
\frac{n_{\mathrm{class,cluster}}}{n_{\mathrm{class}}}.
\]

Figure design:

- One vertical bar per RNA-target superclass.
- Every bar sums to 100% independently.
- Each coloured segment represents a PAM cluster.
- Cluster colours remain fixed within one method × analysis plot.
- The class label includes its protein count, for example `mRNA (n=487)`.
- Segment percentages are printed only above a minimum display threshold to avoid crowding.
- No p-values or enrichment language are used at this stage.

This answers whether a class is concentrated in one or two LCR phenotypes or spread broadly across the available phenotypes.

### Class-distribution summary table

For each RNA-target superclass, record:

| Field | Meaning |
|---|---|
| `rna_class` | RNA-target superclass |
| `n_class` | Number of proteins in the method × analysis population |
| `n_clusters_total` | Number of PAM clusters in the selected solution |
| `representation_threshold` | Minimum \(P(\mathrm{cluster}\mid\mathrm{class})\) treated as represented |
| `n_clusters_represented` | Number of clusters at or above the threshold |
| `clusters_represented` | Cluster IDs and conditional percentages |
| `dominant_cluster` | Cluster with the greatest class contribution |
| `dominant_cluster_n` | Number of class proteins in that cluster |
| `dominant_cluster_fraction` | \(P(\mathrm{dominant\ cluster}\mid\mathrm{class})\) |

The initial representation threshold is:

\[
P(\mathrm{cluster}\mid\mathrm{class}) \geq 0.05.
\]

This threshold is descriptive and configurable. Small RNA classes must always be interpreted with their reported denominator.

Outputs per method × analysis:

```text
<method>_<analysis>_<pam_solution>_rna_class_distribution.png
<method>_<analysis>_<pam_solution>_rna_class_distribution_counts.csv
<method>_<analysis>_<pam_solution>_rna_class_distribution_fractions.csv
<method>_<analysis>_<pam_solution>_rna_class_distribution_summary.csv
```

## Amount-normalized cluster-class comparison

After examining the class-normalized distribution, assess the reverse relationship—RNA-class contribution to each cluster—using abundance normalization.

For each cluster × RNA-class cell, calculate the expected count under independence:

\[
E_{\mathrm{class,cluster}} =
\frac{n_{\mathrm{class}} \times n_{\mathrm{cluster}}}{N}.
\]

Then calculate the descriptive observed/expected ratio:

\[
O/E =
\frac{n_{\mathrm{class,cluster}}}
{E_{\mathrm{class,cluster}}}
=
\frac{P(\mathrm{class}\mid\mathrm{cluster})}
{P(\mathrm{class})}.
\]

Visualize this quantity as:

\[
\log_2(O/E).
\]

Interpretation:

| Value | Meaning |
|---:|---|
| `0` | Observed representation equals the population expectation |
| `+1` | Twofold over-representation |
| `−1` | Twofold under-representation |

The appropriate plot is a cluster × RNA-class heatmap:

- rows: RNA-target superclasses;
- columns: clusters;
- colour: `log2(O/E)`, centred at zero;
- optional text: observed count;
- visibly mark sparse cells, such as cells with observed or expected count below 5.

This is still descriptive. It corrects the visual bias caused by unequal class sizes but does not establish statistical enrichment.

## Later inferential extension

Formal enrichment testing is deferred until after descriptive distribution and normalized heatmap inspection.

For each cluster × RNA-class comparison, construct a cluster-versus-rest 2×2 table and calculate:

- Fisher exact-test p-value;
- Haldane-corrected log2 odds ratio;
- confidence interval for the odds ratio;
- BH-FDR across a predeclared cluster × class test family;
- a small-sample `descriptive_only` flag.

The final inferential heatmap should use the signed log2 odds ratio as colour and FDR status as a separate visual annotation. Raw p-values should not be the main colour scale.

## Completed results

### CAST carriers — PAM (k=5)

**Carrier population:** n=1,055 proteins.
**Primary outputs:** `cluster_card_summary_CAST_carriers_cluster_pam_k5.csv`, `CAST_carriers_cluster_pam_k5_rna_class_distribution_summary.csv`, and the corresponding Gower-distance classical-MDS plot.

| Cluster | n | Medoid | Median coverage | Median NCPR | Median FCR | Dominant position | Dominant signature | Provisional phenotype description |
|---|---:|---|---:|---:|---:|---|---|---|
| C0 | 288 | A0ACI8U7J9 | 0.408 | +0.006 | 0.285 | Distal-terminal, 33.3% | Cation–π-rich neighborhood, 71.9% | High-coverage cation–π-rich mixed-context LCR phenotype |
| C1 | 236 | O15541 | 0.252 | −0.143 | 0.459 | Domain-adjacent, 55.1% | Acidic region, 88.1% | Domain-adjacent acidic LCR phenotype |
| C2 | 173 | Q8WWH5 | 0.204 | +0.014 | 0.134 | Domain-adjacent, 66.5% | Disorder-rich spacer, 70.5% | Domain-adjacent disorder-rich / RG-associated phenotype |
| C3 | 203 | Q8IU60 | 0.322 | 0.000 | 0.250 | Distal-terminal, 78.8% | Aromatic-sticker region, 67.0% | Distal-terminal aromatic-sticker phenotype |
| C4 | 155 | Q13601 | 0.262 | +0.184 | 0.460 | Domain-edge, 35.5% | Basic-enriched region, 56.8% | Basic / arginine-rich, domain-associated phenotype |

#### CAST MDS display

The accompanying scatter plot is a two-dimensional classical-MDS projection of the PAM Gower distance. Classical MDS 1 and 2 explain 12.2% and 11.0% of positive-eigenvalue variance, respectively (23.2% combined). The displayed overlap therefore reflects dimensional projection of a higher-dimensional mixed-feature distance and should not be interpreted as the full cluster separation.

The acidic C1 and basic C4 groups are visibly displaced from the central mixed cloud. C2 is positioned toward higher MDS2 values, C3 has an upper-right concentration, and C0 occupies a broad mixed/intermediate region. CAST therefore captures interpretable phenotype poles together with a more continuous carrier landscape.

#### CAST RNA-class distribution

Distributions are reported as \(P(\mathrm{cluster}\mid\mathrm{RNA\ class})\) within the CAST carrier population.

| RNA class | Carrier n | Dominant cluster | Fraction in dominant cluster | Interpretation |
|---|---:|---|---:|---|
| mRNA | 487 | C0 | 29.0% | Distributed across all five phenotypes; C3 (23.2%), C2 (20.9%), and C1 (19.1%) are also represented |
| tRNA | 51 | C0 | 27.5% | Broad distribution across all five phenotypes |
| pre-rRNA | 82 | C1 | 37.8% | Mainly split between acidic C1 (37.8%) and basic C4 (36.6%) |
| snRNA | 60 | C1 | 33.3% | All five phenotypes represented |
| snoRNA | 29 | C1 | 34.5% | Acidic C1 is the largest component; interpret with small n |
| ncRNA | 69 | C0 | 31.9% | Broad distribution across all five phenotypes |
| ribosomal protein | 46 | C4 | 69.6% | Strong concentration in the basic / arginine-rich phenotype |
| diverse | 15 | C0 | 46.7% | Small carrier subset |
| unknown | 216 | C0 | 32.4% | Broad distribution across all five phenotypes |

### Intermediate-SEG carriers — PAM (k=5)

**Carrier population:** n=622 proteins.
**Primary outputs:** `cluster_card_summary_SEG_intermediate_carriers_cluster_pam_k5.csv` and `SEG_intermediate_carriers_cluster_pam_k5_rna_class_distribution_summary.csv`.

| Cluster | n | Medoid | Median coverage | Median NCPR | Median FCR | Dominant position | Dominant signature | Provisional phenotype description |
|---|---:|---|---:|---:|---:|---|---|---|
| C0 | 250 | Q8IWZ8 | 0.065 | 0.000 | 0.089 | Distal-terminal, 62.8% | Disorder-rich spacer, 79.6% | Distal-terminal disorder-rich-spacer phenotype |
| C1 | 119 | A6NDY0 | 0.048 | −0.516 | 0.643 | Distal-terminal, 54.6% | Acidic region, 99.2% | Strongly acidic LCR phenotype |
| C2 | 99 | Q8WW36 | 0.147 | +0.099 | 0.179 | Domain-adjacent, 65.7% | RG/RGG-repeat region, 66.7% | Domain-adjacent RG/RGG-repeat phenotype |
| C3 | 80 | Q8WXF0 | 0.184 | +0.261 | 0.433 | Domain-adjacent, 60.0% | SR/RS-repeat region, 87.5% | Domain-adjacent SR/RS-repeat phenotype |
| C4 | 74 | Q9UI10 | 0.056 | +0.239 | 0.658 | Distal-terminal, 66.2% | Basic-enriched region, 67.6% | Basic / arginine-rich LCR phenotype |

The intermediate-SEG partition resolves a particularly direct five-phenotype map: disorder-rich spacer, acidic, RG/RGG-repeat, SR/RS-repeat, and basic/arginine-rich. Under the present method-specific feature and distance representation, its cluster identities are more signature-concentrated than the corresponding CAST solution.

#### Intermediate-SEG RNA-class distribution

Distributions are reported as \(P(\mathrm{cluster}\mid\mathrm{RNA\ class})\) within the intermediate-SEG carrier population.

| RNA class | Carrier n | Dominant cluster | Fraction in dominant cluster | Interpretation |
|---|---:|---|---:|---|
| mRNA | 330 | C0 | 43.6% | Broad disorder-spacer component with substantial RG/RGG C2 (17.9%) and SR/RS C3 (19.4%) components |
| tRNA | 17 | C0 | 41.2% | Small carrier subset; C4 also represents 35.3% |
| pre-rRNA | 39 | C1 | 74.4% | Strong concentration in the acidic phenotype |
| snRNA | 39 | C0 | 38.5% | All five phenotypes represented |
| snoRNA | 18 | C1 | 44.4% | Acidic phenotype is the largest component; interpret with small n |
| ncRNA | 36 | C0 | 44.4% | Mainly disorder-spacer and acidic components |
| ribosomal protein | 16 | C4 | 37.5% | Largest component is basic; small carrier subset |
| diverse | 9 | C0 | 66.7% | Small carrier subset |
| unknown | 118 | C0 | 44.1% | Broad distribution with an acidic C1 component (21.2%) |

### AlcoR carriers — PAM (k=5)

**Carrier population:** n=968 proteins.
**Primary outputs:** `cluster_card_summary_AlcoR_carriers_cluster_pam_k5.csv`, `AlcoR_carriers_cluster_pam_k5_rna_class_distribution_summary.csv`, and the corresponding Gower-distance classical-MDS plot.

| Cluster | n | Medoid | Median coverage | Median NCPR | Median FCR | Dominant position | Dominant signature | Provisional phenotype description |
|---|---:|---|---:|---:|---:|---|---|---|
| C0 | 235 | A6NHQ2 | 0.265 | +0.015 | 0.250 | Domain-edge, 100.0% | Aromatic-sticker region, 72.3% | Domain-edge aromatic-sticker phenotype |
| C1 | 269 | Q9BXT6 | 0.119 | 0.000 | 0.238 | Distal-terminal, 27.9% | Aromatic-sticker region, 72.9% | Low-coverage distal-terminal aromatic-sticker phenotype |
| C2 | 253 | Q9H9J2 | 0.946 | +0.008 | 0.263 | Domain-edge, 87.7% | Aromatic-sticker region, 73.5% | Near-complete aromatic-sticker domain-edge phenotype |
| C3 | 84 | J3QT46 | 0.287 | 0.000 | 0.249 | Distal-terminal, 86.9% | Aromatic-sticker region, 92.9% | Distal-terminal aromatic-sticker enriched phenotype |
| C4 | 127 | M0R2C6 | 0.091 | 0.000 | 0.333 | Distal-terminal, 71.7% | Mixed-charge polyampholyte, 47.2% | Low-coverage mixed-charge polyampholyte phenotype |

AlcoR resolves a heterogeneous carrier landscape dominated by aromatic-sticker signatures across all five clusters. The clusters primarily differ in domain-position architecture and LCR coverage rather than physicochemical signature, reflecting AlcoR's tendency to detect longer, more heterogeneous low-complexity segments.

#### AlcoR MDS display

The accompanying scatter plot is a two-dimensional classical-MDS projection of the PAM Gower distance. Classical MDS 1 and 2 explain 14.6% and 11.3% of positive-eigenvalue variance, respectively (25.9% combined). The displayed overlap therefore reflects dimensional projection of a higher-dimensional mixed-feature distance and should not be interpreted as the full cluster separation.

The clusters show substantial mutual overlap, consistent with the shared aromatic-sticker signature dominance. C2, the near-complete coverage phenotype, occupies a distinct region toward higher MDS1 values. The remaining clusters form a more continuous distribution, suggesting that AlcoR's carrier phenotypes represent a gradient of LCR coverage and positional context rather than sharply separated groups.

#### AlcoR RNA-class distribution

Distributions are reported as \(P(\mathrm{cluster}\mid\mathrm{RNA\ class})\) within the AlcoR carrier population.

| RNA class | Carrier n | Dominant cluster | Fraction in dominant cluster | Interpretation |
|---|---:|---|---:|---|
| mRNA | 377 | C1 | 29.2% | Broad distribution across all five phenotypes |
| tRNA | 85 | C0 | 30.6% | Broad distribution across all five phenotypes |
| pre-rRNA | 53 | C2 | 30.2% | Moderate concentration in the near-complete coverage phenotype |
| snRNA | 41 | C1 | 36.6% | All five phenotypes represented |
| snoRNA | 20 | C0 | 30.0% | Broad distribution; interpret with small n |
| ncRNA | 64 | C1 | 40.6% | Moderate concentration in the low-coverage distal-terminal phenotype |
| ribosomal protein | 52 | C2 | 36.5% | Moderate concentration in the near-complete coverage phenotype |
| diverse | 22 | C0 | 27.3% | Broad distribution; small carrier subset |
| unknown | 254 | C2 | 33.1% | Broad distribution across all five phenotypes |

### fLPS carriers — PAM (k=5)

**Carrier population:** n=1,697 proteins.
**Primary outputs:** `cluster_card_summary_FLPS_carriers_cluster_pam_k5.csv`, `FLPS_carriers_cluster_pam_k5_rna_class_distribution_summary.csv`, and the corresponding Gower-distance classical-MDS plot.

| Cluster | n | Medoid | Median coverage | Median NCPR | Median FCR | Dominant position | Dominant signature | Provisional phenotype description |
|---|---:|---|---:|---:|---:|---|---|---|
| C0 | 417 | Q9BQ04 | 0.219 | +0.014 | 0.159 | Distal-terminal, 79.4% | Disorder-rich spacer, 85.1% | Distal-terminal disorder-rich-spacer phenotype |
| C1 | 368 | Q6P6C2 | 0.214 | −0.002 | 0.351 | Distal-terminal, 84.0% | Acidic region, 85.9% | Distal-terminal acidic phenotype |
| C2 | 301 | Q96DP5 | 0.129 | +0.261 | 0.435 | Domain-edge, 36.5% | Arginine-rich RNA-contact region, 61.8% | Basic / arginine-rich RNA-contact phenotype |
| C3 | 338 | Q9BZE1 | 0.104 | +0.020 | 0.172 | Domain-edge, 30.8% | Aromatic-sticker region, 59.8% | Low-coverage aromatic-sticker phenotype |
| C4 | 273 | R4GNH1 | 0.137 | −0.142 | 0.464 | Distal-terminal, 30.8% | Acidic region, 88.6% | Strongly acidic domain-associated phenotype |

fLPS resolves a high-coverage carrier landscape with five broadly distributed phenotypes. The partition separates disorder-rich, acidic, basic/arginine-rich, and aromatic-sticker groups, with acidic and basic clusters forming the most compositionally distinct poles.

#### fLPS MDS display

The accompanying scatter plot is a two-dimensional classical-MDS projection of the PAM Gower distance. Classical MDS 1 and 2 explain 14.7% and 10.9% of positive-eigenvalue variance, respectively (25.6% combined). The displayed overlap therefore reflects dimensional projection of a higher-dimensional mixed-feature distance and should not be interpreted as the full cluster separation.

The basic C2 and strongly acidic C4 groups are visibly displaced from the central cloud. C0 and C1, the disorder-rich and acidic groups, show substantial mutual overlap, consistent with their shared distal-terminal positional context. C3 occupies a lower-left region. fLPS therefore recovers interpretable phenotype poles alongside a more continuous carrier distribution.

#### fLPS RNA-class distribution

Distributions are reported as \(P(\mathrm{cluster}\mid\mathrm{RNA\ class})\) within the fLPS carrier population.

| RNA class | Carrier n | Dominant cluster | Fraction in dominant cluster | Interpretation |
|---|---:|---|---:|---|
| mRNA | 674 | C0 | 35.5% | Broad distribution across all five phenotypes |
| tRNA | 145 | C3 | 36.6% | Moderate concentration in the aromatic-sticker phenotype |
| pre-rRNA | 109 | C1 | 33.9% | Concentrated in the distal-terminal acidic phenotype |
| snRNA | 92 | C4 | 28.3% | All five phenotypes represented |
| snoRNA | 37 | C4 | 37.8% | Moderate concentration in the strongly acidic phenotype; interpret with small n |
| ncRNA | 113 | C3 | 23.9% | Broad distribution across all five phenotypes |
| ribosomal protein | 118 | C2 | 59.3% | Strong concentration in the basic / arginine-rich phenotype |
| diverse | 43 | C3 | 34.9% | Broad distribution |
| unknown | 366 | C1 | 25.7% | Broad distribution across all five phenotypes |

### LCRFinder carriers — PAM (k=5)

**Carrier population:** n=1,008 proteins.
**Primary outputs:** `cluster_card_summary_LCRFinder_carriers_cluster_pam_k5.csv`, `LCRFinder_carriers_cluster_pam_k5_rna_class_distribution_summary.csv`, and the corresponding Gower-distance classical-MDS plot.

| Cluster | n | Medoid | Median coverage | Median NCPR | Median FCR | Dominant position | Dominant signature | Provisional phenotype description |
|---|---:|---|---:|---:|---:|---|---|---|
| C0 | 272 | Q6XE24 | 0.058 | 0.000 | 0.050 | Distal-terminal, 66.5% | Disorder-rich spacer, 90.4% | Distal-terminal disorder-rich-spacer phenotype |
| C1 | 281 | Q03468 | 0.062 | −0.042 | 0.511 | Distal-terminal, 76.5% | Acidic region, 94.3% | Distal-terminal acidic phenotype |
| C2 | 103 | O95373 | 0.029 | −0.612 | 0.800 | Distal-terminal, 62.1% | Acidic region, 98.1% | Strongly acidic low-coverage phenotype |
| C3 | 146 | Q96J94 | 0.038 | +0.420 | 0.500 | Distal-terminal, 52.1% | Arginine-rich RNA-contact region, 71.9% | Basic / arginine-rich RNA-contact phenotype |
| C4 | 206 | Q5T5J6 | 0.030 | 0.000 | 0.076 | Distal-terminal, 56.8% | Hydrophobic region, 86.9% | Distal-terminal hydrophobic phenotype |

LCRFinder resolves a low-coverage carrier landscape with signature-concentrated phenotypes. The partition separates disorder-rich, acidic (two clusters), basic/arginine-rich, and hydrophobic groups. The presence of two acidic clusters (C1 and C2) distinguished by NCPR magnitude suggests LCRFinder detects a continuous acidic spectrum rather than a single acidic phenotype.

#### LCRFinder MDS display

The accompanying scatter plot is a two-dimensional classical-MDS projection of the PAM Gower distance. Classical MDS 1 and 2 explain 16.5% and 12.5% of positive-eigenvalue variance, respectively (29.0% combined). The displayed overlap therefore reflects dimensional projection of a higher-dimensional mixed-feature distance and should not be interpreted as the full cluster separation.

The strongly acidic C2 and basic C3 groups are visibly displaced from the central cloud. C0, C1, and C4 show more mutual overlap, consistent with their shared low-coverage and distal-terminal positional context. LCRFinder therefore recovers interpretable phenotype poles with a more continuous intermediate distribution.

#### LCRFinder RNA-class distribution

Distributions are reported as \(P(\mathrm{cluster}\mid\mathrm{RNA\ class})\) within the LCRFinder carrier population.

| RNA class | Carrier n | Dominant cluster | Fraction in dominant cluster | Interpretation |
|---|---:|---|---:|---|
| mRNA | 465 | C0 | 34.4% | Broad distribution across all five phenotypes |
| tRNA | 51 | C4 | 37.3% | Moderate concentration in the hydrophobic phenotype |
| pre-rRNA | 67 | C1 | 32.8% | Split between acidic C1 (32.8%) and strongly acidic C2 (32.8%) |
| snRNA | 58 | C1 | 37.9% | Moderate concentration in the acidic phenotype |
| snoRNA | 29 | C1 | 51.7% | Strong concentration in the acidic phenotype; interpret with small n |
| ncRNA | 65 | C0 | 29.2% | Broad distribution across all five phenotypes |
| ribosomal protein | 21 | C3 | 52.4% | Strong concentration in the basic / arginine-rich phenotype; small carrier subset |
| diverse | 21 | C4 | 52.4% | Strong concentration in the hydrophobic phenotype; small carrier subset |
| unknown | 231 | C1 | 30.3% | Broad distribution across all five phenotypes |

### Strict-SEG carriers — PAM (k=5)

**Carrier population:** n=296 proteins.
**Primary outputs:** `cluster_card_summary_SEG_carriers_cluster_pam_k5.csv`, `SEG_carriers_cluster_pam_k5_rna_class_distribution_summary.csv`, and the corresponding Gower-distance classical-MDS plot.

| Cluster | n | Medoid | Median coverage | Median NCPR | Median FCR | Dominant position | Dominant signature | Provisional phenotype description |
|---|---:|---|---:|---:|---:|---|---|---|
| C0 | 75 | H0YDK8 | 0.040 | 0.000 | 0.040 | Distal-terminal, 86.7% | Disorder-rich spacer, 81.3% | Distal-terminal disorder-rich-spacer phenotype |
| C1 | 65 | Q15029 | 0.031 | −0.656 | 0.763 | Distal-terminal, 55.4% | Acidic region, 100.0% | Strongly acidic LCR phenotype |
| C2 | 44 | Q99729 | 0.106 | +0.126 | 0.196 | Domain-adjacent, 81.8% | RG/RGG-repeat region, 54.5% | Domain-adjacent RG/RGG-repeat phenotype |
| C3 | 56 | Q9Y2U8 | 0.035 | 0.000 | 0.000 | Domain-adjacent, 33.9% | Disorder-rich spacer, 80.4% | Low-coverage neutral disorder phenotype |
| C4 | 56 | Q8NAV1 | 0.082 | +0.416 | 0.522 | Distal-terminal, 62.5% | SR/RS-repeat region, 62.5% | SR/RS-repeat phenotype |

Strict-SEG resolves a small but signature-concentrated carrier population. The five clusters separate disorder-rich, acidic, RG/RGG-repeat, neutral disorder, and SR/RS-repeat phenotypes. The partition is more sparsely populated than other methods, reflecting SEG's strict detection thresholds.

#### Strict-SEG MDS display

The accompanying scatter plot is a two-dimensional classical-MDS projection of the PAM Gower distance. Classical MDS 1 and 2 explain 25.8% and 13.8% of positive-eigenvalue variance, respectively (39.5% combined). The displayed overlap therefore reflects dimensional projection of a higher-dimensional mixed-feature distance and should not be interpreted as the full cluster separation.

The strongly acidic C1 and SR/RS-repeat C4 groups are visibly displaced from the central cloud. C2, the RG/RGG-repeat phenotype, occupies a distinct upper region. C0 and C3, the two disorder-related clusters, show more mutual overlap. Strict-SEG therefore achieves the highest MDS variance explanation among the methods, consistent with its signature-concentrated partition.

#### Strict-SEG RNA-class distribution

Distributions are reported as \(P(\mathrm{cluster}\mid\mathrm{RNA\ class})\) within the strict-SEG carrier population.

| RNA class | Carrier n | Dominant cluster | Fraction in dominant cluster | Interpretation |
|---|---:|---|---:|---|
| mRNA | 170 | C0 | 29.4% | Broad distribution across all five phenotypes |
| tRNA | 3 | C1 | 66.7% | Small carrier subset; C1 only |
| pre-rRNA | 14 | C1 | 92.9% | Strong concentration in the acidic phenotype |
| snRNA | 19 | C4 | 42.1% | All five phenotypes represented |
| snoRNA | 10 | C1 | 60.0% | Moderate concentration in the acidic phenotype; small carrier subset |
| ncRNA | 14 | C3 | 35.7% | Broad distribution across all five phenotypes |
| ribosomal protein | 4 | C3 | 50.0% | Small carrier subset |
| diverse | 2 | C0 | 50.0% | Small carrier subset |
| unknown | 60 | C1 | 30.0% | Broad distribution across all five phenotypes |

## Cross-method synthesis

Across the six methods, acidic and basic/arginine-rich phenotypes recur but are not recovered uniformly. CAST, intermediate SEG, fLPS, and LCRFinder resolve both phenotype types; strict SEG resolves a strongly acidic but no distinct basic/arginine-rich cluster, whereas AlcoR is dominated by aromatic-sticker and coverage/position-defined phenotypes:

| Method | Carrier n | MDS combined variance | Phenotype separation |
|---|---:|---:|---|
| CAST | 1,055 | 23.2% | Broad cation–π-rich phenotype with continuous landscape |
| Intermediate-SEG | 622 | — | Five signature-concentrated phenotypes |
| AlcoR | 968 | 25.9% | Aromatic-sticker dominated with coverage gradients |
| fLPS | 1,697 | 25.6% | Five broadly distributed phenotypes |
| LCRFinder | 1,008 | 29.0% | Low-coverage signature-concentrated phenotypes |
| Strict-SEG | 296 | 39.5% | Small, highly signature-concentrated partition |

The most direct descriptive class-to-phenotype patterns currently observed are:

- pre-rRNA carriers: acidic phenotype, especially strict-SEG C1 (92.9%), intermediate-SEG C1 (74.4%), and LCRFinder C1+C2 (65.6%)
- mRNA carriers: broad disorder-spacer distribution across all methods
- ribosomal-protein carriers: basic / arginine-rich phenotype, particularly CAST C4 (69.6%) and fLPS C2 (59.3%)
- snoRNA carriers: acidic phenotype, especially LCRFinder C1 (51.7%) and strict-SEG C1 (60.0%)

The approaches are complementary and should be interpreted as method-specific views of LCR sequence space rather than as directly interchangeable partitions.

## Summary 
Across the carrier-only PAM k=5 analyses, LCR phenotypes recur as disorder-rich spacer, acidic, repeat-rich, and basic/arginine-rich architectures, but their resolution depends on the detection method. Intermediate SEG provides the clearest separation of sequence-defined acidic, RG/RGG-repeat, SR/RS-repeat, disorder-spacer, and basic phenotypes. CAST and fLPS retain broader, more continuous LCR landscapes, whereas AlcoR is primarily structured by aromatic-sticker content, coverage, and domain position. The most consistent descriptive RNA-class associations are acidic phenotypes among pre-rRNA and snoRNA carriers, broad disorder/repeat phenotypes among mRNA carriers, and basic/arginine-rich phenotypes among ribosomal-protein carriers.

## Interpretation boundaries

1. A cluster is an LCR phenotype, not an RNA class.
2. A medoid is a central example, not a universal biological prototype.
3. RNA-class plots are initially descriptive; do not call a class enriched before formal testing.
4. Always report class and cluster denominators.
5. Do not compare raw cluster number, coverage, or cluster size across LCR callers as if callers defined equivalent LCR populations.
6. Do not combine methods into one figure before a separate cluster-correspondence analysis.
7. Treat sparse classes and small cells as hypothesis-generating only.
8. Keep every output table so figures can be traced to individual proteins, LCR coordinates, domains, and input annotations.
