# LCR physicochemical behavior taxonomy (v1)

## Overview

This document specifies the rule-based physicochemical annotation system for RNA-binding protein (RBP) low-complexity regions (LCRs). Each LCR receives one or more behavior labels derived from three analysis layers:

1. **Composition** (what properties are present)
2. **Distribution** (how properties are arranged: compact, dispersed, periodic)
3. **Co-occurrence** (which properties appear together locally)

These annotations are **hypothesis-oriented descriptors** of likely physical behavior, not direct functional assignments. They are intended to be combined upstream with:

- RNA-target superclass (e.g., mRNA, rRNA, snRNA, etc.)
- Domain-position class (intrinsic, edge, adjacent, linker, external)

## Amino-acid property encoding

Each amino acid is encoded using a 20×¬5 property table based on Taylor/Livingstone/Jalview overlapping property sets, with an independent charge encoding at pH 7.4.

### Properties

- **Charge (pH 7.4)**:  
  - Positive: R, K  
  - Negative: D, E  
  - Neutral: all others (H treated as neutral)
- **Polar**: binary membership in the “polar” set
- **Hydrophobic**: binary membership in the “hydrophobic” set
- **Aromatic**: binary membership in the “aromatic” set (F, Y, W)
- **Small**: binary membership in the “small” set

These are non-exclusive; e.g., K is both polar and hydrophobic, Y is polar, hydrophobic, and aromatic.

## Feature layers

### 1. Composition

For each LCR sequence:

- **AA composition**: fractions of all 20 amino acids.
- **Property fractions**:
  - \(f_{\text{polar}}\)
  - \(f_{\text{hydrophobic}}\)
  - \(f_{\text{aromatic}}\)
  - \(f_{\text{small}}\)
  - \(f_{+}\) (fraction positive)
  - \(f_{-}\) (fraction negative)
- **Charge metrics**:
  - FCR (fraction charged): \(f_{\text{charged}} = f_{+} + f_{-}\)
  - NCPR (net charge per residue): \(f_{+} - f_{-}\)

### 2. Distribution (per property)

For each property track \(X \in \{\text{polar}, \text{hydrophobic}, \text{aromatic}, \text{small}, \text{charged}\}\):

- Compute:
  - Number of runs \(R_X\)
  - Mean run length \(\overline{\ell}_X\)
  - Max run length \(\ell^{\max}_X\)
  - Mean gap \(\overline{d}_X\)
  - Coefficient of variation of gaps \(\mathrm{CV}_X\)

Assign one of three distribution classes:

- **Compact**: few runs, at least one long block
- **Dispersed**: many short runs, no long blocks
- **Periodic-ish**: many occurrences with relatively uniform spacing (low gap CV)
- **Insufficient**: too few occurrences to assign reliably

### 3. Co-occurrence (property-based, no explicit motifs)

Using a sliding window of size \(w\) (default 5):

- For selected property pairs \(X,Y\), count how often both occur within the same window.
- Normalize by number of windows to get a co-occurrence score \(c_{XY}\) in [0,1].

Pairs considered:

- Positive–aromatic (cation–π potential)
- Hydrophobic–hydrophobic (local cohesive patches)
- Polar–polar (polar clusters)
- Polar–hydrophobic (mixed neighborhoods)
- Positive–negative (local charge mixing)

Scores are later classified as low/medium/high using dataset percentiles (e.g., top 25% = “high”).

## Behavior signature taxonomy

Each LCR can receive multiple labels. Rules are applied in order; a label is assigned if its conditions are met.

### Charge-driven signatures

#### A. Basic-enriched patch

- \(f_{+} \geq 0.25\)
- NCPR ≥ 0.10
- Charge distribution: **compact**

**Annotation**: “Basic patch; strong electrostatic interaction potential, likely to interact with negatively charged partners such as RNA backbones or acidic proteins.”

#### B. Acidic patch

- \(f_{-} \geq 0.25\)
- NCPR ≤ −0.10
- Charge distribution: **compact**

**Annotation**: “Acidic patch; strong negative electrostatics, potential docking site for basic regions and metal ions.”

#### C. Mixed-charge polyampholyte

- FCR ≥ 0.30
- |NCPR| ≤ 0.05
- Charge distribution: **dispersed** or **periodic-ish**

**Annotation**: “Mixed-charge polyampholyte; sequence expected to be salt-sensitive with tunable chain dimensions and electrostatic interactions.”

### Polarity / hydropathy signatures

#### D. Polar linker

- \(f_{\text{polar}} \geq 0.50\)
- \(f_{\text{hydrophobic}} \leq 0.20\)
- Polar distribution: **dispersed**

**Annotation**: “Polar, solvent-exposed, likely flexible linker or interaction surface with low intrinsic cohesion.”

#### E. Hydrophobic patch

- \(f_{\text{hydrophobic}} \geq 0.40\)
- Hydrophobic distribution: **compact**

**Annotation**: “Hydrophobic patch; cohesive region with potential to form local cores, interface surfaces, or aggregation-prone segments.”

#### F. Mixed polar–hydrophobic segment

- \(f_{\text{polar}}\) and \(f_{\text{hydrophobic}}\) both in [0.30, 0.50]
- Hydrophobic and polar distribution both **dispersed** or **periodic-ish**

**Annotation**: “Mixed polar–hydrophobic segment; potential amphipathic-like behavior with solvent-exposed and cohesive patches interleaved.”

### Aromatic / sticker signatures

#### G. Aromatic-enriched patch

- \(f_{\text{aromatic}} \geq 0.08\)
- Aromatic distribution: **compact**

**Annotation**: “Aromatic patch; enriched in sticker residues capable of π–π interactions; potential contributor to cohesion and phase behavior.”

#### H. Cation–π rich neighborhood

- \(c_{\text{positive, aromatic}}\) in top 25% of dataset

**Annotation**: “Cation–π-rich segment; many local positive–aromatic contacts, often associated with strong, transient multivalent interactions.”

### Size / flexibility signatures

#### I. Small-rich flexible segment

- \(f_{\text{small}} \geq 0.50\)
- Small distribution: **dispersed**

**Annotation**: “Small-residue-rich segment; high local backbone flexibility and potential entropy-rich linker behavior.”

#### J. Bulky segment

- \(f_{\text{small}} \leq 0.20\)

**Annotation**: “Bulky side-chain segment; higher steric crowding and potential for packing-mediated interactions.”

### Neutral / bland signature

#### K. Chemically neutral linker

- All property fractions in moderate ranges (no extreme enrichment)
- All properties’ distributions **dispersed**
- Co-occurrence scores for key pairs (positive–aromatic, hydrophobic–hydrophobic) low

**Annotation**: “Chemically neutral linker; lacks strong compositional biases or local property clustering, likely generic flexible connector.”

## Implementation notes

- Thresholds can be refined, but the initial version uses fixed values as above.
- Distribution classes rely on simple run and gap statistics; no SCD/SHD metrics are used in v1.
- Co-occurrence uses a fixed window size (default 5) and property pairs only; no explicit AA motifs (e.g., RGG, SR) are considered in v1.
- Each LCR may carry 1–3 labels; the output includes:
  - `primary_physicochemical_annotation`
  - `secondary_physicochemical_annotations`
  - `annotation_evidence` (features/rules that triggered each label)

## Intended use

These annotations are designed to be used alongside:

- RNA-target superclass (e.g., mRNA vs rRNA binders)
- Domain-position class (intrinsic, edge, adjacent, linker, external)

to explore whether certain physicochemical behaviors are enriched in particular RBP classes or positional contexts.