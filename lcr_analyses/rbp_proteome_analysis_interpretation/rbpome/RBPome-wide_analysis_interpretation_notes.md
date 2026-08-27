# RBPome-wide analysis interpretation notes

## 1. LCR multiplicity per protein

### Method-specific notes

- **AlcoR:** Ribosomal proteins are enriched for exactly one LCR (32/52 = 61.5% vs 37.7%; log2OR +1.39, FDR 0.0098); unknown proteins are depleted for one LCR and enriched for 3+ LCRs.
- **CAST:** tRNA proteins are enriched for one LCR (84.3% vs 53.9%; log2OR +2.13, FDR 0.00013); unknown proteins are depleted for one LCR and enriched for 3+.
- **FLPS/LCRFinder:** mRNA is enriched for 3+ LCRs; ribosomal proteins and tRNA tend toward lower multiplicity.
- **SEG methods:** Sparse or detector-specific support only for several classes.

### Protein-level conclusion

mRNA-associated proteins tend toward multi-LCR architectures in permissive methods, whereas ribosomal-protein and tRNA-associated proteins tend toward lower multiplicity.

## 2. LCR coverage per protein

### What the metric measures

\[
\text{LCR coverage} = \frac{\text{residues assigned to one or more LCRs}}{\text{protein length}}
\]

Coverage differs from multiplicity: one long LCR can yield high coverage, while several short LCRs can yield low coverage.

### Method-specific notes

- **AlcoR:** Overall coverage differs across classes (FDR = 0.0190). Ribosomal proteins have higher coverage (median 0.446; Cliff’s δ = +0.228); ncRNA and snoRNA proteins have lower coverage.
- **CAST:** Overall coverage differs across classes (FDR = 0.00721); use Cliff’s δ for class-specific direction.
- **FLPS/LCRFinder/SEG:** Absolute values differ by detector and must be interpreted within method.

### Protein-level conclusion

Total LCR coverage differs among RNA-target classes; in AlcoR, ribosomal proteins are higher and ncRNA/snoRNA are lower.

## 3. Domain-position presence

### Core patterns across methods

Ribosomal proteins consistently show domain-edge and/or domain-intrinsic LCR enrichment with distal-terminal depletion. mRNA proteins show relatively more domain-adjacent or interdomain LCR presence in permissive methods.

### Method-specific notes

- **AlcoR:** Ribosomal proteins are domain-edge enriched (76.9% vs 54.5%; log2OR +1.44, FDR 0.016) and distal-terminal depleted (11.5% vs 37.2%; log2OR −2.09, FDR 0.0018).
- **CAST:** Strong ribosomal domain-edge/domain-intrinsic enrichment, with domain-adjacent and distal-terminal depletion.
- **FLPS:** Ribosomal domain-edge enrichment and distal-terminal depletion; mRNA domain-adjacent/interdomain-linker enrichment.
- **LCRFinder:** Ribosomal distal-terminal depletion and domain-intrinsic enrichment; mRNA domain-adjacent/interdomain-linker enrichment.
- **Strict SEG:** Too sparse for general inference.
- **Intermediate SEG:** Retains the ribosomal domain-edge/domain-intrinsic pattern.

### Protein-level conclusion

Ribosomal-protein-class LCRs are preferentially domain-edge and/or domain-intrinsic, while distal-terminal LCR presence is depleted.

### LCR-level descriptive conclusion

Detected ribosomal-protein LCRs concentrate around or inside domains. This remains LCR-level descriptive context; protein-level results carry the claim.

## 4. Physicochemical signature presence

### What the metric measures

A protein is positive when it contains at least one LCR assigned to a signature. This tests signature occurrence, not the amount of LCR sequence occupied by it.

### Core patterns across methods

- **mRNA:** RG/RGG-repeat and SR/RS-repeat enrichment; permissive methods also show GS-rich neutral and disorder-rich-spacer enrichment.
- **pre-rRNA/snoRNA:** Acidic and/or basic-enriched signatures.
- **Ribosomal protein:** Arginine-rich RNA-contact enrichment in informative methods; FLPS also finds aromatic-sticker/disorder-rich-spacer depletion.
- **tRNA:** Relative aromatic-sticker and disorder-rich-spacer depletion in FLPS.

### Method-specific notes

- **AlcoR:** mRNA is enriched for RG/RGG and SR/RS carriers.
- **CAST:** Pre-rRNA and ribosomal proteins are enriched for arginine-rich RNA-contact signatures; pre-rRNA, ribosomal protein, and snoRNA show basic-enriched signatures, while mRNA is depleted.
- **FLPS:** mRNA is enriched for SR/RS, RG/RGG, GS-rich neutral, and disorder-rich-spacer signatures; pre-rRNA/snoRNA show acidic/basic patterns; ribosomal proteins are arginine-rich but acidic/aromatic/disorder-spacer depleted.
- **LCRFinder:** mRNA has the same broad SR/RS, RG/RGG, GS-rich, and disorder-spacer pattern; pre-rRNA, snoRNA, and tRNA are depleted for disorder-rich spacer.

### Protein-level conclusion

The reproducible profiles are mRNA-associated RG/RGG and SR/RS repeats, pre-rRNA/snoRNA-associated acidic and basic signatures, and a ribosomal-protein-associated arginine-rich architecture.

## 5. Physicochemical property fraction per protein

### What the metric measures

\[
\text{property fraction}_{s} = \frac{\text{detected LCR residues assigned to signature }s}{\text{all detected LCR residues in the protein}}
\]

This is not whole-protein composition. It measures the relative prominence of a signature within a protein’s detected LCR residues; signature-specific n is the number of contributing proteins.

### Core patterns across methods

- **Basic-enriched fraction:** Strongly elevated in ribosomal proteins in AlcoR, FLPS, and LCRFinder; pre-rRNA and snoRNA also tend upward in FLPS/LCRFinder.
- **Arginine-rich fraction:** Strongly elevated in ribosomal proteins in AlcoR, FLPS, and LCRFinder.
- **Acidic fraction:** Higher in pre-rRNA in AlcoR, FLPS, LCRFinder, and intermediate SEG.
- **Disorder-rich-spacer fraction:** Varies across methods; no universal class ordering.
- **GS-rich neutral fraction:** Elevated in mRNA in AlcoR and FLPS.
- **Hydrophobic/aromatic fractions:** Some ribosomal-protein effects recur, but are more detector-dependent.
- **Other/unmapped fractions:** `Other` is elevated in mRNA and ribosomal proteins in AlcoR, FLPS, and LCRFinder; unmapped composition is not a biological signature.

### Method-specific notes

- **AlcoR:** Higher acidic fraction in pre-rRNA; higher arginine-rich/basic-enriched fractions in ribosomal proteins; higher GS-rich neutral fraction in mRNA.
- **CAST:** Supported global differences primarily for disorder-rich spacer, `other`, and unmapped properties.
- **FLPS:** Broadest supported profile: high basic-enriched/arginine-rich fractions in ribosomal proteins, high acidic fraction in pre-rRNA, and higher GS-rich neutral fraction in mRNA.
- **LCRFinder:** Higher acidic fraction in pre-rRNA; higher arginine-rich, basic-enriched, and hydrophobic fractions in ribosomal proteins; higher disorder-rich-spacer fraction in mRNA.
- **Strict SEG:** Too sparse for cross-method inference.
- **Intermediate SEG:** Supports acidic, basic-enriched, disorder-rich-spacer, `other`, and SR/RS differences, but several ribosomal-protein rows are descriptive-only.

### Protein-level conclusion

Ribosomal-protein LCR material is enriched for basic-enriched and arginine-rich RNA-contact features, while pre-rRNA LCR material tends toward a higher acidic-region fraction. mRNA LCR material recurrently has greater GS-rich-neutral and, in some methods, disorder-related composition.

## 6. LCR length

### What the metric measures

`continuous_lcr → length` is the amino-acid length of each detected LCR. This is an LCR-level descriptive analysis: proteins with multiple LCRs contribute multiple observations.

### Core patterns across methods

No length pattern is as method-consistent as the position or composition patterns. mRNA LCRs tend toward longer segments in FLPS/intermediate SEG; tRNA LCRs tend toward shorter segments in CAST, FLPS, and LCRFinder.

### Method-specific notes

- **AlcoR:** Overall FDR = 0.0360; individual effects are small.
- **CAST:** FDR \(=9.16\times10^{-6}\); ribosomal (median 42.5 aa; δ −0.339) and tRNA LCRs (52.5 aa; δ −0.292) are shorter.
- **FLPS:** FDR \(=3.10\times10^{-20}\); mRNA LCRs are longer (21 aa; δ +0.120), tRNA LCRs shorter (15 aa; δ −0.173).
- **LCRFinder:** FDR \(=2.79\times10^{-7}\); tRNA LCRs are shorter (5 aa; δ −0.148).
- **Strict SEG:** Not FDR-supported (0.112).
- **Intermediate SEG:** FDR \(=4.02\times10^{-7}\); mRNA is longer (30 aa; δ +0.162), pre-rRNA shorter (21 aa; δ −0.346).

### LCR-level descriptive conclusion

LCR length varies by class within several callers, but the class ordering is not sufficiently consistent for a strong method-independent conclusion.

## 7. LCR composition / charge

### What the metric measures

This LCR-level module uses individual-LCR composition metrics from `continuous_lcr`: `frac_polar`, `frac_hydrophobic`, `frac_strong_hydro`, `frac_aromatic`, `frac_disorder`, `frac_GS`, `frac_positive`, `frac_negative`, `fcr`, and `ncpr`.

`fcr` is the fraction of charged residues; `ncpr` is net charge per residue. Kruskal–Wallis tests whether a metric differs across classes within a method, while Cliff’s δ gives the direction of each class against all others. This is descriptive because proteins carrying multiple LCRs contribute repeatedly.

### Core patterns across methods

- **Ribosomal proteins:** The strongest recurring composition profile. Their LCRs tend to have higher positive-residue fraction, higher net charge per residue, lower negative-residue fraction, and often lower GS fraction. High `ncpr` recurs in AlcoR (+0.417), CAST (+0.542), FLPS (+0.402), LCRFinder (+0.324), and intermediate SEG (+0.428).
- **pre-rRNA:** LCRs tend to have higher total charged fraction (`fcr`) and polar/charged composition in several methods, but are not consistently net-positive.
- **mRNA:** LCRs tend toward GS-rich and disorder-associated composition rather than strong charge specialization. GS is positive in AlcoR, CAST, and FLPS; disorder is positive in FLPS and LCRFinder.

### Method-specific notes

- **AlcoR:** Ribosomal-protein LCRs have higher hydrophobic fraction (δ +0.283), positive fraction (+0.276), and `ncpr` (+0.417), with lower negative fraction (−0.311) and GS fraction (−0.182).
- **CAST:** Ribosomal LCRs have higher hydrophobic fraction (+0.363), positive fraction (+0.541), `fcr` (+0.335), and `ncpr` (+0.542), with lower GS fraction (−0.467). Pre-rRNA has high `fcr` (+0.462), consistent with charged but not uniformly cationic LCRs.
- **FLPS:** Ribosomal LCRs have higher hydrophobic (+0.254), positive (+0.373), `fcr` (+0.182), and `ncpr` (+0.402), with lower negative (−0.168) and GS (−0.159) fractions. Pre-rRNA has increased polar, positive, negative, and `fcr` effects; mRNA has lower hydrophobic and mildly higher disorder composition.
- **LCRFinder:** Ribosomal LCRs have strong hydrophobic enrichment (+0.386), positive `ncpr` (+0.324), positive disorder (+0.235), and lower GS (−0.259). Pre-rRNA has high `fcr` (+0.362) but negative `ncpr`.
- **Strict SEG:** Sparse ribosomal-protein and tRNA LCR counts make large apparent effects descriptive-only.
- **Intermediate SEG:** Reinforces the ribosomal positive-charge profile: positive fraction +0.354, `ncpr` +0.428, and lower GS fraction −0.281.

### Other classes

- **snRNA:** LCRs often show higher positive-residue fraction or more charged/polar character, but their direction is less consistent across callers than the ribosomal-protein cationic profile.
- **snoRNA:** LCRs tend toward more polar and/or acidic composition in several methods. The stronger evidence for snoRNA remains the protein-level acidic/basic signature-presence and property-fraction modules.
- **tRNA:** LCR composition is method-dependent. FLPS and LCRFinder identify relatively greater hydrophobic contributions, whereas strict SEG and intermediate SEG have limited class support for broad tRNA composition claims.
- **Unknown:** LCRs in the unknown class lack a stable cross-method composition profile and should not be interpreted as representing a single biological RNA-target class.

### LCR-level descriptive conclusion

Individual ribosomal-protein LCRs show the most reproducible compositional specialization: relatively greater positive-residue content and net charge, less negative charge, and often lower GS content. Pre-rRNA LCRs are comparatively charged and polar but not consistently net-positive, whereas mRNA LCRs tend toward GS- and disorder-associated composition. These results describe detected LCR populations and complement, rather than replace, the protein-level findings.

## 8. Within-class LCR position/signature proportions

### What the module measures

This final module describes the makeup of the detected LCR pool within each RNA-target class, using `categorical_lcr`. It includes positional category proportions, the fraction of individual LCRs assigned to each signature (`signature_sequence_count_lcr`), and signature representation on a residue basis (`signature_residue_count_lcr`).

The module is descriptive-only. Fisher tests and log2 odds ratios can contrast a class with all other classes, but proteins with multiple LCRs contribute multiple observations; protein-level analyses remain the basis for biological claims.

### Position makeup

The LCR-level proportions reinforce the established ribosomal-protein architecture:

- **FLPS:** 37.3% of ribosomal-protein LCRs are domain-edge versus 10.1% in other classes; 28.0% are domain-intrinsic versus 9.5%; and 7.3% are distal-terminal versus 52.6%.
- **CAST:** 68.1% of ribosomal-protein LCRs are domain-edge versus 15.5% elsewhere; 2.1% are distal-terminal versus 44.4%.
- **AlcoR:** 52.6% are domain-edge versus 33.2% elsewhere, while 6.2% are distal-terminal versus 36.5%.
- **LCRFinder:** Ribosomal-protein LCRs are enriched at domain edges and within domains, while distal-terminal LCRs are depleted.

mRNA LCR populations contain relatively more domain-adjacent and interdomain-linker segments in FLPS and LCRFinder, consistent with the protein-level position-presence results.

### Signature makeup

- **mRNA LCR population:** Repeat- and low-complexity-associated signatures are overrepresented. In FLPS, 3.37% of mRNA LCRs are RG/RGG-repeat versus 1.22% elsewhere, 3.82% are SR/RS-repeat versus 0.87%, and 5.31% are GS-rich neutral versus 2.72%. AlcoR and LCRFinder similarly support enrichment for RG/RGG and/or SR/RS LCRs.

- **pre-rRNA and snoRNA LCR populations:** Acidic and basic features are overrepresented. In FLPS, pre-rRNA LCRs are 29.6% acidic versus 16.1% elsewhere and 19.5% basic-enriched versus 8.6%; snoRNA LCRs are 31.4% acidic versus 16.6% and 27.0% basic-enriched versus 8.8%. LCRFinder and intermediate SEG provide convergent acidic/basic patterns.

- **Ribosomal-protein LCR population:** RNA-contact/basic signatures are enriched, whereas GS-rich-neutral and disorder-rich-spacer LCRs are depleted in FLPS. In FLPS, 34.3% of ribosomal-protein LCRs are arginine-rich RNA-contact segments versus 9.6% elsewhere and 19.1% are basic-enriched versus 9.0%. CAST also shows 27.1% arginine-rich versus 3.7% and 33.3% basic-enriched versus 6.2%.

### LCR-level descriptive conclusion

Within-class LCR proportions reinforce the protein-level findings: ribosomal proteins contribute predominantly domain-edge/domain-intrinsic, basic and arginine-rich LCRs; pre-rRNA and snoRNA contribute acidic/basic LCR populations; and mRNA proteins contribute more SR/RS-, RG/RGG-, and GS-associated LCRs.

This is a descriptive consistency check, not a separate source of protein-level inference.
