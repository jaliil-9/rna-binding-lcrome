# Full-proteome background interpretation and validation

## 1. RBPome-wide LCR incidence

### Core patterns across methods

CAST, LCRFinder, strict SEG, and intermediate SEG support higher LCR incidence in RBPs than in the non-RBP proteome. AlcoR is near-null; FLPS is near-saturated.

### Method-specific notes

- **AlcoR:** 43.7% RBPs versus 43.4% background; log2OR +0.015, FDR 0.857.
- **CAST:** 55.7% versus 43.4%; log2OR +0.716, FDR \(3.76\times10^{-22}\).
- **FLPS:** 89.3% versus 87.5%; log2OR +0.247, FDR 0.0360; near-saturated.
- **LCRFinder:** 53.1% versus 43.7%; log2OR +0.546, FDR \(1.45\times10^{-13}\).
- **Strict SEG:** 15.8% versus 8.53%; log2OR +1.013, FDR \(5.72\times10^{-20}\).
- **Intermediate SEG:** 32.6% versus 22.4%; log2OR +0.749, FDR \(5.47\times10^{-20}\).

### Conclusion

Across informative non-saturated detectors, RBPs are more likely than non-RBP proteins to contain detected LCRs.

## 2. RBPome-wide incidence after length control

### Core patterns across methods

CAST, LCRFinder, strict SEG, and intermediate SEG retain RBP enrichment through multiple length strata. The RBPome-wide incidence difference is not explained by longer average RBP length.

### Method-specific notes

- **AlcoR:** No length bin is FDR-supported.
- **CAST:** Enrichment in bins 1–4: log2OR +0.506, +1.029, +0.416, and +0.762.
- **FLPS:** No bin survives FDR correction because of saturation.
- **LCRFinder:** Enrichment in the three central bins: +0.702, +0.477, and +0.482.
- **Strict SEG:** Enrichment in bins 2–5: +1.599, +0.863, +0.764, and +0.890.
- **Intermediate SEG:** Enrichment in bins 2–5: +0.923, +0.629, +0.916, and +0.482.

### Conclusion

The RBPome-wide excess of LCR carriers remains after length stratification in the informative methods.

## 3. RNA-class LCR incidence

### Core patterns across methods

mRNA, pre-rRNA, snRNA, snoRNA, and often unknown proteins are externally LCR-enriched; ribosomal-protein and tRNA classes are externally LCR-depleted. ncRNA is positive only in some methods.

### Method-specific notes

- **AlcoR:** Only ribosomal proteins are depleted: 31.2% versus 43.4%; log2OR −0.751, FDR 0.0185.
- **CAST:** mRNA, ncRNA, pre-rRNA, snRNA, snoRNA, and unknown proteins are enriched; ribosomal protein and tRNA are depleted.
- **FLPS:** mRNA is modestly enriched and ribosomal proteins depleted; saturation limits other contrasts.
- **LCRFinder:** mRNA, ncRNA, pre-rRNA, snRNA, snoRNA, and unknown are enriched; ribosomal protein and tRNA are depleted.
- **Strict SEG:** mRNA, snRNA, snoRNA, and unknown are enriched; ribosomal protein and tRNA are depleted.
- **Intermediate SEG:** mRNA, pre-rRNA, snRNA, snoRNA, and unknown are enriched; ribosomal protein and tRNA are depleted.

### External validation implications

> **mRNA multi-LCR/repeat-spacer architecture is strengthened by external enrichment. Ribosomal-protein and tRNA low-LCR patterns are absolute proteome-level depletion, not merely relative RBPome differences. Ribosomal-protein domain/basic/arginine-rich architecture belongs to a smaller carrier subset.**

### Conclusion

External incidence separates an LCR-enriched group of RNA classes from LCR-depleted ribosomal-protein and tRNA classes.

## 4. RBP LCR architecture among carriers

### Core patterns across methods

Among LCR carriers, RBPs generally have more LCRs per protein. Coverage and mean LCR length are also higher in most methods, except AlcoR.

### Method-specific notes

- **AlcoR:** No architecture difference.
- **CAST:** Higher coverage (δ +0.056), LCR number (+0.081), and mean length (+0.042).
- **FLPS:** Higher coverage (+0.052), LCR number (+0.077), and mean length (+0.085).
- **LCRFinder:** Higher coverage (+0.113), LCR number (+0.124), and mean length (+0.133).
- **Strict SEG:** More LCRs (+0.141) and longer mean LCRs (+0.146); coverage is not supported.
- **Intermediate SEG:** Higher coverage (+0.146), LCR number (+0.188), and mean length (+0.178).

### Conclusion

Conditional on carrying an LCR, RBP proteins generally have more elaborate LCR architecture than non-RBP carriers.

## 5. RBP LCR composition and charge

### Core patterns across methods

RBP LCRs are more polar, disordered, and charge-rich in all six methods. Positive and negative fractions, total charged fraction, and net charge are higher; hydrophobic and strongly hydrophobic fractions are lower.

### Method-specific notes

- **AlcoR:** Polar +0.107, disorder +0.097, positive +0.129, negative +0.072, `fcr` +0.148, and `ncpr` +0.049; hydrophobic properties are lower.
- **CAST:** Polar +0.169, disorder +0.092, positive +0.191, negative +0.096, `fcr` +0.192, and `ncpr` +0.086.
- **FLPS:** Polar +0.154, disorder +0.150, positive +0.122, negative +0.091, `fcr` +0.173, and `ncpr` +0.050.
- **LCRFinder:** Polar +0.152, disorder +0.132, positive +0.162, `fcr` +0.187, and `ncpr` +0.080.
- **Strict SEG:** Polar +0.246, disorder +0.097, positive +0.234, negative +0.099, `fcr` +0.239, and `ncpr` +0.145.
- **Intermediate SEG:** Polar +0.280, positive +0.197, negative +0.128, `fcr` +0.238, and `ncpr` +0.104.

### External validation implications

> **The internal composition/charge findings are validated and generalized.** Internal RNA-class composition patterns occur on top of an RBPome-wide polar, disordered, charge-rich LCR background; later class tests must separate generic from class-specific effects.

### Conclusion

RBP-associated LCRs are polar, disordered, charge-rich, modestly cationic, and depleted in hydrophobic sequence character.

## 6. RBPome-wide signature carrier rates

### Core patterns across methods

RG/RGG-repeat, basic-enriched, arginine-rich RNA-contact, and acidic signatures are enriched across the RBPome. SR/RS and GS-rich-neutral signatures are enriched in most methods; hydrophobic signatures are depleted. Aromatic-sticker and disorder-rich-spacer results are method-dependent.

### Method-specific notes

All six methods support enrichment of acidic, basic-enriched, arginine-rich, and RG/RGG signatures, with strong SR/RS enrichment in five methods. Hydrophobic carriers are depleted across the informative methods.

### External validation implications

> **Repeat-, acidic-, basic-, and arginine-rich signatures are broad RBPome features.** Their internal enrichment in a class is not automatically class-exclusive biology.

### Conclusion

The RBPome is broadly enriched for charge/RNA-contact/repeat signatures and depleted for hydrophobic signature carriers.

## 7. Signature occupancy among carriers

### Core patterns across methods

Carrier enrichment is stronger and more reproducible than occupancy enrichment. Most RBPome-wide signature effects reflect more RBP proteins carrying a signature, not signature expansion across every carrier.

### Method-specific notes

- **AlcoR:** Hydrophobic, aromatic-sticker, and unmapped occupancy are lower.
- **FLPS:** Acidic and basic occupancy are higher; aromatic-sticker, cation–π, and hydrophobic occupancy are lower.
- **Intermediate SEG:** Acidic, basic, and SR/RS occupancy are higher; hydrophobic occupancy is lower.
- **CAST/LCRFinder:** Occupancy effects are weaker than carrier-rate effects.
- **Strict SEG:** Sparse carrier subsets limit inference.

### External validation implications

> **The RBPome-wide signature landscape is mainly a prevalence effect.** Class-specific interpretation must distinguish increased carrier frequency from increased occupancy among carriers.

### Conclusion

Charged, RNA-contact, and repeat signatures are primarily enriched through increased carrier frequency; hydrophobic occupancy is reduced among RBP carriers.

## 8. External validation of class signatures

### Core patterns across methods

The strongest externally validated signatures are mRNA RG/RGG and SR/RS repeats, pre-rRNA acidic/basic signatures, and snoRNA acidic/basic signatures. Many other external enrichments are new external signals rather than internally class-specific patterns.

### Method-specific notes

- **mRNA:** RG/RGG is externally validated in 4/6 methods, SR/RS in 5/6, GS-rich neutral in 3/6, and disorder-rich spacer in 2/6.
- **pre-rRNA:** Acidic enrichment is validated in 4/6 methods and basic enrichment in 3/6.
- **snoRNA:** Acidic and basic patterns are supported but less decisive because several methods are underpowered.
- **Ribosomal protein:** Arginine-rich RNA-contact enrichment is validated in 2 methods; many other signature rows are null or underpowered.
- **snRNA:** Basic and SR/RS signatures are mainly new external signals.
- **tRNA:** Mostly null or underpowered; some signature depletions validate.
- **Unknown:** Acidic signature is a new external signal in all six methods, consistent with class heterogeneity.

### External validation implications

> **Retain:** mRNA RG/RGG and SR/RS repeats; pre-rRNA acidic/basic enrichment; tentative snoRNA acidic/basic enrichment.  
> **Reframe:** many charged, arginine-rich, acidic, and disorder-related class signals as RBP-generic or external-only.  
> **Downgrade:** broad ribosomal-protein signature claims because of recurrent external underpowering.

### Conclusion

The external comparison retains the mRNA repeat-rich and pre-rRNA acidic/basic profiles, while showing that several internal signals are broader RBPome properties rather than exclusive class-specific biology.
