Brief introduction

Epistatic effects refer to non-additive interactions between loci controlling the same trait. Genomic selection (GS) accelerates plant breeding by predicting phenotypic effects in individuals based on genomic marker profiles. However, current mainstream GS models primarily construct prediction equations based on additive genetic effects, failing to fully leverage non-additive effects such as epistasis to enhance prediction accuracy. This study identified significant dominant interaction SNP pairs through pathway-level dominance testing, incorporated them into the GS model, and evaluated the enhancement of genomic prediction accuracy by epistatic effects.

BL_epistasis.R

The Bayesian BL model is used as an example to illustrate how to integrate epistatic interactions into the GS model.

RUN

“run_bridge_epistasis_pipeline.sh” is the code used in this institute for detecting epistasis effects and analyzing pathway interactions.
The script performs the following key stages:
1. GWAS data preprocessing (VCF → PLINK format + QC + LD pruning)
2. BridGE data processing (SNP-to-gene and gene-to-pathway mapping)
3. Pairwise epistasis detection using CASSI on real data
4. Generation of random phenotype permutations
5. CASSI analysis on randomized datasets
6. Conversion of CASSI results to network pickle format
7. Computation of BPM/WPM/PATH statistics (real + permutations)
8. FDR calculation
9. Summarization and export of significant epistatic interactions
