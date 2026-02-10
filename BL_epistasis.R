#!/usr/bin/env Rscript

# Bayesian LASSO with epistasis interactions using BGLR
# 50 repeats × 5-fold cross-validation

# Required packages
library(vcfR)      
library(BGLR)
library(data.table)

# ── Please replace these with your actual local file paths ──
# These are example placeholder names only
vcf_file    <- "your_own_filtered.vcf"              # Your own pre-filtered VCF file
pheno_file  <- "your_own_phenotype.txt"             # Your own centered & scaled phenotype file
wpm_file    <- "WPM.csv"                            # Your WPM interaction list
bpm_file    <- "BPM.csv"                            # Your BPM interaction list

# ── Main analysis ───────────────────────────────────────────────────────────────

cat("=== Bayesian LASSO + Epistasis (BGLR) ===\n")
cat("Loading data...\n")

# 1. Load pre-processed, centered & scaled genotype matrix (0/1/2 coding)
# Expected format: first column = sample IDs, other columns = SNPs
Markers <- as.matrix(fread("genotype_matrix_scaled.csv", data.table = FALSE))
rownames(Markers) <- Markers[, 1]
Markers <- Markers[, -1, drop = FALSE]

# Use generic SNP names to avoid exposing original IDs
colnames(Markers) <- paste0("m", sprintf("%05d", seq_len(ncol(Markers))))

# 2. Load pre-scaled phenotype
pheno <- fread(pheno_file, header = TRUE, data.table = FALSE)
y <- as.numeric(pheno[[2]])         
names(y) <- pheno[[1]]

# 3. Keep only samples present in both genotype and phenotype
common <- intersect(rownames(Markers), names(y))
Markers <- Markers[common, , drop = FALSE]
y <- y[common]

cat(sprintf("Number of samples used: %d\n", length(y)))

# 4. Load and filter epistatic interaction pairs
interactions <- rbind(
  read.csv(wpm_file, stringsAsFactors = FALSE),
  read.csv(bpm_file, stringsAsFactors = FALSE)
)

interactions <- interactions[interactions$effect == "protective" & interactions$GI > 1, ]

snp1 <- as.character(interactions$snp1)
snp2 <- as.character(interactions$snp2)

# Only keep pairs where both SNPs exist in the genotype matrix
keep <- snp1 %in% colnames(Markers) & snp2 %in% colnames(Markers)
snp1 <- snp1[keep]
snp2 <- snp2[keep]

cat(sprintf("Number of epistatic pairs used for modeling: %d\n", length(snp1)))

# 5. Construct epistasis (interaction) matrix
EpiX <- matrix(NA_real_, nrow = nrow(Markers), ncol = length(snp1))

for (i in seq_along(snp1)) {
  EpiX[, i] <- Markers[, snp1[i]] * Markers[, snp2[i]]
}

# Center and scale the interaction terms
EpiX <- scale(EpiX)

# ── Cross-validation ────────────────────────────────────────────────────────────

set.seed(123)   

n_repeats <- 50
n_folds   <- 5
accuracies <- numeric(n_repeats * n_folds)
counter <- 0

cat("Starting 50 × 5 cross-validation...\n\n")

for (rep in 1:n_repeats) {
  folds <- sample(rep(1:n_folds, length.out = length(y)))
  
  for (k in 1:n_folds) {
    counter <- counter + 1
    
    trn <- folds != k
    tst <- !trn
    
    y_trn <- y[trn]
    y_tst <- y[tst]
    
    ETA <- list(
      list(X = EpiX[trn, , drop = FALSE], model = "BRR"),   
      list(X = Markers[trn, , drop = FALSE], model = "BL")  
    )
    
    fm <- BGLR(
      y       = y_trn,
      ETA     = ETA,
      nIter   = 6000,
      burnIn  = 2000,
      thin    = 5,
      verbose = FALSE
    )
    
    pred <- fm$mu +
      EpiX[tst, , drop = FALSE] %*% fm$ETA[[1]]$b +
      Markers[tst, , drop = FALSE] %*% fm$ETA[[2]]$b
    
    r <- cor(pred, y_tst, use = "complete.obs")
    accuracies[counter] <- r
    
    cat(sprintf("Rep %2d | Fold %d   r = %.4f\n", rep, k, r))
  }
}

# ── Final summary ───────────────────────────────────────────────────────────────

mean_r <- mean(accuracies, na.rm = TRUE)
sd_r   <- sd(accuracies, na.rm = TRUE)

cat("\n=== Final result ===\n")
cat(sprintf("Mean prediction correlation: %.4f ± %.4f\n", mean_r, sd_r))
cat(sprintf("Number of CV folds used: %d\n", length(accuracies)))
cat("\nAnalysis completed.\n")