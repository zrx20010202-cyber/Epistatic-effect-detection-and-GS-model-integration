#!/usr/bin/env bash
# =============================================================================
# BridGE Pipeline Script: Pathway-based Genetic Interaction Analysis
# This script organizes and executes the core steps of the BridGE analysis workflow using BridGE-Python and CASSI, following the methodology described in the protocol paper: https://doi.org/10.1038/s41596-024-00954-8
# Required reference files and related executable scripts should be obtained from: https://github.com/csbio/BridGE-Python
# All original file names and command parameters are preserved for consistency with the described methods.
# =============================================================================

set -euo pipefail

# ── Configuration Section ──────────────────────────────────────────────────────
# Adjust these variables according to your project setup

PROJECT_DIR="data"
INPUT_VCF="input/final_filtered.vcf.gz"
PREFIX_RAW="gwas_data"
PREFIX_FINAL="gwas_data_final"

CASSI_REAL="ssM_cassi_LR_R0.txt"
PICKLE_REAL="ssM_cassi_LR_R0.pkl"

N_RANDOM=10
SNP_PERMS=10000
SAMPLE_PERMS=10
FDR_CUT=0.25
DENSITY_CUT=0.1

GENE_ANNO="gene_annotation.txt"
GENE_SET="gene_set"
SNP_PATH_PKL="snp_pathway_min10_max300.pkl"

N_WORKER=10

# ── Execution ──────────────────────────────────────────────────────────────────

echo "=== BridGE Pathway Genetic Interaction Analysis Pipeline ==="
echo "Project directory: ${PROJECT_DIR}"
echo "Input VCF:         ${INPUT_VCF}"
echo ""

# Stage 1 - Preprocessing GWAS data
echo "[1] Preprocessing: VCF → PLINK binary format + QC + LD pruning"
plink --vcf ${INPUT_VCF} --double-id --make-bed --out ${PREFIX_RAW}

bash data_checkpopulation.sh --plinkFile=${PREFIX_RAW}
bash preprocessgwas.sh --plinkFile=${PREFIX_RAW} --ldR2=0.02
bash LD_Dprime_filter.sh --plinkFile=${PROJECT_DIR}/intermediate/${PREFIX_FINAL}

[ ! -f "${PREFIX_FINAL}.bed" ] && { echo "Error: ${PREFIX_FINAL}.bed was not generated"; exit 1; }

# Stage 2 - BridGE data processing (SNP to gene/pathway mapping)
echo "[2] BridGE data processing"
python bridge.py --projectDir ${PROJECT_DIR} --job=DataProcess --plinkFile=${PREFIX_FINAL} --geneAnnotation=${GENE_ANNO} --genesets=${GENE_SET}

# Stage 3 - CASSI pairwise interaction detection on real data
echo "[3] CASSI pairwise interaction detection (real data)"
cassi -i ${PREFIX_FINAL}.bed -lin -lin-th 0.1 -max 0 -o ${CASSI_REAL}

# Stage 4 - Generate random phenotypes and corresponding PLINK files
echo "[4] Generating ${N_RANDOM} sets of random phenotypes"
plink --bfile ${PREFIX_FINAL} --make-perm-pheno ${N_RANDOM} --allow-no-sex --out rand_pheno

for j in $(seq 1 ${N_RANDOM}); do
    plink --bfile ${PREFIX_FINAL} --pheno rand_pheno.pphe --mpheno ${j} --make-bed --out ${PREFIX_FINAL}_R${j}
done

# Stage 5 - CASSI on randomized datasets
echo "[5] CASSI pairwise interaction detection on randomized sets"
for j in $(seq 1 ${N_RANDOM}); do
    echo "  Processing random set R${j}"
    cassi -i ${PREFIX_FINAL}_R${j}.bed -lin -lin-th 0.1 -max 0 -o ssM_cassi_LR_R${j}.txt
done

# Stage 6 - Convert CASSI outputs to pickle format
echo "[6] Converting CASSI output to pickle network format"
python -m corefuns.cassissm ${PREFIX_FINAL}.bim ${CASSI_REAL} lin ${PICKLE_REAL}

for j in $(seq 1 ${N_RANDOM}); do
    python -m corefuns.cassissm ${PREFIX_FINAL}.bim ssM_cassi_LR_R${j}.txt lin ssM_cassi_LR_R${j}.pkl
done

# Stage 7 - Compute BPM/WPM/PATH statistics (real + random)
echo "[7] Computing BPM/WPM/PATH statistics"
for j in $(seq 0 ${N_RANDOM}); do
    echo "  Processing R${j}"
    python bridge.py --projectDir ${PROJECT_DIR} --job=ComputeStats --ssmfile=ssM_cassi_LR_R${j}.pkl --nWorker=${N_WORKER} --snpPerms=${SNP_PERMS} --minPath=10
done

# Stage 8 - Compute FDR
echo "[8] Computing FDR"
python bridge.py --projectDir ${PROJECT_DIR} --job=ComputeFDR --ssmfile=${PICKLE_REAL} --pvalueCutoff=0.05 --minPath=10 --samplePerms=${SAMPLE_PERMS}

# Stage 9 - Summarize and export results
echo "[9] Summarizing and exporting significant results"
python bridge.py --projectDir ${PROJECT_DIR} --job=Summarize --ssmfile=${PICKLE_REAL} --fdrcut=${FDR_CUT} --snpPathFile=${SNP_PATH_PKL} --densityCutoff=${DENSITY_CUT}

echo ""
echo "=== Analysis pipeline completed ==="
echo "Check output files in the ${PROJECT_DIR} directory, especially FDR-corrected significant interaction lists"
echo "If you encounter any issues, please refer to the original code and detailed protocol in https://doi.org/10.1038/s41596-024-00954-8 for adjustments"