# =============================================================================
# Script:      15_gsva_run.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Hernandez
# Description: Runs Gene Set Variation Analysis (GSVA) using the GO Biological
#              Process gene sets built in script 16. GSVA estimates per-sample
#              gene set enrichment scores in a non-parametric, unsupervised
#              manner, producing a matrix of enrichment scores (rows = GO
#              terms, columns = samples) that is then passed to limma for
#              differential enrichment testing (script 18).
#
#              Parameters:
#                kcdf = "Gaussian"  : appropriate for log2-normalized microarray
#                minSize = 10       : minimum gene set size (same as filter)
#                maxSize = 300      : maximum gene set size
# Inputs:      data/processed/GSVA_analysis/exp_1to1.rds
#              data/processed/GSVA_analysis/genesbygo_filtered.rds
# Outputs:     data/processed/GSVA_analysis/gsva_results.rds  (ExpressionSet)
#              data/processed/GSVA_analysis/gsva_r.rds         (plain data.frame)
# Dependencies: 00_setup.R, 15_gsva_ortholog_mapping.R, 16_gsva_go_genesets.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(GSVA)
library(Biobase)

exp_1to1           <- readRDS(file.path(PATH_GSVA_PROC, "exp_1to1.rds"))
genesbygo_filtered <- readRDS(file.path(PATH_GSVA_PROC, "genesbygo_filtered.rds"))

# ---- 1. Wrap expression matrix in an ExpressionSet ----
eset <- ExpressionSet(assayData = as.matrix(exp_1to1))

# ---- 2. Define GSVA parameters ----
param <- gsvaParam(
  expr     = eset,
  geneSets = genesbygo_filtered,
  kcdf     = "Gaussian",   # suitable for continuous (log2 microarray) data
  minSize  = 10,
  maxSize  = 300
)

# ---- 3. Run GSVA ----
cat("Running GSVA over", length(genesbygo_filtered), "GO BP gene sets",
    "and", ncol(exp_1to1), "samples...\n")
set.seed(123)
gsva_results <- gsva(param, verbose = TRUE)

cat("\nGSVA complete. Score matrix dimensions:", dim(gsva_results), "\n")

# ---- 4. Extract plain data.frame for downstream limma ----
gsva_r <- as.data.frame(exprs(gsva_results))
cat("Score matrix: rows (GO terms) =", nrow(gsva_r),
    " | columns (samples) =", ncol(gsva_r), "\n")

# ---- 5. Save ----
saveRDS(gsva_results, file.path(PATH_GSVA_PROC, "gsva_results.rds"))
saveRDS(gsva_r,       file.path(PATH_GSVA_PROC, "gsva_r.rds"))
cat("Saved: gsva_results.rds (ExpressionSet) | gsva_r.rds (data.frame)\n")
