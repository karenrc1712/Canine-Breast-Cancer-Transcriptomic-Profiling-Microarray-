# =============================================================================
# Script:      12_gsva_sva_batch_estimation.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Fernández
# Description: Estimates the number and identity of surrogate variables (SVs)
#              in the expression matrix using the SVA package (Leek method).
#              SVs capture unwanted technical variation (batch effects, sample
#              handling differences) not explained by the biological groups.
#
#              NOTE: This script is DIAGNOSTIC. Estimated SVs are inspected
#              but are not automatically included in downstream models. The
#              researcher should evaluate svobj$sv and decide whether to
#              incorporate SVs into the GSVA design matrix (script 18).
# Inputs:      data/processed/GSVA_analysis/expresion_final.rds
#              data/processed/targets1.rds
# Outputs:     data/processed/GSVA_analysis/svobj.rds
# Dependencies: 00_setup.R, 13_gsva_prepare_expression.R, 03_load_targets.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(sva)

expresion_final <- readRDS(file.path(PATH_GSVA_PROC, "expresion_final.rds"))
targets1        <- readRDS(file.path(PATH_PROCESSED,  "targets1.rds"))

# ---- 1. Align column names with targets ----
# Column names in expresion_final are numeric strings; match them to Renomear
colnames(expresion_final) <- as.numeric(colnames(expresion_final))

targets_aligned <- targets1[match(colnames(expresion_final), targets1$Renomear), ]
stopifnot("Column-target mismatch" =
            all(colnames(expresion_final) == targets_aligned$Renomear))

# ---- 2. Build full and null models ----
group <- factor(targets_aligned$Grupos)
mod   <- model.matrix(~ group)   # full model: includes biological group
mod0  <- model.matrix(~ 1)       # null model: intercept only

# ---- 3. Estimate number of surrogate variables ----
n.sv <- num.sv(expresion_final, mod, method = "leek")
cat("Estimated number of surrogate variables:", n.sv, "\n")

# ---- 4. Compute surrogate variables ----
svobj <- sva(expresion_final, mod, mod0, n.sv = n.sv)

cat("SVA complete. Surrogate variable matrix dimensions:",
    dim(svobj$sv), "\n")
cat("First few SV values:\n")
print(head(svobj$sv))

# ---- 5. Save ----
saveRDS(svobj, file.path(PATH_GSVA_PROC, "svobj.rds"))
cat("Saved: data/processed/GSVA_analysis/svobj.rds\n")
cat("\nNext step: inspect svobj$sv and decide whether to include SVs\n")
cat("in the limma design matrix inside 18_gsva_differential_enrichment.R\n")

