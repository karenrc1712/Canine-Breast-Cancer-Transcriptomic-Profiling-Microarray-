# =============================================================================
# Script:      11_gsva_prepare_expression.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Fernández
# Description: Loads the annotated RMA expression table produced in script 05
#              and prepares a gene-level matrix suitable for GSVA:
#                - Separates metadata columns from sample intensity columns
#                - Resolves duplicate gene symbols by averaging probe
#                  intensities across all probesets mapping to the same gene
#              The resulting matrix has one row per unique gene symbol and
#              one column per sample.
# Inputs:      data/processed/probes_named_rma.rds
# Outputs:     data/processed/GSVA_analysis/expresion_final.rds
# Dependencies: 00_setup.R, 05_gene_annotation.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(dplyr)

probes_named_rma <- readRDS(file.path(PATH_PROCESSED, "probes_named_rma.rds"))

# ---- 1. Separate metadata columns from sample columns ----
metadata_cols <- c("transcript_cluster_id", "gene_assignment", "n_probesets", "gene_symbol")
sample_cols   <- setdiff(colnames(probes_named_rma), metadata_cols)

cat("Metadata columns (", length(metadata_cols), "):",
    paste(metadata_cols, collapse = ", "), "\n")
cat("Sample columns   (", length(sample_cols),   "):", length(sample_cols), "samples\n")

# ---- 2. Build raw expression matrix ----
expresion_matrix <- as.matrix(probes_named_rma[, sample_cols])
rownames(expresion_matrix) <- probes_named_rma$transcript_cluster_id

# ---- 3. Collapse duplicate gene symbols by averaging probeset intensities ----
# Multiple probesets can map to the same gene. GSVA requires one row per gene,
# so we take the mean across all probesets for each duplicate symbol.
dup_genes <- names(which(table(probes_named_rma$gene_symbol) > 1))
cat("Gene symbols with >1 probeset:", length(dup_genes), "\n")

exp_df <- as.data.frame(expresion_matrix)
exp_df$gene_symbol <- probes_named_rma$gene_symbol

expresion_final <- exp_df %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  group_by(gene_symbol) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

expresion_final_mat <- as.matrix(expresion_final[, -1])
rownames(expresion_final_mat) <- expresion_final$gene_symbol

cat("Final expression matrix dimensions:", dim(expresion_final_mat), "\n")
cat("  Rows (unique gene symbols):", nrow(expresion_final_mat), "\n")
cat("  Columns (samples):         ", ncol(expresion_final_mat), "\n")

# ---- 4. Save ----
saveRDS(expresion_final_mat,
        file.path(PATH_GSVA_PROC, "expresion_final.rds"))

