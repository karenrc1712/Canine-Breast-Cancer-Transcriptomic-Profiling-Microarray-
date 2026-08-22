# =============================================================================
#
# Script:      04_pca_exploratory.R
#
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
#
# Author:      Karen Rodriguez
#
# Advisor:     Geysson Javier Fernández
#
# Description: Exploratory PCA using the 5000 most variable genes from the
#              entire normalized matrix (not yet filtered by significance).
#
#              This serves as a quality control step: if the groups do not
#              separate at all here, batch effects should be investigated
#              before proceeding.
#
# Inputs:      data/processed/exprs_rma.rds
#              data/processed/targets1.rds
#
# Outputs:     Plot displayed on screen (not exported to a file; this step is
#              performed later in 08_pca_degs.R using the already filtered genes)
#
# Dependencies: 00_setup.R, 02_normalization_rma.R, 03_load_targets.R
#
# =============================================================================

source(here::here("R", "00_setup.R"))
library(matrixStats)
library(ggplot2)

exprs.rma <- readRDS(file.path(PATH_PROCESSED, "exprs_rma.rds"))
targets1  <- readRDS(file.path(PATH_PROCESSED, "targets1.rds"))

# ---- 1. Filter the 5000 genes with highest variance ----
vars           <- rowVars(exprs.rma)
top_var_genes  <- order(vars, decreasing = TRUE)[1:5000]
exprs_top      <- exprs.rma[top_var_genes, ]

# ---- 2. PCA ----
pca_post_rma <- prcomp(t(exprs_top), scale. = TRUE)

pca_data <- data.frame(
  Sample = colnames(exprs_top),
  PC1    = pca_post_rma$x[, 1],
  PC2    = pca_post_rma$x[, 2],
  Group  = targets1$Grupos
)

var_explained <- round(100 * pca_post_rma$sdev^2 / sum(pca_post_rma$sdev^2), 1)

# ---- 3. Plot ----
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("CS" = "red", "CTM" = "blue", "N" = "green")) +
  labs(
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)"),
    title = "Exploratory PCA - RMA data (top 5000 most variable genes)",
    color = "Group"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

