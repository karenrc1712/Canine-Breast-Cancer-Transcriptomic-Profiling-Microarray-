# =============================================================================
#
# Script:      08_pca_degs.R
#
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
#
# Author:      Karen Rodriguez
#
# Advisor:     Geysson Javier Hernandez
#
# Description: Builds an expression matrix filtered to include only DEGs
#              (combining N vs CS and N vs CTM, retaining the probeset with
#              the highest variance per gene) and performs PCA: one combined
#              analysis (3 groups) and two individual analyses (N vs CS,
#              N vs CTM). Exports the coordinates for plotting in Prism as
#              well.
#
# Inputs:      data/processed/probes_named_rma.rds
#              data/processed/degs_rma.rds
#              data/processed/targets1.rds
#              data/processed/exprs_matrix_rma.rds (for sample_cols_rma)
#
# Outputs:     results/tables/PCA_NCS_RMA.csv
#              results/tables/PCA_NCTM_RMA.csv
#              results/tables/pca_data_combinado_rma.csv
#              data/processed/expr_matrix_degs.rds (reused for heatmaps)
#
# Dependencies: 00_setup.R, 05_gene_annotation.R, 07_extract_degs.R
#
# =============================================================================

source(here::here("R", "00_setup.R"))
library(dplyr)
library(ggplot2)
library(ggrepel)

probes_named.rma <- readRDS(file.path(PATH_PROCESSED, "probes_named_rma.rds"))
degs             <- readRDS(file.path(PATH_PROCESSED, "degs_rma.rds"))
targets1         <- readRDS(file.path(PATH_PROCESSED, "targets1.rds"))
sample_cols_rma  <- readRDS(file.path(PATH_PROCESSED, "exprs_matrix_rma.rds"))$sample_cols

DEGs_NCS_rma  <- degs$NCS$total
DEGs_NCTM_rma <- degs$NCTM$total

# ---- 1. Combine DEG genes from both contrasts (remove duplicates) ----
combined_genes_rma <- bind_rows(DEGs_NCS_rma, DEGs_NCTM_rma) %>%
  group_by(gene_symbol) %>%
  slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
  ungroup()

selected_genes_rma <- unique(combined_genes_rma$gene_symbol)
cat("Unique combined genes (DEGs):", length(selected_genes_rma), "\n")

# ---- 2. Filter matrix, keeping the probeset with the highest variance per gene ----
expr_filtered_rma <- probes_named.rma %>%
  filter(gene_symbol %in% selected_genes_rma, !is.na(gene_symbol), gene_symbol != "") %>%
  mutate(across(all_of(sample_cols_rma), as.numeric)) %>%
  rowwise() %>%
  mutate(row_variance = var(c_across(all_of(sample_cols_rma)), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(row_variance)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  dplyr::select(-row_variance)

expr_matrix_rma <- as.matrix(expr_filtered_rma[, sample_cols_rma])
rownames(expr_matrix_rma) <- expr_filtered_rma$gene_symbol
cat("Final DEG matrix:", dim(expr_matrix_rma), "\n")

numeros_columna_rma <- as.integer(sub("_.*", "", colnames(expr_matrix_rma)))

# ---- 3. Combined PCA (3 groups: N, CS, CTM) ----
pca_combinado_rma <- prcomp(t(expr_matrix_rma), scale. = TRUE)
numeros_muestras_rma <- as.integer(sub("_.*", "", rownames(pca_combinado_rma$x)))

pca_data_combinado_rma <- data.frame(
  sample = rownames(pca_combinado_rma$x),
  PC1 = pca_combinado_rma$x[, 1],
  PC2 = pca_combinado_rma$x[, 2],
  Group = targets1$Grupos[match(numeros_muestras_rma, targets1$Renomear)],
  row.names = NULL
)

pca_var_per_comb_rma <- round(pca_combinado_rma$sdev^2 / sum(pca_combinado_rma$sdev^2) * 100, 1)
cat("PC1+PC2 variance (combined):", sum(pca_var_per_comb_rma[1:2]), "%\n")

ggplot(pca_data_combinado_rma, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 15, show.legend = FALSE) +
  scale_color_manual(
    name = "Groups", values = c("CS" = "red", "CTM" = "blue", "N" = "green"),
    labels = c("CS", "CTM", "Normal")
  ) +
  xlab(paste0("PC1 (", pca_var_per_comb_rma[1], "% of variance)")) +
  ylab(paste0("PC2 (", pca_var_per_comb_rma[2], "% of variance)")) +
  ggtitle("PCA RMA: CS vs CTM vs N") +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))

write.csv(pca_data_combinado_rma, file.path(PATH_TABLES, "pca_data_combinado_rma.csv"), row.names = FALSE)

# ---- Helper function: Individual PCA for a pair of groups ----
run_pca_pair <- function(degs_table, group_a, group_b, color_map, title) {
  genes_pair <- unique(degs_table$gene_symbol)
  
  expr_pair <- expr_filtered_rma %>% filter(gene_symbol %in% genes_pair)
  expr_matrix_pair <- as.matrix(expr_pair[, sample_cols_rma])
  rownames(expr_matrix_pair) <- expr_pair$gene_symbol
  
  muestras_pair  <- targets1$Renomear[targets1$Grupos %in% c(group_a, group_b)]
  columnas_pair  <- which(numeros_columna_rma %in% muestras_pair)
  expr_matrix_pair_filt <- expr_matrix_pair[, columnas_pair, drop = FALSE]
  
  pca_pair <- prcomp(t(expr_matrix_pair_filt), scale. = TRUE)
  numeros_muestras_pair <- as.integer(sub("_.*", "", rownames(pca_pair$x)))
  
  pca_data_pair <- data.frame(
    sample = rownames(pca_pair$x),
    PC1 = pca_pair$x[, 1],
    PC2 = pca_pair$x[, 2],
    Group = targets1$Grupos[match(numeros_muestras_pair, targets1$Renomear)],
    row.names = NULL
  )
  
  var_pair <- round(pca_pair$sdev^2 / sum(pca_pair$sdev^2) * 100, 1)
  
  p <- ggplot(pca_data_pair, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(aes(label = sample), size = 3, max.overlaps = 20) +
    scale_color_manual(values = color_map) +
    xlab(paste0("PC1 (", var_pair[1], "% of variance)")) +
    ylab(paste0("PC2 (", var_pair[2], "% of variance)")) +
    ggtitle(title) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
  list(pca_data = pca_data_pair, expr_matrix_filtrado = expr_matrix_pair_filt)
}

# ---- 4. PCA N vs CS ----
result_ncs <- run_pca_pair(
  DEGs_NCS_rma, "N", "CS",
  c("CS" = "red", "N" = "green"),
  "PCA RMA: N vs CS"
)
write.csv(result_ncs$pca_data, file.path(PATH_TABLES, "PCA_NCS_RMA.csv"), row.names = FALSE, na = "")

# ---- 5. PCA N vs CTM ----
result_nctm <- run_pca_pair(
  DEGs_NCTM_rma, "N", "CTM",
  c("CTM" = "blue", "N" = "green"),
  "PCA RMA: N vs CTM"
)
write.csv(result_nctm$pca_data, file.path(PATH_TABLES, "PCA_NCTM_RMA.csv"), row.names = FALSE, na = "")

# ---- 6. Save matrices for the heatmap script ----
saveRDS(
  list(
    expr_matrix_rma = expr_matrix_rma,
    numeros_columna_rma = numeros_columna_rma,
    expr_matrix_NCS_filtrado_rma = result_ncs$expr_matrix_filtrado,
    expr_matrix_NCTM_filtrado_rma = result_nctm$expr_matrix_filtrado
  ),
  file.path(PATH_PROCESSED, "expr_matrix_degs.rds")
)
cat("Saved: data/processed/expr_matrix_degs.rds\n")

