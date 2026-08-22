# =============================================================================
#
# Script:      09_heatmaps.R
#
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
#
# Author:      Karen Rodriguez
#
# Advisor:     Geysson Javier Fernández
#
# Description: Generates 3 heatmaps of the DEGs (row-scaled / z-score):
#
#              1. N vs CS
#              2. N vs CTM
#              3. Combined N-CS-CTM (with hierarchical clustering of genes,
#                 cut into 4 clusters for export to Prism)
#
# Inputs:      data/processed/expr_matrix_degs.rds
#              data/processed/targets1.rds
#
# Outputs:     results/tables/matriz_prisma_heatmap2_rma.csv
#
# Dependencies: 00_setup.R, 03_load_targets.R, 08_pca_degs.R
#
# =============================================================================

source(here::here("R", "00_setup.R"))
library(pheatmap)
library(RColorBrewer)

degs_matrices <- readRDS(file.path(PATH_PROCESSED, "expr_matrix_degs.rds"))
targets1      <- readRDS(file.path(PATH_PROCESSED, "targets1.rds"))

expr_matrix_rma                <- degs_matrices$expr_matrix_rma
expr_matrix_NCS_filtrado_rma   <- degs_matrices$expr_matrix_NCS_filtrado_rma
expr_matrix_NCTM_filtrado_rma  <- degs_matrices$expr_matrix_NCTM_filtrado_rma

heatmap_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# ---- 1. Heatmap N vs CS ----
numeros_ncs <- as.integer(colnames(expr_matrix_NCS_filtrado_rma))
grupos_ncs  <- targets1$Grupos[match(numeros_ncs, targets1$Renomear)]
matriz_ncs_escalada <- t(scale(t(expr_matrix_NCS_filtrado_rma)))
annotation_ncs <- data.frame(Group = grupos_ncs)
rownames(annotation_ncs) <- colnames(matriz_ncs_escalada)

pheatmap(
  matriz_ncs_escalada,
  annotation_col = annotation_ncs,
  annotation_colors = list(Group = c("N" = "green", "CS" = "red")),
  show_rownames = FALSE, color = heatmap_palette,
  main = "Heatmap RMA - DEGs N vs CS", fontsize_col = 8,
  cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA
)

# ---- 2. Heatmap N vs CTM ----
numeros_nctm <- as.integer(colnames(expr_matrix_NCTM_filtrado_rma))
grupos_nctm  <- targets1$Grupos[match(numeros_nctm, targets1$Renomear)]
matriz_nctm_escalada <- t(scale(t(expr_matrix_NCTM_filtrado_rma)))
annotation_nctm <- data.frame(Group = grupos_nctm)
rownames(annotation_nctm) <- colnames(matriz_nctm_escalada)

pheatmap(
  matriz_nctm_escalada,
  annotation_col = annotation_nctm,
  annotation_colors = list(Group = c("N" = "green", "CTM" = "blue")),
  show_rownames = FALSE, color = heatmap_palette,
  main = "Heatmap RMA - DEGs N vs CTM", fontsize_col = 8,
  cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA
)

# ---- 3. Combined heatmap (3 groups), ordered by group ----
numeros_combinado <- as.integer(colnames(expr_matrix_rma))
grupos_combinado  <- targets1$Grupos[match(numeros_combinado, targets1$Renomear)]

orden_grupos       <- order(factor(grupos_combinado, levels = c("N", "CS", "CTM")))
expr_ordenado       <- expr_matrix_rma[, orden_grupos]
grupos_ordenados    <- grupos_combinado[orden_grupos]
matriz_combinada_escalada <- t(scale(t(expr_ordenado)))

annotation_combinado <- data.frame(Group = grupos_ordenados)
rownames(annotation_combinado) <- colnames(matriz_combinada_escalada)
colores_3grupos <- c("N" = "green", "CS" = "red", "CTM" = "blue")

set.seed(123)
heatmap_NCSCTM_rma <- pheatmap(
  matriz_combinada_escalada,
  annotation_col = annotation_combinado,
  annotation_colors = list(Group = colores_3grupos),
  show_rownames = FALSE, color = heatmap_palette,
  main = "Heatmap RMA - DEGs N vs CS vs CTM", fontsize_col = 8,
  cluster_rows = TRUE, cluster_cols = FALSE,  # keeps group order
  border_color = NA
)

# ---- 4. Extract hierarchical clusters (k=4) for export to Prism ----
dendrograma_filas      <- heatmap_NCSCTM_rma$tree_row
clusters                <- cutree(dendrograma_filas, k = 4)
orden_dendrograma       <- heatmap_NCSCTM_rma$tree_row$order
matriz_ordenada_clusters <- matriz_combinada_escalada[orden_dendrograma, ]
clusters_ordenados       <- clusters[orden_dendrograma]

matriz_prisma_rma <- data.frame(
  Gene = rownames(matriz_ordenada_clusters),
  Cluster = clusters_ordenados,
  matriz_ordenada_clusters,
  row.names = NULL,
  check.names = FALSE
)

cat("Dimensions of matrix for Prism:", dim(matriz_prisma_rma), "\n")
write.csv(
  matriz_prisma_rma,
  file.path(PATH_TABLES, "matriz_prisma_heatmap2_rma.csv"),
  fileEncoding = "UTF-8", row.names = FALSE
)

