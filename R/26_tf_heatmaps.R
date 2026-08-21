# =============================================================================
# Script:      26_tf_heatmaps.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Hernandez
# Description: Generates two z-score scaled pheatmap figures for TFs that are
#              differentially expressed in either the NCS or NCTM contrast:
#                - UP heatmap:   TFs upregulated in CS and/or CTM vs N
#                - DOWN heatmap: TFs downregulated in CS and/or CTM vs N
#              Samples are ordered N → CS → CTM (no column clustering).
#              Row clustering is hierarchical (Euclidean / complete linkage);
#              genes are cut into 4 clusters for export to Prism.
#              Both heatmaps exported as vectorial PDFs and Prism-ready CSVs.
# Inputs:      data/processed/GSVA_analysis/expresion_final.rds
#              data/processed/targets1.rds
#              data/processed/TFs/tf_secreted_intersections.rds
# Outputs:     results/figures/TFs/heatmap_TF_{up,down}.pdf
#              results/tables/TFs/prism_heatmap_TF_{up,down}.csv
# Dependencies: 00_setup.R, 13_gsva_prepare_expression.R, 26 (self)
# =============================================================================

source(here::here("R", "00_setup.R"))
library(dplyr)
library(pheatmap)
library(RColorBrewer)

PATH_TF_PROC <- file.path(PATH_PROCESSED, "TFs")
PATH_TF_TAB  <- file.path(PATH_TABLES,    "TFs")
PATH_TF_FIG  <- file.path(PATH_FIGURES,   "TFs")

expresion_final <- readRDS(file.path(PATH_GSVA_PROC, "expresion_final.rds"))
targets1        <- readRDS(file.path(PATH_PROCESSED, "targets1.rds"))
tf_int          <- readRDS(file.path(PATH_TF_PROC,   "tf_intersections.rds"))

# ---- 1. Align column names between expression matrix and targets ----
# expresion_final columns are numeric IDs; targets1$Renomear uses zero-padded strings
colnames(expresion_final) <- as.character(colnames(expresion_final))
targets1$Renomear         <- sprintf("%02d", as.numeric(targets1$Renomear))
targets1                  <- targets1 %>% filter(Renomear %in% colnames(expresion_final))

# Sample IDs per group (in the order they appear in expresion_final)
samples_by_group <- function(grupo)
  targets1 %>% filter(Grupos == grupo) %>% pull(Renomear)

s_N   <- samples_by_group("N")
s_CS  <- samples_by_group("CS")
s_CTM <- samples_by_group("CTM")
col_order <- c(s_N, s_CS, s_CTM)   # canonical left-to-right order

group_labels <- c(rep("N", length(s_N)),
                  rep("CS", length(s_CS)),
                  rep("CTM", length(s_CTM)))
names(group_labels) <- col_order

# ---- 2. Heatmap factory function ----
# Receives a TF gene list, builds the z-score matrix, plots, and exports.
make_tf_heatmap <- function(tf_genes, direction, n_clusters = 4) {

  # Keep TFs present in the expression matrix
  tf_present <- intersect(tf_genes, rownames(expresion_final))
  cat(direction, ": using", length(tf_present), "of", length(tf_genes),
      "TFs found in expression matrix\n")
  if (length(tf_present) < 2) {
    warning("Too few TFs for a heatmap (", length(tf_present), "). Skipping.")
    return(invisible(NULL))
  }

  # Build ordered matrix
  expr_mat <- expresion_final[tf_present, col_order]
  z_mat    <- t(scale(t(expr_mat)))   # z-score per gene (row-wise)

  # Annotation
  ann_col <- data.frame(Group = group_labels)
  rownames(ann_col) <- col_order

  # Plot
  hm <- pheatmap::pheatmap(
    z_mat,
    annotation_col    = ann_col,
    annotation_colors = list(Group = c("N" = "#2ECC71", "CS" = "#E74C3C", "CTM" = "#2980B9")),
    show_rownames     = FALSE,
    color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    main              = paste0("TFs ", direction, "-regulated — N vs CS vs CTM (z-score)"),
    fontsize_col      = 8,
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    border_color      = NA,
    silent            = TRUE
  )

  # Export vectorial PDF
  pdf_path <- file.path(PATH_TF_FIG, paste0("heatmap_TF_", tolower(direction), ".pdf"))
  cairo_pdf(pdf_path, width = 8, height = max(5, length(tf_present) * 0.22),
            family = "Helvetica")
  grid::grid.newpage()
  grid::grid.draw(hm$gtable)
  dev.off()
  cat("  Saved:", pdf_path, "\n")

  # Export Prism CSV with cluster assignments
  ord     <- hm$tree_row$order
  clust   <- cutree(hm$tree_row, k = n_clusters)
  prisma  <- data.frame(
    Gene    = rownames(z_mat)[ord],
    Cluster = clust[ord],
    z_mat[ord, ],
    row.names = NULL, check.names = FALSE
  )
  csv_path <- file.path(PATH_TF_TAB, paste0("prism_heatmap_TF_", tolower(direction), ".csv"))
  write.csv(prisma, csv_path, row.names = FALSE, fileEncoding = "UTF-8")
  cat("  Saved:", csv_path, "\n")

  invisible(hm)
}

# ---- 3. UP heatmap (TFs upregulated in CS or CTM vs N) ----
tf_up_combined <- unique(c(tf_int$NCS$TF_up, tf_int$NCTM$TF_up))
make_tf_heatmap(tf_up_combined, "UP")

# ---- 4. DOWN heatmap (TFs downregulated in CS or CTM vs N) ----
tf_down_combined <- unique(c(tf_int$NCS$TF_down, tf_int$NCTM$TF_down))
make_tf_heatmap(tf_down_combined, "DOWN")

cat("\nTF heatmaps complete.\n")

