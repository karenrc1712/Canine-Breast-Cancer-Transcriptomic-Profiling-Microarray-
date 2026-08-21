# =============================================================================
#
# Script:      07_extract_degs.R
#
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
#
# Author:      Karen Rodriguez
#
# Advisor:     Geysson Javier Hernandez
#
# Description: Starting from the top tables, filters differentially expressed
#              genes (DEGs) for each contrast using the following thresholds:
#
#              |logFC| >= 1  and  adj.P.Val < 0.05
#
#              separating them into overexpressed (up) and underexpressed
#              (down) genes.
#
#              It also identifies "contradictory" genes (those appearing as
#              both up- and downregulated due to having multiple probes).
#
# Inputs:      data/processed/toptables_rma.rds
#              data/processed/exprs_matrix_rma.rds  (for genes_info)
#
# Outputs:     results/tables/DEGs_*.csv (8 files: total/up/down x 3
#              contrasts + total for CSCTM)
#              data/processed/degs_rma.rds (list containing all DEGs)
#
# Dependencies: 00_setup.R, 06_differential_expression.R
#
# =============================================================================

source(here::here("R", "00_setup.R"))
library(dplyr)

toptables  <- readRDS(file.path(PATH_PROCESSED, "toptables_rma.rds"))
exprs_info <- readRDS(file.path(PATH_PROCESSED, "exprs_matrix_rma.rds"))
genes_info_rma <- exprs_info$genes_info

# ---- Helper function: filters up/down/total DEGs for a contrast ----
extract_degs <- function(toptable, genes_info, logfc_threshold = 1, padj_threshold = 0.05) {
  base <- toptable %>%
    dplyr::mutate(transcript_cluster_id = rownames(.)) %>%
    dplyr::left_join(genes_info, by = "transcript_cluster_id") %>%
    dplyr::relocate(transcript_cluster_id, gene_symbol)
  
  list(
    up    = base %>% dplyr::filter(logFC >=  logfc_threshold, adj.P.Val < padj_threshold),
    down  = base %>% dplyr::filter(logFC <= -logfc_threshold, adj.P.Val < padj_threshold),
    total = base %>% dplyr::filter(abs(logFC) >= logfc_threshold, adj.P.Val < padj_threshold)
  )
}

# ---- Apply to the 3 contrasts ----
degs_NCS   <- extract_degs(toptables$NCS,   genes_info_rma)
degs_NCTM  <- extract_degs(toptables$NCTM,  genes_info_rma)
degs_CSCTM <- extract_degs(toptables$CSCTM, genes_info_rma)

cat("=== N vs CS ===\n  Up:", nrow(degs_NCS$up), " Down:", nrow(degs_NCS$down), " Total:", nrow(degs_NCS$total), "\n")
cat("=== N vs CTM ===\n  Up:", nrow(degs_NCTM$up), " Down:", nrow(degs_NCTM$down), " Total:", nrow(degs_NCTM$total), "\n")
cat("=== CS vs CTM ===\n  Up:", nrow(degs_CSCTM$up), " Down:", nrow(degs_CSCTM$down), " Total:", nrow(degs_CSCTM$total), "\n")

# ---- Contradictory genes (appear as both up AND down due to multiple probes) ----
contradictorios_NCS  <- intersect(degs_NCS$up$gene_symbol,  degs_NCS$down$gene_symbol)
contradictorios_NCTM <- intersect(degs_NCTM$up$gene_symbol, degs_NCTM$down$gene_symbol)
cat("Contradictory genes N vs CS:",  length(contradictorios_NCS),  "\n")
cat("Contradictory genes N vs CTM:", length(contradictorios_NCTM), "\n")

# ---- Export CSVs ----
write.csv(degs_NCS$total,   file.path(PATH_TABLES, "DEGs_NCS_rma.csv"),       row.names = FALSE)
write.csv(degs_NCS$up,      file.path(PATH_TABLES, "DEGs_up_NCS_rma.csv"),    row.names = FALSE)
write.csv(degs_NCS$down,    file.path(PATH_TABLES, "DEGs_down_NCS_rma.csv"),  row.names = FALSE)
write.csv(degs_NCTM$total,  file.path(PATH_TABLES, "DEGs_NCTM_rma.csv"),      row.names = FALSE)
write.csv(degs_NCTM$up,     file.path(PATH_TABLES, "DEGs_up_NCTM_rma.csv"),   row.names = FALSE)
write.csv(degs_NCTM$down,   file.path(PATH_TABLES, "DEGs_down_NCTM_rma.csv"), row.names = FALSE)
write.csv(degs_CSCTM$total, file.path(PATH_TABLES, "DEGs_CSCTM_rma.csv"),     row.names = FALSE)

# ---- Save everything together for subsequent scripts ----
saveRDS(
  list(NCS = degs_NCS, NCTM = degs_NCTM, CSCTM = degs_CSCTM),
  file.path(PATH_PROCESSED, "degs_rma.rds")
)
cat("Saved: data/processed/degs_rma.rds\n")

