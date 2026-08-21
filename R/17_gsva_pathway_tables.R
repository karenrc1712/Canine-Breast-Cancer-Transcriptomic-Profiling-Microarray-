# =============================================================================
# Script:      17_gsva_pathway_tables.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Hernandez
# Description: Selects the top differentially enriched GO BP terms for both
#              upregulated (logFC > 0) and downregulated (logFC < 0) pathways
#              across the CS vs N and CTM vs N contrasts. For each direction:
#                - Identifies GO terms significant in either contrast
#                  (adj.P.Val < 0.05)
#                - Selects the top 30 terms per contrast ranked by |logFC|
#                - Annotates with GO term names, constituent gene lists, and
#                  mean group GSVA scores
#              All tables are exported as CSV files.
# Inputs:      data/processed/GSVA_analysis/gsva_topTables.rds
#              data/processed/GSVA_analysis/genesbygo_filtered.rds
#              data/processed/GSVA_analysis/gsva_r.rds
#              data/processed/GSVA_analysis/gsva_targets_aligned.rds
# Outputs:     results/tables/GSVA_analysis/{go_gene_table_up/down.csv,
#                tabla_stats_up/down.csv, gsva_group_means_up/down.csv}
#              data/processed/GSVA_analysis/go_gene_tables.rds
# Dependencies: 00_setup.R, 16_gsva_go_genesets.R, 18_gsva_differential_enrichment.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(GO.db)
library(dplyr)

top_tables     <- readRDS(file.path(PATH_GSVA_PROC, "gsva_topTables.rds"))
genesbygo      <- readRDS(file.path(PATH_GSVA_PROC, "genesbygo_filtered.rds"))
gsva_r         <- readRDS(file.path(PATH_GSVA_PROC, "gsva_r.rds"))
targets_gsva   <- readRDS(file.path(PATH_GSVA_PROC, "gsva_targets_aligned.rds"))

res_CS_vs_N  <- top_tables$CS_vs_N
res_CTM_vs_N <- top_tables$CTM_vs_N

# ---- Helper: build annotated GO table for a set of GO IDs ----
build_go_annotation_table <- function(go_ids, res_cs, res_ctm, genesets) {
  # Keep only gene sets present in the filtered list
  valid_ids <- intersect(go_ids, names(genesets))
  gene_sets  <- genesets[valid_ids]
  
  data.frame(
    GO_ID     = names(gene_sets),
    GO_Term   = Term(GOTERM[names(gene_sets)]),
    N_genes   = sapply(gene_sets, length),
    Genes     = sapply(gene_sets, paste, collapse = ", "),
    FC_CS     = res_cs[names(gene_sets),   "logFC"],
    adjpval_CS  = res_cs[names(gene_sets), "adj.P.Val"],
    FC_CTM    = res_ctm[names(gene_sets),  "logFC"],
    adjpval_CTM = res_ctm[names(gene_sets),"adj.P.Val"],
    stringsAsFactors = FALSE, row.names = NULL
  )
}

# ---- Helper: select top 30 terms per contrast by |logFC| ----
top_go_ids <- function(res_cs, res_ctm, direction = c("up", "down"), n = 30) {
  direction <- match.arg(direction)
  if (direction == "up") {
    ids_cs  <- rownames(res_cs  %>% filter(adj.P.Val < 0.05, logFC >  0) %>%
                          arrange(desc(logFC))    %>% head(n))
    ids_ctm <- rownames(res_ctm %>% filter(adj.P.Val < 0.05, logFC >  0) %>%
                          arrange(desc(logFC))    %>% head(n))
  } else {
    ids_cs  <- rownames(res_cs  %>% filter(adj.P.Val < 0.05, logFC <  0) %>%
                          arrange(logFC)          %>% head(n))
    ids_ctm <- rownames(res_ctm %>% filter(adj.P.Val < 0.05, logFC <  0) %>%
                          arrange(logFC)          %>% head(n))
  }
  unique(c(ids_cs, ids_ctm))
}

# ---- Helper: compute per-group mean GSVA scores ----
group_means <- function(go_ids, gsva_matrix, grupos) {
  mat    <- gsva_matrix[intersect(go_ids, rownames(gsva_matrix)), , drop = FALSE]
  result <- sapply(unique(grupos), function(g)
    rowMeans(mat[, grupos == g, drop = FALSE]))
  as.data.frame(result)
}

grupos_named <- targets_gsva$Grupos
names(grupos_named) <- colnames(gsva_r)

# ============================================================
# UPREGULATED pathways
# ============================================================
cat("--- Processing UPREGULATED GO terms ---\n")

up_go_ids   <- unique(c(
  rownames(res_CS_vs_N  %>% filter(adj.P.Val < 0.05, logFC > 0)),
  rownames(res_CTM_vs_N %>% filter(adj.P.Val < 0.05, logFC > 0))
))
cat("  GO terms upregulated in CS or CTM vs N:", length(up_go_ids), "\n")

go_table_up   <- build_go_annotation_table(up_go_ids,  res_CS_vs_N, res_CTM_vs_N, genesbygo)
means_up      <- group_means(up_go_ids, gsva_r, grupos_named)
top_ids_up    <- top_go_ids(res_CS_vs_N, res_CTM_vs_N, "up")
tabla_stats_up <- build_go_annotation_table(top_ids_up, res_CS_vs_N, res_CTM_vs_N, genesbygo)

cat("  Top selected terms (union of top 30 per contrast):", nrow(tabla_stats_up), "\n")

write.csv(go_table_up,   file.path(PATH_GSVA_TAB, "go_gene_table_up.csv"),  row.names = FALSE)
write.csv(means_up,      file.path(PATH_GSVA_TAB, "gsva_group_means_up.csv"))
write.csv(tabla_stats_up,file.path(PATH_GSVA_TAB, "tabla_stats_up.csv"),    row.names = FALSE)

# ============================================================
# DOWNREGULATED pathways
# ============================================================
cat("--- Processing DOWNREGULATED GO terms ---\n")

down_go_ids   <- unique(c(
  rownames(res_CS_vs_N  %>% filter(adj.P.Val < 0.05, logFC < 0)),
  rownames(res_CTM_vs_N %>% filter(adj.P.Val < 0.05, logFC < 0))
))
cat("  GO terms downregulated in CS or CTM vs N:", length(down_go_ids), "\n")

go_table_down   <- build_go_annotation_table(down_go_ids, res_CS_vs_N, res_CTM_vs_N, genesbygo)
means_down      <- group_means(down_go_ids, gsva_r, grupos_named)
top_ids_down    <- top_go_ids(res_CS_vs_N, res_CTM_vs_N, "down")
tabla_stats_down <- build_go_annotation_table(top_ids_down, res_CS_vs_N, res_CTM_vs_N, genesbygo)

cat("  Top selected terms (union of top 30 per contrast):", nrow(tabla_stats_down), "\n")

write.csv(go_table_down,   file.path(PATH_GSVA_TAB, "go_gene_table_down.csv"),  row.names = FALSE)
write.csv(means_down,      file.path(PATH_GSVA_TAB, "gsva_group_means_down.csv"))
write.csv(tabla_stats_down,file.path(PATH_GSVA_TAB, "tabla_stats_down.csv"),    row.names = FALSE)

# ---- Save combined object ----
saveRDS(
  list(
    up   = list(go_table = go_table_up,   top_stats = tabla_stats_up,
                means = means_up,         go_ids = up_go_ids),
    down = list(go_table = go_table_down, top_stats = tabla_stats_down,
                means = means_down,       go_ids = down_go_ids)
  ),
  file.path(PATH_GSVA_PROC, "go_gene_tables.rds")
)
cat("Saved: go_gene_tables.rds\n")

