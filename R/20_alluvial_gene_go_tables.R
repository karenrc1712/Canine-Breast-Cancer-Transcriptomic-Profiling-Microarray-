# =============================================================================
# Script:      20_alluvial_gene_go_tables.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Hernandez
# Description: Identifies "node" genes -- genes shared across more than one
#              enriched GO BP term -- among the top selected pathways
#              (tabla_stats_up / tabla_stats_down) and builds two long-format
#              tables (gene -> collapsed list of GO processes) suitable as
#              input for alluvial/Sankey plots.
# Inputs:      data/processed/GSVA_analysis/go_gene_tables.rds
#              (uses $up$top_stats and $down$top_stats, i.e. tabla_stats_up/down)
# Outputs:     results/tables/GSVA_analysis/{alluvial_info_up.csv,
#                alluvial_info_down.csv}
#              data/processed/GSVA_analysis/alluvial_info.rds
# Dependencies: 00_setup.R, 19_gsva_pathway_tables.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(dplyr)
library(tidyr)

go_tables <- readRDS(file.path(PATH_GSVA_PROC, "go_gene_tables.rds"))

tabla_stats_up   <- go_tables$up$top_stats
tabla_stats_down <- go_tables$down$top_stats

# ---- Helper: explode "Genes" column, find node genes, collapse GO terms ----
# A "node gene" is a gene that appears in more than one GO_Term row within
# the same direction (up or down) table. For each node gene, its associated
# GO terms are collapsed into a single semicolon-separated string, ready for
# use as the two columns of an alluvial/Sankey diagram (Gene -> Process).
build_alluvial_table <- function(tabla_stats) {
  
  gene_go_long <- tabla_stats %>%
    dplyr::select(GO_ID, GO_Term, Genes) %>%
    tidyr::separate_rows(Genes, sep = ",\\s*") %>%
    dplyr::rename(Gene = Genes) %>%
    dplyr::filter(Gene != "") %>%
    dplyr::distinct(Gene, GO_ID, GO_Term)
  
  node_genes <- gene_go_long %>%
    dplyr::count(Gene, name = "n_process") %>%
    dplyr::filter(n_process > 1) %>%
    dplyr::pull(Gene)
  
  cat("Node genes found:", length(node_genes), "\n")
  
  gene_go_long %>%
    dplyr::filter(Gene %in% node_genes) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      N_process = dplyr::n(),
      Process = paste(sort(unique(GO_Term)), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(N_process))
}
# ============================================================
# UPREGULATED node genes
# ============================================================
cat("--- Building alluvial table: UPREGULATED ---\n")
alluvial_info_up <- build_alluvial_table(tabla_stats_up)

# ============================================================
# DOWNREGULATED node genes
# ============================================================
cat("--- Building alluvial table: DOWNREGULATED ---\n")
alluvial_info_down <- build_alluvial_table(tabla_stats_down)

# ---- Export ----
write.csv(alluvial_info_up,   file.path(PATH_GSVA_TAB, "alluvial_info_up.csv"),   row.names = FALSE)
write.csv(alluvial_info_down, file.path(PATH_GSVA_TAB, "alluvial_info_down.csv"), row.names = FALSE)

saveRDS(
  list(up = alluvial_info_up, down = alluvial_info_down),
  file.path(PATH_GSVA_PROC, "alluvial_info.rds")
)

cat("Saved: alluvial_info_up.csv, alluvial_info_down.csv, alluvial_info.rds\n")

