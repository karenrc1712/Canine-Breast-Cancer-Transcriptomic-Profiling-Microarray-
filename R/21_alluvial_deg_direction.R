# =============================================================================
# Script:      21_alluvial_deg_direction.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Fernández
# Description: For each node gene in alluvial_info_up / alluvial_info_down,
#              checks whether it is DE (and in which contrast: NCS and/or
#              NCTM) consistent with its direction -- i.e. up-node genes are
#              cross-checked ONLY against DEGs_up_NCS_rma / DEGs_up_NCTM_rma,
#              and down-node genes ONLY against DEGs_down_NCS_rma /
#              DEGs_down_NCTM_rma. Adds contrast membership + logFC/adj.P.Val
#              per contrast to the alluvial tables.
# Inputs:      data/processed/GSVA_analysis/alluvial_info.rds
#              PATH_TABLES/DEGs_up_NCS_rma.csv
#              PATH_TABLES/DEGs_up_NCTM_rma.csv
#              PATH_TABLES/DEGs_down_NCS_rma.csv
#              PATH_TABLES/DEGs_down_NCTM_rma.csv
# Outputs:     results/tables/GSVA_analysis/{alluvial_info_up.csv,
#                alluvial_info_down.csv}  (overwritten, now annotated)
#              data/processed/GSVA_analysis/alluvial_info_annotated.rds
# Dependencies: 00_setup.R, 20_alluvial_gene_go_tables.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(dplyr)

alluvial_info <- readRDS(file.path(PATH_GSVA_PROC, "alluvial_info.rds"))
alluvial_info_up   <- alluvial_info$up
alluvial_info_down <- alluvial_info$down

# ---- Helper: read a DEGs_* file and standardize the gene-id column ----
# These files are CSV ("Excel comma separated"). Depending on how they were
# exported from limma::topTable(), the gene symbol may either have its own
# header (e.g. "Gene", "Symbol", "Gene_symbol") or be an unnamed row-name
# column (imported by read.csv as "X"). This helper detects and standardizes
# it to a column called "Gene".
read_deg_file <- function(path, contrast_label) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  
  stopifnot("gene_symbol" %in% colnames(df))
  
  df <- df %>%
    dplyr::rename(Gene = gene_symbol) %>%
    dplyr::select(Gene, logFC, adj.P.Val)
  
  colnames(df)[colnames(df) == "logFC"]     <- paste0("logFC_", contrast_label)
  colnames(df)[colnames(df) == "adj.P.Val"] <- paste0("adjP_", contrast_label)
  df
}
# ---- Helper: annotate an alluvial table with contrast membership ----
annotate_direction <- function(alluvial_tbl, deg_cs, deg_ctm) {
  annotated <- alluvial_tbl %>%
    left_join(deg_cs,  by = "Gene") %>%
    left_join(deg_ctm, by = "Gene")
  
  logfc_cs_col  <- grep("^logFC_", colnames(deg_cs),  value = TRUE)
  logfc_ctm_col <- grep("^logFC_", colnames(deg_ctm), value = TRUE)
  
  annotated <- annotated %>%
    mutate(
      DE_in_CS  = !is.na(.data[[logfc_cs_col]]),
      DE_in_CTM = !is.na(.data[[logfc_ctm_col]]),
      Contrast  = case_when(
        DE_in_CS & DE_in_CTM  ~ "CS_and_CTM",
        DE_in_CS              ~ "CS_only",
        DE_in_CTM             ~ "CTM_only",
        TRUE                  ~ "neither"   # shouldn't happen if genes came from these tables
      )
    )
  
  n_neither <- sum(annotated$Contrast == "neither")
  if (n_neither > 0) {
    warning(n_neither, " node gene(s) not matched to any DEG file -- ",
            "check gene ID formatting (symbol vs alias vs case).")
  }
  
  annotated
}

# ============================================================
# UPREGULATED node genes -- match against UP DEG files only
# ============================================================
cat("--- Annotating UPREGULATED node genes with DEG direction ---\n")

deg_up_cs  <- read_deg_file(file.path(PATH_TABLES, "DEGs_up_NCS_rma.csv"),  "CS")
deg_up_ctm <- read_deg_file(file.path(PATH_TABLES, "DEGs_up_NCTM_rma.csv"), "CTM")

alluvial_info_up <- annotate_direction(alluvial_info_up, deg_up_cs, deg_up_ctm)

# ============================================================
# DOWNREGULATED node genes -- match against DOWN DEG files only
# ============================================================
cat("--- Annotating DOWNREGULATED node genes with DEG direction ---\n")

deg_down_cs  <- read_deg_file(file.path(PATH_TABLES, "DEGs_down_NCS_rma.csv"),  "CS")
deg_down_ctm <- read_deg_file(file.path(PATH_TABLES, "DEGs_down_NCTM_rma.csv"), "CTM")

alluvial_info_down <- annotate_direction(alluvial_info_down, deg_down_cs, deg_down_ctm)

# ---- Export (overwrite previous, now with DE annotation) ----
write.csv(alluvial_info_up,   file.path(PATH_GSVA_TAB, "alluvial_info_up.csv"),   row.names = FALSE)
write.csv(alluvial_info_down, file.path(PATH_GSVA_TAB, "alluvial_info_down.csv"), row.names = FALSE)

saveRDS(
  list(up = alluvial_info_up, down = alluvial_info_down),
  file.path(PATH_GSVA_PROC, "alluvial_info_annotated.rds")
)

cat("Saved: alluvial_info_up.csv, alluvial_info_down.csv, alluvial_info_annotated.rds\n")

