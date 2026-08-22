# =============================================================================
#
# Script:      05_gene_annotation.R
#
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
#
# Author:      Karen Rodriguez
#
# Advisor:     Geysson Javier Fernández
#
# Description: Maps the transcript_cluster_id (probe identifiers from the
#              CanGene-1_0-st array) to gene symbols (gene_symbol) using the
#              official Affymetrix/Thermo annotation file, and merges this
#              annotation with the normalized expression matrix (RMA).
#
# Inputs:      data/raw/CanGene-1_0-st-v1.na36.canfam2.probeset.csv
#              data/processed/exprs_rma.rds
#
# Outputs:     results/tables/probes_named_rma_full.csv
#              data/processed/probes_named_rma.rds
#              data/processed/microarray_annot_unique.rds (clean annotation,
#              reused in 10_compare_normalization_methods.R)
#
# Dependencies: 00_setup.R, 02_normalization_rma.R
#
# =============================================================================

source(here::here("R", "00_setup.R"))
library(dplyr)
library(tidyr)

exprs.rma <- readRDS(file.path(PATH_PROCESSED, "exprs_rma.rds"))

# ---- 1. Load official array annotation ----
microarray_annot <- read.delim(
  file.path(PATH_RAW, "CanGene-1_0-st-v1.na36.canfam2.probeset.csv"),
  header = TRUE, sep = ",", stringsAsFactors = FALSE, comment.char = "#"
)
microarray_annot$transcript_cluster_id <- as.character(microarray_annot$transcript_cluster_id)

# ---- 2. Keep one row per transcript_cluster_id + extract gene_symbol ----
microarray_annot_unique <- microarray_annot %>%
  dplyr::group_by(transcript_cluster_id) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    microarray_annot %>% dplyr::count(transcript_cluster_id, name = "n_probesets"),
    by = "transcript_cluster_id"
  ) %>%
  dplyr::mutate(
    gene_symbol = ifelse(
      !is.na(gene_assignment) & gene_assignment != "---",
      sub(".* // ([A-Z0-9]+).*", "\\1", gene_assignment),
      NA
    )
  ) %>%
  dplyr::select(transcript_cluster_id, gene_assignment, n_probesets, gene_symbol, dplyr::everything())

cat("Unique transcript clusters:", length(unique(microarray_annot_unique$transcript_cluster_id)), "\n")

# ---- 3. Merge annotation with expression matrix ----
exprs_df.rma <- as.data.frame(exprs.rma)

# Clean column names (remove everything after "_" and leading zeros)
nombres_limpios          <- sub("_.*", "", colnames(exprs_df.rma))
colnames(exprs_df.rma)   <- as.character(as.numeric(nombres_limpios))
exprs_df.rma$transcript_cluster_id <- rownames(exprs.rma)

probes_named.rma <- merge(
  exprs_df.rma, microarray_annot_unique,
  by = "transcript_cluster_id", all = FALSE  # inner join
)
probes_named.rma <- as_tibble(probes_named.rma)

sample_cols <- names(exprs_df.rma)
sample_cols <- sample_cols[sample_cols != "transcript_cluster_id"]

probes_named.rma <- probes_named.rma %>%
  dplyr::select(transcript_cluster_id, gene_symbol, gene_assignment, n_probesets, all_of(sample_cols)) %>%
  dplyr::filter(!is.na(gene_symbol) & gene_symbol != "")

dim(probes_named.rma)

# ---- 4. Save results ----
write.csv(probes_named.rma, file.path(PATH_TABLES, "probes_named_rma_full.csv"), row.names = FALSE)
saveRDS(probes_named.rma, file.path(PATH_PROCESSED, "probes_named_rma.rds"))
saveRDS(microarray_annot_unique, file.path(PATH_PROCESSED, "microarray_annot_unique.rds"))

cat("Saved: probes_named_rma.rds and microarray_annot_unique.rds\n")

