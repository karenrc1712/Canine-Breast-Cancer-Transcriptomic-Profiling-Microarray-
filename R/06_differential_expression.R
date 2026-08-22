# =============================================================================
#
# Script:      06_differential_expression.R
#
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
#
# Author:      Karen Rodriguez
#
# Advisor:     Geysson Javier Fernández
#
# Description: Fits a linear model (limma) to the annotated matrix and
#              calculates differential expression for 3 contrasts:
#
#              NCS   = CS - N     (simple carcinoma vs normal)
#              NCTM  = CTM - N    (tubulomedullary carcinoma vs normal)
#              CSCTM = CS - CTM   (between the two carcinoma types)
#
#              Applies empirical Bayes moderation (eBayes) and exports the
#              top tables.
#
# Inputs:      data/processed/probes_named_rma.rds
#              data/processed/targets1.rds
#
# Outputs:     results/tables/toptable_{NCS,NCTM,CSCTM}_rma.csv
#              data/processed/fit2_rma.rds
#              data/processed/toptables_rma.rds (list containing the 3 top tables)
#              data/processed/exprs_matrix_rma.rds (numeric matrix + genes_info)
#
# Dependencies: 00_setup.R, 03_load_targets.R, 05_gene_annotation.R
#
# =============================================================================

source(here::here("R", "00_setup.R"))
library(limma)
library(dplyr)

probes_named.rma <- readRDS(file.path(PATH_PROCESSED, "probes_named_rma.rds"))
targets1         <- readRDS(file.path(PATH_PROCESSED, "targets1.rds"))
microarray_annot_unique <- readRDS(file.path(PATH_PROCESSED, "microarray_annot_unique.rds"))

# ---- 1. Identify sample columns and verify against targets ----
sample_cols_rma <- names(probes_named.rma)[grepl("^[0-9]+$", names(probes_named.rma))]
cel_ids_rma     <- as.integer(sample_cols_rma)
targets_ordered <- targets1[match(cel_ids_rma, targets1$Renomear), ]

stopifnot(
  "Samples do not match"  = all(targets_ordered$Renomear == cel_ids_rma),
  "Missing targets"       = nrow(targets_ordered) == length(sample_cols_rma)
)

# ---- 2. Define groups and design matrix ----
G <- factor(targets1$Grupos)
cat("Defined groups:", levels(G), "\n")

MD_rma <- model.matrix(~0 + G, data = targets1)
colnames(MD_rma) <- levels(G)
print(MD_rma)

# ---- 3. Prepare numeric expression matrix ----
exprs_mat_rma <- probes_named.rma %>%
  dplyr::select(gene_symbol, transcript_cluster_id, all_of(sample_cols_rma))

genes_info_rma <- exprs_mat_rma %>%
  dplyr::select(transcript_cluster_id, gene_symbol) %>%
  dplyr::mutate(transcript_cluster_id = as.character(transcript_cluster_id))

exprs_matrix_rma <- exprs_mat_rma %>%
  dplyr::select(all_of(sample_cols_rma)) %>%
  as.matrix()
rownames(exprs_matrix_rma) <- genes_info_rma$transcript_cluster_id

stopifnot(nrow(exprs_matrix_rma) == nrow(genes_info_rma))

# ---- 4. Linear fit + contrasts + eBayes ----
fit_rma <- lmFit(exprs_matrix_rma, MD_rma)

cont.matrix_rma <- makeContrasts(
  NCS   = CS - N,
  NCTM  = CTM - N,
  CSCTM = CS - CTM,
  levels = MD_rma
)
print(cont.matrix_rma)

fit2_rma <- contrasts.fit(fit_rma, cont.matrix_rma)
fit2_rma <- eBayes(fit2_rma)

fit_results_rma <- data.frame(
  transcript_cluster_id = rownames(exprs_matrix_rma),
  fit2_rma$coefficients
) %>%
  dplyr::left_join(genes_info_rma, by = "transcript_cluster_id") %>%
  dplyr::relocate(transcript_cluster_id, gene_symbol)

# ---- 5. Extract toptables per contrast ----
toptable_NCS_rma   <- topTable(fit2_rma, coef = 1, number = Inf, sort.by = "none")
toptable_NCTM_rma  <- topTable(fit2_rma, coef = 2, number = Inf, sort.by = "none")
toptable_CSCTM_rma <- topTable(fit2_rma, coef = 3, number = Inf, sort.by = "none")

for (tbl_name in c("toptable_NCS_rma", "toptable_NCTM_rma", "toptable_CSCTM_rma")) {
  tbl <- get(tbl_name)
  rownames(tbl) <- rownames(exprs_matrix_rma)
  tbl$transcript_cluster_id <- rownames(tbl)
  assign(tbl_name, tbl)
}

toptable_NCS_final   <- toptable_NCS_rma   %>% dplyr::inner_join(microarray_annot_unique, by = "transcript_cluster_id") %>% dplyr::select(transcript_cluster_id, logFC, adj.P.Val, gene_symbol)
toptable_NCTM_final  <- toptable_NCTM_rma  %>% dplyr::inner_join(microarray_annot_unique, by = "transcript_cluster_id") %>% dplyr::select(transcript_cluster_id, logFC, adj.P.Val, gene_symbol)
toptable_CSCTM_final <- toptable_CSCTM_rma %>% dplyr::inner_join(microarray_annot_unique, by = "transcript_cluster_id") %>% dplyr::select(transcript_cluster_id, logFC, adj.P.Val, gene_symbol)

write.csv(toptable_NCS_final,   file.path(PATH_TABLES, "toptable_NCS_rma.csv"),   row.names = FALSE)
write.csv(toptable_NCTM_final,  file.path(PATH_TABLES, "toptable_NCTM_rma.csv"),  row.names = FALSE)
write.csv(toptable_CSCTM_final, file.path(PATH_TABLES, "toptable_CSCTM_rma.csv"), row.names = FALSE)
write.csv(fit_results_rma, file.path(PATH_TABLES, "fit_results_rma.csv"), row.names = FALSE)

# ---- 6. Save objects for subsequent scripts ----
saveRDS(fit2_rma, file.path(PATH_PROCESSED, "fit2_rma.rds"))
saveRDS(
  list(NCS = toptable_NCS_rma, NCTM = toptable_NCTM_rma, CSCTM = toptable_CSCTM_rma),
  file.path(PATH_PROCESSED, "toptables_rma.rds")
)
saveRDS(
  list(matrix = exprs_matrix_rma, genes_info = genes_info_rma, sample_cols = sample_cols_rma),
  file.path(PATH_PROCESSED, "exprs_matrix_rma.rds")
)

cat("Top table dimensions -> NCS:", dim(toptable_NCS_rma),
    " | NCTM:", dim(toptable_NCTM_rma),
    " | CSCTM:", dim(toptable_CSCTM_rma), "\n")

