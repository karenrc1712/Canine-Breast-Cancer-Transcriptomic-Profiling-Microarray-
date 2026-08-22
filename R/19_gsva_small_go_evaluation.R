# =============================================================================
# Script:      19_gsva_small_go_evaluation.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Fernández
# Description: GO terms with few genes after data intersection (< 10) can
#              produce GSVA scores dominated by one or two highly variable
#              genes, making the score unreliable as a pathway-level estimate.
#              This script evaluates "gene dominance" for all small GO terms
#              in the upregulated set by:
#                - Computing the Pearson correlation of each gene's expression
#                  with the GO term's GSVA score, within each group
#                - Flagging terms where one gene contributes > 70% of the
#                  signal (max |r| > 0.7 or max/mean ratio > 2)
#              Results are exported as CSVs for manual review; problematic
#              terms should be interpreted with caution.
# Inputs:      data/processed/GSVA_analysis/go_gene_tables.rds
#              data/processed/GSVA_analysis/gsva_r.rds
#              data/processed/GSVA_analysis/gsva_targets_aligned.rds
#              data/processed/GSVA_analysis/expresion_final.rds
# Outputs:     results/tables/GSVA_analysis/summary_dominance_CS.csv
#              results/tables/GSVA_analysis/summary_dominance_CTM.csv
# Dependencies: 00_setup.R, 19_gsva_pathway_tables.R
# =============================================================================

source(here::here("R", "00_setup.R"))

go_gene_tables <- readRDS(file.path(PATH_GSVA_PROC, "go_gene_tables.rds"))
gsva_r         <- readRDS(file.path(PATH_GSVA_PROC, "gsva_r.rds"))
targets_gsva   <- readRDS(file.path(PATH_GSVA_PROC, "gsva_targets_aligned.rds"))
expresion_final<- readRDS(file.path(PATH_GSVA_PROC, "expresion_final.rds"))

go_table_up <- go_gene_tables$up$go_table
rownames(go_table_up) <- go_table_up$GO_ID

# ---- 1. Identify small GO terms (< 10 genes after intersection) ----
go_small_ids <- go_table_up$GO_ID[go_table_up$N_genes < 10]
cat("Small GO terms (< 10 genes in data):", length(go_small_ids), "\n")

grupos_named <- targets_gsva$Grupos
names(grupos_named) <- colnames(gsva_r)

# ---- 2. Gene dominance evaluation function ----
# Returns: per-gene Pearson r vs GSVA score, dominant gene, and a problem flag
evaluate_gene_dominance <- function(go_id, group_label) {
  genes_str <- go_table_up[go_id, "Genes"]
  genes     <- intersect(unlist(strsplit(genes_str, ", ")),
                         rownames(expresion_final))
  
  if (length(genes) < 2) return(NULL)
  
  samples <- names(grupos_named)[grupos_named == group_label]
  samples <- intersect(samples, colnames(expresion_final))
  
  if (length(samples) < 3) return(NULL)
  
  expr_mat <- as.matrix(expresion_final[genes, samples, drop = FALSE])
  gsva_vec <- as.numeric(gsva_r[go_id, samples])
  
  correlations <- apply(expr_mat, 1, function(x)
    cor(x, gsva_vec, use = "complete.obs", method = "pearson"))
  
  max_corr  <- max(abs(correlations))
  mean_corr <- mean(abs(correlations))
  
  list(
    go_id        = go_id,
    group        = group_label,
    n_samples    = length(samples),
    genes        = genes,
    correlations = correlations,
    dominant_gene = names(which.max(abs(correlations))),
    max_corr     = max_corr,
    mean_corr    = mean_corr,
    dominance_ratio = max_corr / mean_corr,
    flag_problem = max_corr > 0.7 | (max_corr / mean_corr) > 2
  )
}

# ---- 3. Run evaluation for CS and CTM ----
cat("Evaluating gene dominance across", length(go_small_ids), "small GO terms...\n")

results_CS  <- Filter(Negate(is.null), lapply(go_small_ids, evaluate_gene_dominance, "CS"))
results_CTM <- Filter(Negate(is.null), lapply(go_small_ids, evaluate_gene_dominance, "CTM"))

names(results_CS)  <- sapply(results_CS,  `[[`, "go_id")
names(results_CTM) <- sapply(results_CTM, `[[`, "go_id")

# ---- 4. Build summary tables ----
build_summary <- function(results) {
  do.call(rbind, lapply(results, function(r) {
    data.frame(
      GO_ID           = r$go_id,
      N_Genes         = length(r$genes),
      N_Samples       = r$n_samples,
      Max_Correlation  = round(r$max_corr, 3),
      Mean_Correlation = round(r$mean_corr, 3),
      Dominance_Ratio  = round(r$dominance_ratio, 2),
      Dominant_Gene    = r$dominant_gene,
      Flag_Problem     = ifelse(r$flag_problem, "YES", "NO"),
      stringsAsFactors = FALSE
    )
  }))
}

summary_CS  <- build_summary(results_CS)  [order(build_summary(results_CS)$Max_Correlation,  decreasing = TRUE), ]
summary_CTM <- build_summary(results_CTM) [order(build_summary(results_CTM)$Max_Correlation, decreasing = TRUE), ]

cat("Flagged problematic terms (CS):  ", sum(summary_CS$Flag_Problem  == "YES"), "\n")
cat("Flagged problematic terms (CTM): ", sum(summary_CTM$Flag_Problem == "YES"), "\n")

# ---- 5. Export ----
write.csv(summary_CS,  file.path(PATH_GSVA_TAB, "summary_dominance_CS.csv"),  row.names = FALSE)
write.csv(summary_CTM, file.path(PATH_GSVA_TAB, "summary_dominance_CTM.csv"), row.names = FALSE)

cat("\nGO terms flagged as problematic should be interpreted cautiously.\n")
cat("Consider excluding them from final figures or reporting N_genes\n")
cat("alongside the GSVA score in any publication.\n")
