# =============================================================================
# Script:      13_gsva_ortholog_mapping.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Fernández
# Description: Maps canine gene symbols (Canis lupus familiaris) to their
#              human (Homo sapiens) orthologs using Ensembl BioMart. Only
#              strict one-to-one orthologs are retained to ensure unambiguous
#              cross-species mapping. This is required to use human GO gene
#              sets (org.Hs.eg.db) in the subsequent GSVA analysis.
#
#              Mapping statistics are printed to the console so the researcher
#              can assess how many genes were retained after filtering.
# Inputs:      data/processed/GSVA_analysis/expresion_final.rds
# Outputs:     data/processed/GSVA_analysis/exp_1to1.rds
#              data/processed/GSVA_analysis/orthologs_1to1.rds
# Dependencies: 00_setup.R, 13_gsva_prepare_expression.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(biomaRt)
library(dplyr)

expresion_final <- readRDS(file.path(PATH_GSVA_PROC, "expresion_final.rds"))

# ---- 1. Connect to Ensembl BioMart ----
# Uses the useast mirror; change to "uswest" or "www" if connection fails
options(timeout = 200)
cat("Connecting to Ensembl BioMart (useast mirror)...\n")
ensembl <- useEnsembl(biomart = "genes", mirror = "useast")
dog     <- useDataset("clfamiliaris_gene_ensembl", mart = ensembl)
ensembl@version
# ---- 2. Retrieve human orthologs for all canine genes ----
cat("Querying orthologs for", nrow(expresion_final), "canine genes...\n")
orthologs <- getBM(
  attributes = c(
    "external_gene_name",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_type"
  ),
  filters = "external_gene_name",
  values  = rownames(expresion_final),
  mart    = dog
)

cat("Orthology types found:\n")
print(table(orthologs$hsapiens_homolog_orthology_type))

# ---- 3. Filter to strict one-to-one orthologs only ----
orthologs_1to1 <- orthologs %>%
  filter(
    hsapiens_homolog_orthology_type == "ortholog_one2one",
    hsapiens_homolog_associated_gene_name != ""
  )

cat("\nOne-to-one ortholog mapping:\n")
cat("  Canine genes queried:             ", nrow(expresion_final), "\n")
cat("  Unique canine genes with 1:1 hit: ",
    length(unique(orthologs_1to1$external_gene_name)), "\n")
cat("  Unique human gene symbols:        ",
    length(unique(orthologs_1to1$hsapiens_homolog_associated_gene_name)), "\n")

# ---- 4. Subset expression matrix to genes with a 1:1 ortholog ----
exp_1to1 <- expresion_final[
  rownames(expresion_final) %in% orthologs_1to1$external_gene_name, ]

# ---- 5. Rename rows: canine symbol → human symbol ----
map_vec <- setNames(
  orthologs_1to1$hsapiens_homolog_associated_gene_name,
  orthologs_1to1$external_gene_name
)
map_vec         <- map_vec[rownames(exp_1to1)]
rownames(exp_1to1) <- map_vec

# Remove any remaining duplicate human symbols (keeps first occurrence)
exp_1to1 <- exp_1to1[!duplicated(rownames(exp_1to1)), ]
exp_1to1 <- exp_1to1[!is.na(rownames(exp_1to1)), ]

cat("\nFinal matrix with human gene symbols:", dim(exp_1to1), "\n")

# ---- 6. Save ----
saveRDS(exp_1to1,       file.path(PATH_GSVA_PROC, "exp_1to1.rds"))
saveRDS(orthologs_1to1, file.path(PATH_GSVA_PROC, "orthologs_1to1.rds"))
cat("Saved: exp_1to1.rds | orthologs_1to1.rds\n")

