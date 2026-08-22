# =============================================================================
# Script:      10_export_gene_lists_waterfall.R
# Project:     Transcriptomic profile in canine breast cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Fernández
# Description: Prepares ordered logFC data (for plotting a waterfall plot in
#              Prism) and exports up/down gene lists for each contrast as .txt
#              files, ready to upload to Enrichr (https://maayanlab.cloud/Enrichr/)
#              and obtain Gene Ontology terms (input for the next script).
# Inputs:      data/processed/degs_rma.rds
# Outputs:     results/tables/lf_{NCS,NCTM}_rma.csv
#              results/tables/genes_{up,down}_{NCS,NCTM}_rma.txt
# Dependencies: 00_setup.R, 07_extract_degs.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(dplyr)

degs <- readRDS(file.path(PATH_PROCESSED, "degs_rma.rds"))

# ---- 1. Data ordered by logFC for waterfall plot ----
lf_NCS_rma <- data.frame(
  gene_symbol = degs$NCS$total$gene_symbol,
  logFC = degs$NCS$total$logFC
) %>% arrange(desc(logFC))

lf_NCTM_rma <- data.frame(
  gene_symbol = degs$NCTM$total$gene_symbol,
  logFC = degs$NCTM$total$logFC
) %>% arrange(desc(logFC))

write.csv(lf_NCS_rma,  file.path(PATH_TABLES, "lf_NCS_rma.csv"),  row.names = FALSE)
write.csv(lf_NCTM_rma, file.path(PATH_TABLES, "lf_NCTM_rma.csv"), row.names = FALSE)

# ---- 2. Gene lists for Gene Ontology (Enrichr) ----
writeLines(degs$NCTM$down$gene_symbol, file.path(PATH_TABLES, "genes_down_NCTM_rma.txt"))
writeLines(degs$NCS$down$gene_symbol,  file.path(PATH_TABLES, "genes_down_NCS_rma.txt"))
writeLines(degs$NCTM$up$gene_symbol,   file.path(PATH_TABLES, "genes_up_NCTM_rma.txt"))
writeLines(degs$NCS$up$gene_symbol,    file.path(PATH_TABLES, "genes_up_NCS_rma.txt"))

cat("Gene lists exported to results/tables/. Upload them to Enrichr and\n")
cat("save the resulting .txt files in data/raw/enrichr/ for script 12.\n")

