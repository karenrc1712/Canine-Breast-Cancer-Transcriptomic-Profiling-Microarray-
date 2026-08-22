# =============================================================================
# Script:      14_gsva_go_genesets.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Fernández
# Description: Builds GO Biological Process (BP) gene sets from the human
#              annotation database (org.Hs.eg.db). Applies three filters:
#                1. Retains only current (non-obsolete) GO terms
#                2. Restricts gene set size to 10-300 genes (broad enough
#                   for GSVA, small enough to be interpretable)
#                3. Intersects each gene set with the genes actually present
#                   in the expression matrix (exp_1to1), keeping sets with
#                   at least 7 genes after intersection
#              The filtered gene sets are the direct input to GSVA.
# Inputs:      data/processed/GSVA_analysis/exp_1to1.rds
# Outputs:     data/processed/GSVA_analysis/genesbygo_filtered.rds
#              results/tables/GSVA_analysis/genes_per_GO_term.csv
# Dependencies: 00_setup.R, 15_gsva_ortholog_mapping.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)
library(dplyr)

exp_1to1 <- readRDS(file.path(PATH_GSVA_PROC, "exp_1to1.rds"))

# ---- 1. Retrieve GO BP annotations for all human genes ----
cat("Retrieving GO Biological Process annotations...\n")
goannot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = keys(org.Hs.eg.db, keytype = "ENTREZID"),
  columns = c("GO", "ONTOLOGY", "SYMBOL"),
  keytype = "ENTREZID"
)

# ---- 2. Filter: Biological Process + valid (non-obsolete) GO terms ----
valid_go_ids <- keys(GO.db)

goannot_bp <- goannot %>%
  filter(ONTOLOGY == "BP",
         !is.na(GO),
         GO %in% valid_go_ids) %>%
  distinct(GO, SYMBOL)

cat("GO BP terms after removing obsolete:", length(unique(goannot_bp$GO)), "\n")

# ---- 3. Build gene sets; filter by size (10-300 genes) ----
genesbygo_bp <- split(goannot_bp$SYMBOL, goannot_bp$GO)
genesbygo_bp <- genesbygo_bp[
  sapply(genesbygo_bp, length) >= 10 &
    sapply(genesbygo_bp, length) <= 300
]
cat("GO terms with 10-300 genes:", length(genesbygo_bp), "\n")

# ---- 4. Intersect each gene set with genes present in the data ----
genes_in_data      <- rownames(exp_1to1)
genesbygo_filtered <- lapply(genesbygo_bp, function(gs) intersect(gs, genes_in_data))

# Keep sets with at least 7 genes after intersection
genesbygo_filtered <- genesbygo_filtered[
  sapply(genesbygo_filtered, length) >= 7
]

cat("GO terms retained after intersection (>= 7 genes in data):",
    length(genesbygo_filtered), "\n")
cat("Gene set size summary:\n")
print(summary(sapply(genesbygo_filtered, length)))

# ---- 5. Export gene-per-GO table ----
df_genes_go <- stack(genesbygo_filtered)
colnames(df_genes_go) <- c("gene_symbol", "GO_ID")
write.csv(df_genes_go,
          file.path(PATH_GSVA_TAB, "genes_per_GO_term.csv"),
          row.names = FALSE)

# ---- 6. Save gene sets ----
saveRDS(genesbygo_filtered,
        file.path(PATH_GSVA_PROC, "genesbygo_filtered.rds"))
cat("Saved: genesbygo_filtered.rds\n")

