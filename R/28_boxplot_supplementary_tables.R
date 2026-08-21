# =============================================================================
# Script:      28_boxplot_supplementary_tables.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Hernandez
# Description: Builds the input tables for the supplementary boxplot figure.
#              For each direction (UP / DOWN) independently:
#                1. Takes the GO terms + Hallmark assignment from
#                   GO_selected_up/down_gsva_4.
#                2. Pulls the gene list per GO term from tabla_stats_up/down.
#                3. Intersects those genes with the DEG gene universe for
#                   that SAME direction (DEGs_up/down_NCS_rma and
#                   DEGs_up/down_NCTM_rma) -- i.e. UP genes are only checked
#                   against UP DEG files, DOWN genes only against DOWN DEG
#                   files.
#                4. For each surviving (GO term, gene) pair, attaches the
#                   gene's own logFC from the DEG files as "Combined score"
#                   for CS and CTM respectively (NA if that gene isn't DE in
#                   that particular contrast).
#              Output is long-format: one row per GO term x gene, ready to
#              feed a boxplot of logFC distribution per GO term / Hallmark.
# Inputs:      PATH_GSVA_TAB/GO_selected_up_gsva_4.csv
#              PATH_GSVA_TAB/GO_selected_down_gsva_4.csv
#              PATH_GSVA_TAB/tabla_stats_up.csv
#              PATH_GSVA_TAB/tabla_stats_down.csv
#              PATH_TABLES/DEGs_up_NCS_rma.csv
#              PATH_TABLES/DEGs_up_NCTM_rma.csv
#              PATH_TABLES/DEGs_down_NCS_rma.csv
#              PATH_TABLES/DEGs_down_NCTM_rma.csv
# Outputs:     results/tables/boxplots/boxplot_input_up.csv
#              results/tables/boxplots/boxplot_input_down.csv
#              data/processed/GSVA_analysis/boxplot_input.rds
# Dependencies: 00_setup.R, 19_gsva_pathway_tables.R, 21_alluvial_deg_direction.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

source(here::here("R", "00_setup.R"))
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

# ---- 0. Load expression matrix + sample groups, create output folder ----
expresion_final <- readRDS(file.path(PATH_GSVA_PROC, "expresion_final.rds"))
targets1         <- readRDS(file.path(PATH_PROCESSED, "targets1.rds"))

# Align column names between expression matrix and targets
colnames(expresion_final) <- as.character(colnames(expresion_final))
targets1$Renomear         <- sprintf("%02d", as.numeric(targets1$Renomear))
targets1                  <- targets1 %>% filter(Renomear %in% colnames(expresion_final))

muestras_N   <- targets1 %>% filter(Grupos == "N")   %>% pull(Renomear)
muestras_CS  <- targets1 %>% filter(Grupos == "CS")  %>% pull(Renomear)
muestras_CTM <- targets1 %>% filter(Grupos == "CTM") %>% pull(Renomear)

PATH_BOXPLOTS_TAB <- file.path(PATH_TABLES, "boxplots")
dir.create(PATH_BOXPLOTS_TAB, showWarnings = FALSE, recursive = TRUE)
# ---- 1. Read source tables ----
go_selected_up   <- read_excel(file.path(PATH_GSVA_TAB, "GO_selected_up_gsva_4.xlsx")) %>%
  as.data.frame()
go_selected_down <- read_excel(file.path(PATH_GSVA_TAB, "GO_selected_down_gsva_4.xlsx")) %>%
  as.data.frame()

tabla_stats_up   <- read.csv(file.path(PATH_GSVA_TAB, "go_gene_table_up.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
tabla_stats_down <- read.csv(file.path(PATH_GSVA_TAB, "go_gene_table_down.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)

# ---- 2. Helper: read a DEGs_* file, standardize to Gene / logFC ----
read_deg_file <- function(path, contrast_label) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot("gene_symbol" %in% colnames(df))
  
  df <- df %>%
    dplyr::rename(Gene = gene_symbol) %>%
    dplyr::select(Gene, logFC)
  
  colnames(df)[colnames(df) == "logFC"] <- paste0("logFC_", contrast_label)
  df
}

deg_up_cs    <- read_deg_file(file.path(PATH_TABLES, "DEGs_up_NCS_rma.csv"),    "CS")
deg_up_ctm   <- read_deg_file(file.path(PATH_TABLES, "DEGs_up_NCTM_rma.csv"),   "CTM")
deg_down_cs  <- read_deg_file(file.path(PATH_TABLES, "DEGs_down_NCS_rma.csv"),  "CS")
deg_down_ctm <- read_deg_file(file.path(PATH_TABLES, "DEGs_down_NCTM_rma.csv"), "CTM")

# ---- 3. Helper: build the long GO-term x gene boxplot input table ----
# go_selected  : GO_ID, GO_Term, HALLMARK, logFC_CS, adj.P.Val_CS, logFC_CTM, adj.P.Val_CTM
# tabla_stats  : GO_ID, GO_Term, N_genes, Genes, FC_CS, adjpval_CS, FC_CTM, adjpval_CTM
# deg_cs/deg_ctm : Gene, logFC_CS / logFC_CTM  (DEG-level, one row per gene)
build_boxplot_table <- function(go_selected, tabla_stats, deg_cs, deg_ctm) {
  
  # ---- 3a. Hallmark + gene list per GO term (join by GO_ID) ----
  go_genes <- tabla_stats %>%
    dplyr::select(GO_ID, GO_Term, Genes) %>%
    left_join(
      go_selected %>%
        dplyr::select(GO_ID, HALLMARK),
      by = "GO_ID"
    )
  
  n_missing <- sum(is.na(go_genes$Genes))
  if (n_missing > 0) {
    warning(n_missing, " GO term(s) in GO_selected have no match in tabla_stats ",
            "(check GO_ID formatting -- e.g. stray whitespace like 'GO:0048304 ').")
  }
  
  # ---- 3b. Explode to one row per (GO term, gene) ----
  go_gene_long <- go_genes %>%
    filter(!is.na(Genes)) %>%
    mutate(Genes = str_trim(Genes)) %>%
    separate_rows(Genes, sep = ",\\s*") %>%
    dplyr::rename(Gene = Genes) %>%
    mutate(Gene = str_trim(Gene)) %>%
    filter(Gene != "") %>%
    distinct(GO_ID, GO_Term, HALLMARK, Gene)
  
  # ---- 3c. Intersect with the DEG universe for this direction ----
  deg_gene_universe <- union(deg_cs$Gene, deg_ctm$Gene)
  
  n_before <- n_distinct(go_gene_long$Gene)
  go_gene_long <- go_gene_long %>%
    filter(Gene %in% deg_gene_universe)
  n_after <- n_distinct(go_gene_long$Gene)
  
  cat("  Genes before DEG intersection:", n_before, "\n")
  cat("  Genes after DEG intersection: ", n_after, "\n")
  
  # ---- 3d. Attach the gene's own logFC ("Combined score") per contrast ----
  go_gene_long %>%
    left_join(deg_cs,  by = "Gene") %>%
    left_join(deg_ctm, by = "Gene") %>%
    dplyr::rename(
      logfc_CS  = starts_with("logFC_CS"),
      logfc_CTM = starts_with("logFC_CTM")
    ) %>%
    arrange(HALLMARK, GO_Term, Gene)
}

# ============================================================
# UPREGULATED
# ============================================================
cat("--- Building boxplot input table: UPREGULATED ---\n")
boxplot_up <- build_boxplot_table(go_selected_up, tabla_stats_up, deg_up_cs, deg_up_ctm)

# ============================================================
# DOWNREGULATED
# ============================================================
cat("--- Building boxplot input table: DOWNREGULATED ---\n")
boxplot_down <- build_boxplot_table(go_selected_down, tabla_stats_down, deg_down_cs, deg_down_ctm)

# ---- 4. Export ----
write.csv(boxplot_up,   file.path(PATH_BOXPLOTS_TAB, "boxplot_input_up.csv"),   row.names = FALSE)
write.csv(boxplot_down, file.path(PATH_BOXPLOTS_TAB, "boxplot_input_down.csv"), row.names = FALSE)

saveRDS(
  list(up =    boxplot_up, down = boxplot_down),
  file.path(PATH_GSVA_PROC, "boxplot_input.rds")
)

cat("Saved: boxplot_input_up.csv, boxplot_input_down.csv, boxplot_input.rds\n")


####GENES FOR SELECTED FOR BOXPLOT
############################
### GENES FOR BOXPLOTS
############################
genes_metabolism <- c(
  # Amino acid & Glutathione
  "GCLM", "SLC7A11", "SLC1A1", "SLC1A2", "SLC16A10", "AADAT",
  
  # Glucose & Lipid metabolism
  "PFKP", "SQLE", "CYP1B1", "CYB561",
  
  # Ion & Proton transport
  "ATP1B1", "ATP6V1A", "ATP6V1B1", 
  
  # Glycan & Receptor transport
  "ST6GAL2", "LRP2"
)

genes_cell_death <- c(
  # Ferroptosis & Redox
  "SLC7A11", "GPX2", "ALOX5AP", "CP",
  
  # Apoptosis signaling
  "THBS1", "CLU", "S100A9", "NR4A2",
  
  # Stress & DNA repair
  "PDK3", "RAD18", "MAL",
  
  # Survival & Receptor signaling
  "BMPR1B", "SPP1", "ADGRF1"
)

genes_invasion <- c(
  # ECM remodeling
  "MMP9", "MMP13", "FN1", "TNC", "LTBP2",
  
  # Cell adhesion
  "FERMT1", "ITGA2", "CDH3", "CDH2",
  
  # WNT-mediated migration
  "WNT5A", "SFRP2", "SFRP1",
  
  # Ephrin signaling
  "EFNA5", "NGEF",
  
  # Vascular migration
  "IGFBP5", "CCL28"
)

genes_microenvironment <- c(
  # Immune & Inflammation
  "APOE", "LCN2", "LTF", "SPP1",
  
  # Stromal crosstalk & Hypoxia
  "NRG1", "STC1",
  
  # Cell adhesion & Checkpoint
  "IGSF11",
  
  # Niche & Homeostasis
  "SLC26A4", "RDH10", "DLX5"
)

genes_proliferation <- c(
  # Growth factors & Ligands
  "EREG", "NRG1", "TGFB2", "MMP7",
  
  # Transcription & Cell cycle
  "FOS", "ID4", "SIX1", "ASPM", "DLX5",
  
  # WNT signaling
  "WNT5A", "SFRP1", "SFRP2",
  
  # Adhesion & Membrane receptors
  "FN1", "CSPG5", "SERPINB5", "GABRA4"
)
genes_all <- unique(c(
  genes_metabolism,
  genes_cell_death,
  genes_invasion,
  genes_microenvironment,
  genes_proliferation
))

# Remove empty and NA
genes_all <- genes_all[!is.na(genes_all) & genes_all != ""]

genes_present <- genes_all[genes_all %in% boxplot_up$Gene]
genes_missing <- genes_all[!genes_all %in% boxplot_up$Gene]

###########################################################
# 1. CALCULATION OF RELATIVE LOG2FC (NORMALIZED BY N GROUP)
###########################################################
# 1. log2FC MATRIX
###########################################################

# Calculate the average expression for the N group
mean_log_N_down <- rowMeans(expresion_final[, muestras_N], na.rm = TRUE)

# Subtract that mean from the entire matrix
# N is centered around 0
matriz_log2FC_down <- sweep(
  expresion_final,
  1,
  mean_log_N_down,
  "-"
)

# Order columns: N, CS, CTM
matriz_log2FC_down <- matriz_log2FC_down[
  ,
  c(muestras_N, muestras_CS, muestras_CTM)
]

# Convert to dataframe and put genes as a column
matriz_log2FC_down <- data.frame(
  gene_symbol = rownames(matriz_log2FC_down),
  matriz_log2FC_down,
  row.names = NULL,
  check.names = FALSE
)


###########################################################
# 2. EXPORT FUNCTION FOR PRISM
###########################################################

crear_matriz_prism <- function(lista_genes, nombre_archivo) {
  
  # Filter only the genes of interest
  matriz <- matriz_log2FC_down %>%
    dplyr::filter(gene_symbol %in% lista_genes) %>%
    dplyr::mutate(
      gene_symbol = factor(
        gene_symbol,
        levels = lista_genes
      )
    ) %>%
    dplyr::arrange(gene_symbol)
  
  # Build path within PATH_BOXPLOTS_TAB
  ruta_salida <- file.path(
    PATH_BOXPLOTS_TAB,
    nombre_archivo
  )
  
  # Export
  write.csv(
    matriz,
    ruta_salida,
    row.names = FALSE
  )
  
  cat(
    "File exported successfully:",
    ruta_salida,
    "\n"
  )
}
############EXPORT TO PRISM
############################
### EXPORT CSV FOR PRISM
############################

crear_matriz_prism(
  genes_proliferation,
  "boxplot_proliferative_signaling.csv"
)

crear_matriz_prism(
  genes_metabolism,
  "boxplot_cellular_energetics.csv"
)

crear_matriz_prism(
  genes_cell_death,
  "boxplot_cell_death_ros.csv"
)

crear_matriz_prism(
  genes_invasion,
  "boxplot_invasion_migration.csv"
)

crear_matriz_prism(
  genes_microenvironment,
  "boxplot_tumor_microenvironment.csv"
)

library(dplyr)

###########################################################
# 3. GENE LISTS BY HALLMARK (DOWNREGULATED)
###########################################################

genes_proliferative_down <- c(
  # Loss of Differentiation & Cell Fate Specification
  "PROX1", "MAB21L1", "TBX5", "ZFPM2",
  
  # Loss of Growth Inhibitory & Tumor Suppressor Signaling
  "TGFBR3", "HPGD", "BMP4",
  
  # Cell Surface Receptors & Signal Transduction Downregulation
  "S1PR1", "F2RL2", "GRIA1", "RCAN2",
  
  # Growth Factors & MicroRNA Regulation
  "FGF10", "MIR10A"
  
  
)

# 2. Deregulating Cellular Energetics
genes_energetics_down <- c(
  # Master Nuclear Regulators of Lipid Homeostasis
  "PPARG", "NR1H3", "NR5A2",
  
  # Fatty Acid Transport & Lipid Clearance
  "CD36", "ABCA8",
  
  # Triglyceride Synthesis & Lipid Droplet Dynamics
  "DGAT2", "CIDEC", "ACVR1C",
  
  # Adipokine & Endocrine Metabolic Homeostasis
  "ADIPOQ", "IGF1"
)

# 3. Activating Invasion and Metastasis
genes_invasion_down <- c(
  # Loss of ECM Integrity & Basement Membrane Architecture
  "LAMA2", "COL6A1", "VTN", "ECM1", "MMRN1",
  
  # Cell-Cell Junctions & Intercellular Adhesion Loss
  "JAM2", "VCAM1",
  
  # Loss of EMT Suppression & Motility Control
  "CHRDL1", "RASGRF2", "TGFBR3", "PDPN"
  
)
# 4. Inducing Angiogenesis
genes_angiogenesis_down <- c(
  # Loss of Endothelial Adherens Junctions & Vascular Barrier
  "CDH5", "PECAM1", "PTPRB", "RAMP2",
  
  # Endothelial Receptors & Tie2/VEGF Signaling Disruption
  "TEK", "FLT4", "CD34", "FZD4",
  
  # Pericyte Coverage & Vascular Smooth Muscle Dysfunction
  "PDGFRB", "ABCC9", "S1PR1",
  
  # Endothelial Motility & Vessel Remodeling Loss
  "ANGPT2", "RHOJ", "MIR126", "FYN", "RNF213", "MEIS1",
  
  # Loss of Endogenous Angiostatic Inhibitors & Brakes
  "DCN", "CXCL10"
)

# 5. Evading Immune Destruction
genes_immune_down <- c(
  # Impaired NK & CD8+ T-Cell Cytotoxic Function
  "KLRK1", "KLRD1", "SH2D1A", "CRTAM",
  
  # T-Cell Receptor (TCR) Signaling & Co-stimulation Loss
  "CD3D", "LCK", "PRKCQ", "CD2", "CD6", "CD69",
  
  # Loss of Lymphocyte Trafficking, Homing & Adhesion
  "CCL5", "CXCL12", "ITGAL", "ITGA4", "DOCK2",
  
  # Innate Immune Regulation & Purine Metabolism Loss
  "ADA", "COLEC12", "BCL6B",
  
  # Loss of Immunosuppressive Checkpoint (Enzymatic Control)
  "IDO1"
)


###########################################################
# 4. CSV
###########################################################

crear_matriz_prism(genes_proliferative_down, "Box_1_Proliferative_Signaling.csv")
crear_matriz_prism(genes_energetics_down,    "Box_2_Cellular_Energetics.csv")
crear_matriz_prism(genes_invasion_down,      "Box_3_Invasion_Metastasis.csv")
crear_matriz_prism(genes_angiogenesis_down,  "Box_4_Inducing_Angiogenesis.csv")
crear_matriz_prism(genes_immune_down,        "Box_5_Immune_Destruction.csv")

