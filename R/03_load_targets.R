# =============================================================================
# Script:      03_load_targets.R

# Project:     Transcriptomic profiling in canine mammary cancer (microarray)

# Author:      Karen Rodriguez

# Advisor:     Geysson Javier Fernández

# Description: Loads the metadata table (targets: which sample belongs to
#              which experimental group: N = normal, CS = simple carcinoma,
#              CTM = tubulomedullary carcinoma) and matches it to the actual
#              order of the .CEL files. This object (targets1) is used by
#              almost all subsequent scripts, so it is saved separately.

# Inputs:      data/raw/targets1.xlsx - Sheet1.tsv
#              data/raw/*.CEL (only to read file names, not intensity values)

# Outputs:     data/processed/targets1.rds

# Dependencies: 00_setup.R

# =============================================================================

source(here::here("R", "00_setup.R"))
library(limma)
library(oligo)

# ---- 1. Load metadata table ----
targets <- readTargets(file.path(PATH_RAW, "targets1.xlsx - Sheet1.tsv"))

# ---- 2. Match with the actual order of .CEL files ----
cel_files2 <- list.celfiles(path = PATH_RAW, full.names = FALSE)
cel_ids2   <- as.integer(sub("_.*", "", cel_files2))

targets1 <- targets[targets$Renomear %in% cel_ids2, ]
targets1 <- targets1[match(cel_ids2, targets1$Renomear), ]

head(targets1)
cat("Groups found:", paste(unique(targets1$Grupos), collapse = ", "), "\n")

# ---- 3. Save for subsequent scripts ----
saveRDS(targets1, file.path(PATH_PROCESSED, "targets1.rds"))
cat("Saved: data/processed/targets1.rds\n")

