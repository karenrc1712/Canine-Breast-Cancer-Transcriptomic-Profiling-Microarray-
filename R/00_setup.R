# =============================================================================
# Script:      00_setup.R
# Project:     Canine Breast Cancer Transcriptomic Profiling (Microarray)
# Author:      Karen Rodriguez
# Supervisor:  Geysson Javier Fernández
#
# Description:
# Installs and loads all packages required by the analysis pipeline.
# It also defines the standard project directories.
# This script should be run ONCE at the beginning of every R session.
#
# Inputs:
# None
#
# Outputs:
# None (creates project folders and loads required packages)
#
# Dependencies:
# None (base script)
# =============================================================================

# ---- 1. Bioconductor packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c("limma", "oligo", "pd.cangene.1.0.st")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# ---- 2. CRAN packages ----
cran_pkgs <- c(
  "here", "reshape2", "preprocessCore", "ggplot2", "ggrepel",
  "dplyr", "tidyr", "matrixStats", "pheatmap", "RColorBrewer",
  "DT", "openxlsx", "readr", "ggalluvial", "showtext", "sysfonts")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# ---- 3. Load libraries used across the pipeline ----
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- 4. Standard project paths ----
# here() automatically locates the project root (where the .Rproj is),
# thus avoiding absolute paths like "C:/Users/Karen/...".

PATH_RAW        <- here::here("data", "raw")          # .CEL, annotation, targets
PATH_PROCESSED  <- here::here("data", "processed")     # intermediate .rds objects
PATH_TABLES     <- here::here("results", "tables")      # final .csv / .xlsx
PATH_FIGURES    <- here::here("results", "figures")     # .pdf figures
PATH_GSVA_PROC <- here::here("data", "processed", "GSVA")
PATH_GSVA_TAB <- here::here("results", "tables", 'GSVA')
PATH_GSVA_FIG <- here::here("results", "figures", 'GSVA')
PATH_BOXPLOTS_TAB <- here::here("results", "tables", "boxplots")

dir.create(PATH_BOXPLOTS_TAB, showWarnings = FALSE, recursive = TRUE)
dir.create(PATH_PROCESSED, showWarnings = FALSE, recursive = TRUE)
dir.create(PATH_TABLES,    showWarnings = FALSE, recursive = TRUE)
dir.create(PATH_FIGURES,   showWarnings = FALSE, recursive = TRUE)
dir.create(PATH_GSVA_PROC,   showWarnings = FALSE, recursive = TRUE)
dir.create(PATH_GSVA_FIG,   showWarnings = FALSE, recursive = TRUE)
dir.create(PATH_GSVA_TAB,   showWarnings = FALSE, recursive = TRUE)
cat("Setup complete. Project paths:\n")
cat("  Raw:        ", PATH_RAW, "\n")
cat("  Processed:  ", PATH_PROCESSED, "\n")
cat("  Tables:     ", PATH_TABLES, "\n")
cat("  Figures:    ", PATH_FIGURES, "\n")
cat("  Tables_Boxplots:    ", PATH_BOXPLOTS_TAB, "\n")
