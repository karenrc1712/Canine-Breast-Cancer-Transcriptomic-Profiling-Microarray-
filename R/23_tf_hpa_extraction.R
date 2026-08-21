# =============================================================================
# Script:      23_tf_hpa_extraction.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Hernandez
# Description: Downloads the Human Protein Atlas (HPA) full dataset and
#              extracts the transcription factor (TF) reference gene list:
#              genes with "Transcription factor" in the Protein.class field.
#              The download date is recorded for reproducibility, since HPA
#              releases are updated periodically.
# Inputs:      https://www.proteinatlas.org/download/proteinatlas.tsv.zip
# Outputs:     data/processed/TFs/TFs.rds
#              data/processed/TFs/proteinatlas.tsv   (raw database, large file)
#              results/tables/TFs/hpa_download_metadata.txt
# Dependencies: 00_setup.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(dplyr)

# ---- Define TF-specific paths ----
PATH_TF_PROC <- file.path(PATH_PROCESSED, "TFs")
PATH_TF_TAB  <- file.path(PATH_TABLES,    "TFs")
PATH_TF_FIG  <- file.path(PATH_FIGURES,   "TFs")
dir.create(PATH_TF_PROC, showWarnings = FALSE, recursive = TRUE)
dir.create(PATH_TF_TAB,  showWarnings = FALSE, recursive = TRUE)
dir.create(PATH_TF_FIG,  showWarnings = FALSE, recursive = TRUE)

# ---- 1. Download HPA ----
hpa_zip  <- file.path(PATH_TF_PROC, "proteinatlas.tsv.zip")
hpa_tsv  <- file.path(PATH_TF_PROC, "proteinatlas.tsv")
dl_date  <- Sys.Date()

cat("Downloading Human Protein Atlas (this may take a few minutes)...\n")
download.file(
  url      = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip",
  destfile = hpa_zip,
  mode     = "wb"
)
unzip(hpa_zip, exdir = PATH_TF_PROC)
cat("Download complete:", format(dl_date), "\n")

# ---- 2. Load HPA ----
hpa <- read.delim(hpa_tsv, sep = "\t", stringsAsFactors = FALSE)
cat("HPA loaded:", nrow(hpa), "genes\n")
cat("Columns:", paste(names(hpa), collapse = ", "), "\n")

# ---- 3. Extract transcription factors ----
TFs <- hpa %>%
  filter(grepl("Transcription factor", Protein.class, ignore.case = TRUE)) %>%
  pull(Gene) %>%
  unique()
cat("Transcription factors found:", length(TFs), "\n")

# ---- 4. Save download metadata for reproducibility ----
meta <- c(
  paste("Download date:", format(dl_date)),
  paste("Source: https://www.proteinatlas.org/download/proteinatlas.tsv.zip"),
  paste("HPA genes in database:", nrow(hpa)),
  paste("Transcription factors extracted:", length(TFs)),
  paste("R version:", R.version$version.string)
)
writeLines(meta, file.path(PATH_TF_TAB, "hpa_download_metadata.txt"))
cat(paste(meta, collapse = "\n"), "\n")

# ---- 5. Save object ----
saveRDS(TFs, file.path(PATH_TF_PROC, "TFs.rds"))
cat("Saved: TFs.rds\n")

