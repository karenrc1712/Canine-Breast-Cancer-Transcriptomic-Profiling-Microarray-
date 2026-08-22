# =============================================================================
# Script:      01_qc_raw_data.R
#
# Project:
# Canine Breast Cancer Transcriptomic Profiling (Microarray)
#
# Author:
# Karen Rodriguez
#
# Supervisor:
# Geysson Javier Fernández
#
# Description:
# Loads raw .CEL files and performs quality control BEFORE normalization.
# This includes chip pseudo-images, raw intensity boxplots, density plots,
# and Relative Log Expression (RLE) using a Probe Level Model.
#
# Input:
# data/raw/*.CEL
#
# Output:
# results/figures/QC_pseudoimages_rank.pdf
# data/processed/data_raw.rds
#
# Dependencies:
# 00_setup.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(oligo)

# ---- 1. Load .CEL files ----
cel_files <- list.celfiles(path = PATH_RAW, full.names = TRUE)
data.nc   <- read.celfiles(cel_files)

# Clean sample names (remove everything after the first "_")
sampleNames(data.nc) <- sub("_.*", "", sampleNames(data.nc))

cat("Samples loaded:", length(sampleNames(data.nc)), "\n")
head(data.nc)

# ---- 2. Pseudoimages (detection of physical chip defects) ----
pdf(file.path(PATH_FIGURES, "QC_pseudoimages_rank.pdf"), width = 12, height = 10)
par(mfrow = c(5, 5), mar = c(1, 1, 2, 1))
for (i in seq_len(ncol(data.nc))) {
  image(data.nc[, i], transfo = rank, main = sampleNames(data.nc)[i])
}
dev.off()

# ---- 3. Raw intensity boxplot (log2) ----
exprs_raw <- log2(exprs(data.nc))

boxplot(
  exprs_raw,
  outline  = FALSE,
  las      = 2,
  cex.axis = 0.7,
  main     = "Raw log2 probe intensities"
)

# ---- 4. Density plots per sample ----
plot(density(exprs_raw[, 1]), main = "Raw intensity density")
for (i in 2:ncol(exprs_raw)) {
  lines(density(exprs_raw[, i]))
}

# ---- 5. RLE (Relative Log Expression) via Probe Level Model ----
plm <- fitProbeLevelModel(data.nc)
par(mar = c(10, 4, 4, 2))
RLE(
  plm,
  main     = "RLE plot",
  outline  = FALSE,
  las      = 2,
  cex.axis = 0.7
)

# NUSE is available if sample-wise precision assessment is required:
# nuse_values <- NUSE(plm)
# boxplot(nuse_values, main = "NUSE - Raw Data", las = 2)

# ---- 6. Save raw object for the next script ----
saveRDS(data.nc, file.path(PATH_PROCESSED, "data_raw.rds"))
cat("Saved: data/processed/data_raw.rds\n")
