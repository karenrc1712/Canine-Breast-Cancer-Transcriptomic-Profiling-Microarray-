# =============================================================================
# Script:      02_normalization_rma.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Hernandez
# Description: Applies RMA normalization (background correction + quantile
#              normalization + summarization at transcript cluster level) and
#              generates the post-normalization boxplot to verify that the
#              distributions across samples are aligned.
# Inputs:      data/processed/data_raw.rds
# Outputs:     results/tables/rma_normalized.txt
#              data/processed/exprs_rma.rds   (normalized matrix, log2)
# Dependencies: 00_setup.R, 01_qc_raw_data.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(oligo)
library(ggplot2)
library(reshape2)
library(pd.cangene.1.0.st)

# ---- 1. Load raw data ----
data.nc <- readRDS(file.path(PATH_PROCESSED, "data_raw.rds"))

# ---- 2. RMA normalization (Gene ST arrays, "core" level) ----
data.rma  <- oligo::rma(data.nc, target = "core")
exprs.rma <- exprs(data.rma)  # Normalized and log2-transformed matrix

dim(exprs.rma)
head(exprs.rma)

# ---- 3. Export normalized matrix ----
write.table(
  exprs.rma,
  file.path(PATH_TABLES, "rma_normalized.txt"),
  sep = "\t", quote = FALSE, col.names = NA
)

# ---- 4. Post-RMA boxplot ----
df_rma        <- reshape2::melt(exprs.rma)
df_rma$Var2   <- as.factor(df_rma$Var2)

boxplot_rma <- ggplot(df_rma, aes(x = Var2, y = value, fill = Var2)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +
  labs(
    title = "Boxplot after RMA normalization",
    x = "Samples",
    y = "Intensity (log2)"
  )

print(boxplot_rma)

# ---- 5. Save normalized matrix for subsequent scripts ----
saveRDS(exprs.rma, file.path(PATH_PROCESSED, "exprs_rma.rds"))
getwd()
cat("Saved: data/processed/exprs_rma.rds\n")

