# =============================================================================
# Script:      16_gsva_differential_enrichment.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Fernández
# Description: Applies limma to the GSVA enrichment score matrix to identify
#              GO Biological Process terms that are significantly and
#              differentially enriched between experimental groups.
#              Three contrasts are tested:
#                CS_vs_N   : simple carcinoma vs normal
#                CTM_vs_N  : tubulomedullary carcinoma vs normal
#                CS_vs_CTM : comparison between the two carcinoma subtypes
#              Multiple testing is corrected with Benjamini-Hochberg (FDR).
# Inputs:      data/processed/GSVA_analysis/gsva_r.rds
#              data/processed/targets1.rds
# Outputs:     results/tables/GSVA_analysis/res_{CS_vs_N, CTM_vs_N,
#                CS_vs_CTM}.csv
#              data/processed/GSVA_analysis/gsva_topTables.rds
#              data/processed/GSVA_analysis/gsva_targets_aligned.rds
# Dependencies: 00_setup.R, 03_load_targets.R, 17_gsva_run.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(limma)

gsva_r   <- readRDS(file.path(PATH_GSVA_PROC, "gsva_r.rds"))
targets1 <- readRDS(file.path(PATH_PROCESSED,  "targets1.rds"))

# ---- 1. Align GSVA column names with targets ----
sample_ids              <- as.numeric(colnames(gsva_r))
targets1$Renomear       <- as.character(as.numeric(targets1$Renomear))
targets_gsva            <- targets1[targets1$Renomear %in% sample_ids, ]
targets_gsva            <- targets_gsva[match(sample_ids, targets_gsva$Renomear), ]

stopifnot("Sample-target alignment failed" =
            all(targets_gsva$Renomear == sample_ids))

group <- factor(targets_gsva$Grupos)
cat("Group distribution:\n")
print(table(group))

# ---- 2. Design matrix and contrasts ----
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
print(design)

contr.matrix <- makeContrasts(
  CS_vs_N   = CS  - N,
  CTM_vs_N  = CTM - N,
  CS_vs_CTM = CS  - CTM,
  levels = design
)
print(contr.matrix)

# ---- 3. Linear model fitting ----
fit  <- lmFit(gsva_r, design)
fit2 <- contrasts.fit(fit, contr.matrix)
fit2 <- eBayes(fit2)

# ---- 4. Extract results for each contrast ----
res_CS_vs_N   <- topTable(fit2, coef = "CS_vs_N",   number = Inf, adjust.method = "BH")
res_CTM_vs_N  <- topTable(fit2, coef = "CTM_vs_N",  number = Inf, adjust.method = "BH")
res_CS_vs_CTM <- topTable(fit2, coef = "CS_vs_CTM", number = Inf, adjust.method = "BH")

cat("\nSignificant GO terms (adj.P.Val < 0.05):\n")
cat("  CS vs N:   ",
    sum(res_CS_vs_N$adj.P.Val   < 0.05, na.rm = TRUE), "\n")
cat("  CTM vs N:  ",
    sum(res_CTM_vs_N$adj.P.Val  < 0.05, na.rm = TRUE), "\n")
cat("  CS vs CTM: ",
    sum(res_CS_vs_CTM$adj.P.Val < 0.05, na.rm = TRUE), "\n")

# ---- 5. Export CSVs ----
write.csv(res_CS_vs_N,   file.path(PATH_GSVA_TAB, "res_CS_vs_N.csv"),   row.names = TRUE)
write.csv(res_CTM_vs_N,  file.path(PATH_GSVA_TAB, "res_CTM_vs_N.csv"),  row.names = TRUE)
write.csv(res_CS_vs_CTM, file.path(PATH_GSVA_TAB, "res_CS_vs_CTM.csv"), row.names = TRUE)

# ---- 6. Save R objects ----
saveRDS(
  list(CS_vs_N = res_CS_vs_N, CTM_vs_N = res_CTM_vs_N, CS_vs_CTM = res_CS_vs_CTM),
  file.path(PATH_GSVA_PROC, "gsva_topTables.rds")
)
saveRDS(targets_gsva, file.path(PATH_GSVA_PROC, "gsva_targets_aligned.rds"))
cat("Saved: gsva_topTables.rds | gsva_targets_aligned.rds\n")

