# =============================================================================
# Script:      24_tf_intersections_volcanos.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Hernandez
# Description: (1) Extracts upregulated and downregulated DEG gene lists from
#              the limma toptables (|logFC| >= 1, FDR < 0.05). (2) Intersects
#              DEG lists with the HPA-derived TF reference list.
#              (3) Produces volcano plots highlighting TFs.
# Inputs:      results/tables/toptable_{NCS,NCTM,CSCTM}_rma.csv
#              data/processed/TFs/TFs.rds
#
# Outputs:     data/processed/TFs/deg_lists.rds
#              data/processed/TFs/tf_intersections.rds
#              results/tables/TFs/tf_summary.csv
#              results/figures/TFs/volcano_TF_{CS,CTM}.pdf
# Dependencies: 00_setup.R, 07_extract_degs.R, 23_tf_hpa_extraction.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(dplyr)
library(ggplot2)

PATH_TF_PROC <- file.path(PATH_PROCESSED, "TFs")
PATH_TF_TAB  <- file.path(PATH_TABLES,    "TFs")
PATH_TF_FIG  <- file.path(PATH_FIGURES,   "TFs")

TFs <- readRDS(file.path(PATH_TF_PROC, "TFs.rds"))

# ---- 1. Load toptables ----
read_top <- function(contrast)
  read.csv(file.path(PATH_TABLES, paste0("toptable_", contrast, "_rma.csv")),
           stringsAsFactors = FALSE)

top_NCS   <- read_top("NCS")
top_NCTM  <- read_top("NCTM")

# ---- 2. Extract DEG gene lists (shared thresholds) ----
LFC_THR  <- 1
PADJ_THR <- 0.05

get_degs <- function(toptable) {
  list(
    up   = toptable %>% filter(logFC  >  LFC_THR,  adj.P.Val < PADJ_THR) %>% pull(gene_symbol),
    down = toptable %>% filter(logFC  < -LFC_THR,  adj.P.Val < PADJ_THR) %>% pull(gene_symbol)
  )
}

degs <- list(NCS = get_degs(top_NCS), NCTM = get_degs(top_NCTM))

cat("DEG counts:\n")
for (ctr in names(degs))
  cat(" ", ctr, "- Up:", length(degs[[ctr]]$up),
      "| Down:", length(degs[[ctr]]$down), "\n")

# ---- 3. Intersect with TFs ----
# Reusable function: returns named list of intersections for one contrast
intersect_lists <- function(up, down, reference, ref_name) {
  list(
    up   = intersect(up,   reference),
    down = intersect(down, reference)
  ) |> setNames(paste0(ref_name, c("_up", "_down")))
}

tf_int <- list(
  NCS  = intersect_lists(degs$NCS$up,  degs$NCS$down,  TFs, "TF"),
  NCTM = intersect_lists(degs$NCTM$up, degs$NCTM$down, TFs, "TF")
)

cat("\nIntersection counts:\n")
for (ctr in names(tf_int))
  for (nm in names(tf_int[[ctr]]))
    cat(" ", ctr, nm, ":", length(tf_int[[ctr]][[nm]]), "\n")

# ---- 4. Summary table ----
summary_df <- data.frame(
  Contrast   = rep(names(tf_int), each = 2),
  Category   = rep(c("TF_up", "TF_down"), 2),
  N_genes    = c(
    length(tf_int$NCS$TF_up),   length(tf_int$NCS$TF_down),
    length(tf_int$NCTM$TF_up),  length(tf_int$NCTM$TF_down)
  )
)
write.csv(summary_df, file.path(PATH_TF_TAB, "tf_summary.csv"), row.names = FALSE)

# ---- 5. Volcano plots highlighting TFs ----
# Shared function avoids duplicated ggplot code between NCS and NCTM

make_volcano <- function(toptable, tf_vec, label) {
  
  df <- toptable %>%
    mutate(
      logP = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
      regulation = case_when(
        gene_symbol %in% tf_vec &
          logFC > LFC_THR &
          adj.P.Val < PADJ_THR ~ "TF_up",
        
        gene_symbol %in% tf_vec &
          logFC < -LFC_THR &
          adj.P.Val < PADJ_THR ~ "TF_down",
        
        TRUE ~ "Other"
      )
    )
  
  # ============================================================
  # Top 5 TFs with lowest logFC + Top 5 TFs with highest logFC
  # ============================================================
  
  label_df <- df %>%
    dplyr::filter(
      gene_symbol %in% tf_vec,
      adj.P.Val < PADJ_THR,
      abs(logFC) > LFC_THR
    ) %>%
    dplyr::arrange(logFC) %>%
    dplyr::slice_head(n = 5) %>%
    dplyr::bind_rows(
      df %>%
        dplyr::filter(
          gene_symbol %in% tf_vec,
          adj.P.Val < PADJ_THR,
          abs(logFC) > LFC_THR
        ) %>%
        dplyr::arrange(dplyr::desc(logFC)) %>%
        dplyr::slice_head(n = 5)
    ) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)
  
  # ============================================================
  # Volcano
  # ============================================================
  
  ggplot(df, aes(x = logFC, y = logP)) +
    
    geom_point(
      color = "grey85",
      size = 1.5,
      alpha = 0.6
    ) +
    
    geom_point(
      data = dplyr::filter(df, regulation == "TF_up"),
      color = "#C0392B",
      size = 2.7
    ) +
    
    geom_point(
      data = dplyr::filter(df, regulation == "TF_down"),
      color = "#2980B9",
      size = 2.7
    ) +
    
    geom_vline(
      xintercept = c(-LFC_THR, LFC_THR),
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.3
    ) +
    
    geom_hline(
      yintercept = -log10(PADJ_THR),
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.3
    ) +
    
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = gene_symbol),
      size = 3.2,
      family = "Helvetica",
      fontface = "bold",
      box.padding = 0.4,
      point.padding = 0.3,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.color = "grey50",
      segment.linewidth = 0.3
    ) +
    
    labs(
      title = paste("TF volcano —", label),
      x = expression(log[2] * " Fold Change"),
      y = expression(-log[10] * " (adjusted p-value)")
    ) +
    
    theme_classic(base_size = 11) +
    
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        size = 12
      )
    )
}
p_cs  <- make_volcano(top_NCS,  TFs, "CS vs Normal")
p_ctm <- make_volcano(top_NCTM, TFs, "CTM vs Normal")
p_cs
p_ctm

# ---- 6. Save all objects ----
saveRDS(degs,   file.path(PATH_TF_PROC, "deg_lists.rds"))
saveRDS(tf_int, file.path(PATH_TF_PROC, "tf_intersections.rds"))
cat("Saved: deg_lists.rds | tf_intersections.rds\n")

