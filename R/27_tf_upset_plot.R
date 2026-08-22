# =============================================================================
# Script:      27_tf_upset_plot.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Fernández
# Description: UpSet plot showing which DEG subsets (up/down in CS and CTM)
#              overlap with the HPA transcription factor reference list.
#              Five sets are intersected:
#                - Upregulated CS (N vs CS, logFC >= 1, FDR < 0.05)
#                - Downregulated CS
#                - Upregulated CTM (N vs CTM)
#                - Downregulated CTM
#                - All HPA transcription factors
#              Intersection bars colored navy; set labels colored per set;
#              exported as vectorial PDF.
# Inputs:      data/processed/TFs/deg_lists.rds
#              data/processed/TFs/TFs.rds
# Outputs:     results/figures/TFs/upset_tf_degs.pdf
#              results/figures/TFs/upset_tf_degs.png
# Dependencies: 00_setup.R, 23_tf_hpa_extraction.R, 24_tf_intersections_volcanos.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(ComplexUpset)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble)
library(scales)

PATH_TF_PROC <- file.path(PATH_PROCESSED, "TFs")
PATH_TF_FIG  <- file.path(PATH_FIGURES,   "TFs")

degs <- readRDS(file.path(PATH_TF_PROC, "deg_lists.rds"))
TFs  <- readRDS(file.path(PATH_TF_PROC, "TFs.rds"))

# ---- 1. Design tokens ----
COL_ACTIVE   <- "#1B3A5C"
COL_INACTIVE <- "#DEE4EB"
FONT         <- "Helvetica"

SET_COLORS <- c(
  "Up CS"    = "#C0392B",
  "Down CS"  = "#2874A6",
  "Up CTM"   = "#E67E22",
  "Down CTM" = "#1E8449",
  "TF"       = "#7D3C98"
)
SETS <- names(SET_COLORS)

# ---- 2. Binary presence/absence matrix ----
all_genes <- unique(c(degs$NCS$up, degs$NCS$down,
                      degs$NCTM$up, degs$NCTM$down, TFs))

plot_data <- tibble(gene = all_genes) %>%
  mutate(
    `Up CS`    = gene %in% degs$NCS$up,
    `Down CS`  = gene %in% degs$NCS$down,
    `Up CTM`   = gene %in% degs$NCTM$up,
    `Down CTM` = gene %in% degs$NCTM$down,
    `TF`       = gene %in% TFs
  )

cat("Set sizes:\n")
plot_data %>%
  summarise(across(all_of(SETS), sum)) %>%
  pivot_longer(everything(), names_to = "Set", values_to = "N") %>%
  print()

# ---- 3. UpSet plot ----
int_ann <- (
  intersection_size(
    counts              = TRUE,
    bar_number_threshold = 1,
    text = list(size = 2.3, family = FONT, color = "#2D2D2D", vjust = -0.35)
  ) +
  scale_fill_manual(
    values = c("TRUE" = COL_ACTIVE, "FALSE" = COL_ACTIVE),
    guide  = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = label_comma()) +
  labs(y = "Intersection\nsize", x = NULL) +
  theme_minimal(base_family = FONT, base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "#EAECEF", linewidth = 0.2),
    axis.line.y        = element_line(color = "grey55", linewidth = 0.3),
    axis.ticks.y       = element_line(color = "grey55", linewidth = 0.3),
    axis.text          = element_text(size = 7.5, color = "#2D2D2D"),
    axis.title.y       = element_text(size = 8.5, margin = margin(r = 6))
  )
)

set_ann <- (
  upset_set_size(geom = geom_bar(fill = "grey75", width = 0.62)) +
  labs(x = "Set size") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0)),
                     labels = label_comma()) +
  theme_minimal(base_family = FONT, base_size = 9) +
  theme(
    panel.grid     = element_blank(),
    axis.line.x    = element_line(color = "grey55", linewidth = 0.3),
    axis.ticks.x   = element_line(color = "grey55", linewidth = 0.3),
    axis.title.x   = element_text(size = 8.5, margin = margin(t = 4)),
    axis.text      = element_text(size = 8, color = "#2D2D2D"),
    # Color each set label — rev() because ggplot y-axis runs bottom→top
    # while SETS is defined top→bottom. Adjust if colors appear misaligned.
    axis.text.y    = element_text(colour = rev(unname(SET_COLORS)),
                                  face = "bold", size = 9)
  )
)

p <- upset(
  data               = plot_data,
  intersect          = SETS,
  width_ratio        = 0.24,
  height_ratio       = 0.5,
  sort_sets          = FALSE,
  sort_intersections = "descending",
  base_annotations   = list("Intersection size" = int_ann),
  set_sizes          = set_ann
) +
  scale_color_manual(
    values = c("TRUE" = COL_ACTIVE, "FALSE" = COL_INACTIVE),
    guide  = "none"
  ) +
  scale_fill_manual(
    values = c("TRUE" = COL_ACTIVE, "FALSE" = "white"),
    guide  = "none"
  ) +
  ggtitle("DEG–Transcription Factor Intersections") +
  theme_minimal(base_family = FONT, base_size = 9) +
  theme(
    panel.grid  = element_blank(),
    axis.text   = element_text(size = 8.5, color = "#2D2D2D"),
    plot.title  = element_text(size = 11, face = "bold", hjust = 0.5,
                               color = "#1A1A1A", margin = margin(b = 6)),
    plot.margin = margin(8, 10, 8, 8)
  )

# ---- 4. Export ----
cairo_pdf(
  filename = file.path(PATH_TF_FIG, "upset_tf_degs.pdf"),
  width    = 10, height = 6.5,
  family   = FONT, onefile = FALSE
)
p
dev.off()

ggsave(
  filename = file.path(PATH_TF_FIG, "upset_tf_degs.png"),
  plot     = p, width = 10, height = 6.5, dpi = 600, bg = "white"
)

cat("Saved: upset_tf_degs.pdf / .png →", PATH_TF_FIG, "\n")
