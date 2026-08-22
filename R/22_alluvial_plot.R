# =============================================================================
# Script:      22_alluvial_plot_hub_genes.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Fernández
# Description: Publication-quality 3-axis alluvial plots (Hub Gene -> GO
#              Biological Process -> Hallmark of Cancer), split into two
#              SEPARATE figures: one for UP-regulated hub genes/processes,
#              one for DOWN-regulated. Each figure colors ribbons by GO
#              process (own palette per direction) for maximum readability
#              within a single-direction panel. Panoramic layout, target
#              size per panel: 6.16in x 2.18in.
# Inputs:      results/tables/GSVA_analysis/alluvial_info_annotated.rds
#              (or the hand-curated hub_data table below)
# Outputs:     results/figures/alluvial_hub_genes_UP.pdf / .png
#              results/figures/alluvial_hub_genes_DOWN.pdf / .png
# Dependencies: ggplot2, ggalluvial, dplyr, tidyr, stringr, showtext, sysfonts
# =============================================================================

# ---- 0. Packages ----
required_pkgs <- c("ggplot2", "ggalluvial", "dplyr", "tidyr",
                   "stringr", "showtext", "sysfonts")
missing_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(missing_pkgs)) install.packages(missing_pkgs)

library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyr)
library(stringr)
library(showtext)
library(sysfonts)

# ---- 1. Font setup ----
sysfonts::font_add_google("Arimo", "Arial")
showtext_auto()
showtext_opts(dpi = 600)

# ---- 2. Hub-gene source table -------------------------------------------
hub_data <- tibble::tribble(
  ~Regulation, ~GO_process_raw,
  ~Gene, ~Hallmark_raw,
  
  "DOWN", "natural killer cell activation; regulation of T cell differentiation",
  "CD2", "Immune evasion",
  
  "DOWN", "regulation of calcium ion transport; T cell chemotaxis",
  "CXCL12", "Immune evasion; Invasion & metastasis",
  
  "DOWN", "cGMP-mediated signaling; negative regulation of vascular permeability",
  "PDE2A", "Angiogenesis; Proliferation",
  
  "DOWN", "leukocyte chemotaxis; T cell migration",
  "S1PR1", "Immune evasion",
  
  "UP", "establishment of planar polarity; Wnt signaling pathway, planar cell polarity pathway",
  "WNT5A", "Invasion & metastasis; Proliferation",
  
  "UP", "intracellular transport",
  "APOE", "Metabolism",
  
  "UP", "glutathione metabolic process; zinc ion transmembrane transport",
  "SLC1A1", "Death resistance; Plasticity",
  
  "UP", "glutathione metabolic process; negative regulation of ferroptosis",
  "SLC7A11", "Death resistance"
)

# ---- 3. Label-simplification lookups -------------------------------------
go_label_map <- c(
  "natural killer cell activation"                     = "NK cell activation",
  "regulation of T cell differentiation"                = "T-cell differentiation",
  "regulation of calcium ion transport"                 = "Ca\u00b2\u207a transport reg.",
  "T cell chemotaxis"                                   = "T-cell chemotaxis",
  "cGMP-mediated signaling"                             = "cGMP signaling",
  "negative regulation of vascular permeability"        = "Vascular permeability",
  "leukocyte chemotaxis"                                = "Leukocyte chemotaxis",
  "T cell migration"                                    = "T-cell migration",
  "establishment of planar polarity"                    = "Planar polarity",
  "Wnt signaling pathway, planar cell polarity pathway" = "Wnt/PCP signaling",
  "intracellular transport"                             = "Intracellular transport",
  "glutathione metabolic process"                       = "Glutathione metabolism",
  "zinc ion transmembrane transport"                    = "Zn\u00b2\u207a transport",
  "negative regulation of ferroptosis"                  = "Ferroptosis inhibition"
)

hallmark_label_map <- c(
  "Immune evasion"        = "Immune Evasion",
  "Invasion & metastasis" = "Invasion & Metastasis",
  "Angiogenesis"           = "Angiogenesis",
  "Proliferation"          = "Proliferation",
  "Metabolism"             = "Metabolism",
  "Death resistance"       = "Death Resistance",
  "Plasticity"              = "Plasticity"
)

# ---- 4. Expand semicolon-separated cells into long many-to-many format ---
hub_long_all <- hub_data %>%
  separate_rows(GO_process_raw, sep = ";\\s*") %>%
  separate_rows(Hallmark_raw,   sep = ";\\s*") %>%
  mutate(
    GO_process = recode(GO_process_raw, !!!go_label_map, .default = GO_process_raw),
    Hallmark   = recode(Hallmark_raw,   !!!hallmark_label_map, .default = Hallmark_raw)
  ) %>%
  distinct(Regulation, GO_process, Gene, Hallmark)

# ---- 5. Per-GO-term color palettes (warm = UP, cool = DOWN) --------------
warm_pal <- c("#B2182B", "#D6604D", "#E9967A", "#C97B63",
              "#A15843", "#DDA15E", "#8C2F39")
cool_pal <- c("#2166AC", "#4393C3", "#67A9CF", "#5B7C99",
              "#3E5C76", "#7FA6B5", "#4C6785", "#89A3B2")

# ---- 6. Reusable function: build one direction's alluvial plot ----------
build_direction_plot <- function(direction, palette_colors) {
  
  hub_long <- hub_long_all %>%
    dplyr::filter(Regulation == direction) %>%
    droplevels()
  
  gene_order     <- sort(unique(as.character(hub_long$Gene)))
  go_order       <- sort(unique(as.character(hub_long$GO_process)))
  hallmark_order <- hub_long %>% dplyr::count(Hallmark, sort = TRUE) %>% pull(Hallmark)
  
  hub_long <- hub_long %>%
    dplyr::mutate(
      Gene       = factor(Gene,       levels = gene_order),
      GO_process = factor(GO_process, levels = go_order),
      Hallmark   = factor(Hallmark,   levels = hallmark_order)
    )
  
  go_palette <- setNames(palette_colors[seq_along(go_order)], go_order)
  
  ggplot(hub_long,
         aes(axis1 = Gene, axis2 = GO_process, axis3 = Hallmark)) +
    geom_alluvium(
      aes(fill = GO_process),
      width = 1/40,            # thin tick, not a wide box
      alpha = 0.6,
      linewidth = 0,           # no ribbon border -> softer look
      knot.pos = 0.5,
      curve_type = "sigmoid"
    ) +
    geom_stratum(
      width = 1/40,
      fill = "grey15",         # solid thin black tick, like the reference
      color = "grey15",
      linewidth = 0
    ) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      family = "Arial", size = 1.9, lineheight = 0.75, color = "grey15",
      nudge_x = 0.12, hjust = 0     # push label off the thin tick, like reference
    ) +
    scale_x_discrete(
      limits = c("Hub Gene", "GO Biological Process", "Hallmark of Cancer"),
      expand = c(0.1, 0.1)
    ) +
    scale_fill_manual(values = go_palette) +
    theme_void(base_family = "Arial") +
    theme(legend.position = "none")
}
# ---- 7. Build both plots ----
p_up   <- build_direction_plot("UP",   warm_pal)
p_down <- build_direction_plot("DOWN", cool_pal)

p_up
p_down

# ---- 8. Export both, same panoramic target size ----
ggsave(
  filename = file.path(PATH_GSVA_FIG, "alluvial_hub_genes_UP.pdf"),
  plot = p_up, width = 6.16, height = 2.18, units = "in",
  device = cairo_pdf
)
ggsave(
  filename = file.path(PATH_GSVA_FIG, "alluvial_hub_genes_UP.png"),
  plot = p_up, width = 6.16, height = 2.18, units = "in",
  dpi = 600, bg = "white"
)

ggsave(
  filename = file.path(PATH_GSVA_FIG, "alluvial_hub_genes_DOWN.pdf"),
  plot = p_down, width = 6.16, height = 2.18, units = "in",
  device = cairo_pdf
)
ggsave(
  filename = file.path(PATH_GSVA_FIG, "alluvial_hub_genes_DOWN.png"),
  plot = p_down, width = 6.16, height = 2.18, units = "in",
  dpi = 600, bg = "white"
)

cat("Saved: alluvial_hub_genes_UP.pdf/.png, alluvial_hub_genes_DOWN.pdf/.png\n")
