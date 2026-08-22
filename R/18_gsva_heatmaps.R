# =============================================================================
# Script:      18_gsva_heatmaps.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Fernández
# Description: Produces publication-quality ComplexHeatmap figures for GSVA
#              pathway enrichment scores (logFC from limma on GSVA matrix).
#              Two separate heatmaps are generated:
#                - UP:   GO terms upregulated in CS and/or CTM vs N
#                - DOWN: GO terms downregulated in CS and/or CTM vs N
#              Both are exported as vectorial PDFs to results/figures/GSVA_analysis/
#
#              Design principles:
#                - Helvetica font throughout
#                - Heatmap matrix maximizes figure width; labels are compact
#                - Custom Nature-style color gradients (no RColorBrewer defaults)
#                - Nearly invisible cell borders (#F0F0F0, lwd = 0.2)
#                - Subtle inter-hallmark spacing via row_split (titles hidden)
#                - Minimal, publication-style legend
#                - Identical dimensions for UP and DOWN for panel assembly
#
# Inputs:      results/tables/GSVA_analysis/GO_selected_up_gsva_3.xlsx
#              results/tables/GSVA_analysis/GO_selected_down_gsva_3.xlsx
# Outputs:     results/figures/GSVA_analysis/heatmap_UP_pathways.pdf
#              results/figures/GSVA_analysis/heatmap_DOWN_pathways.pdf
# Dependencies: 00_setup.R, 19_gsva_pathway_tables.R
# =============================================================================

source(here::here("R", "00_setup.R"))

library(ComplexHeatmap)
library(circlize)
library(readxl)
library(tibble)
library(dplyr)

# =============================================================================
# SECTION 1: GLOBAL HEATMAP THEME
# Applied once before any Heatmap() call; affects all heatmaps in the session
# =============================================================================

ht_opt(
  TITLE_PADDING           = unit(c(1.5, 1.5), "mm"),
  heatmap_border          = FALSE,
  annotation_border       = FALSE,
  legend_border           = FALSE,
  legend_gap              = unit(4, "mm")
)

# =============================================================================
# SECTION 2: GO TERM NAME DICTIONARIES
# All label simplifications in one place — change here, reflects everywhere
# =============================================================================

go_names_UP <- c(
  
  "Golgi lumen acidification" =
    "Golgi acidification",
  
  "protein folding in endoplasmic reticulum" =
    "Protein folding in ER",
  
  "pentose-phosphate shunt" =
    "Pentose phosphate pathway",
  
  "endoplasmic reticulum mannose trimming" =
    "ER mannose trimming",
  
  "retrograde protein transport, ER to cytosol" =
    "ER-to-cytosol transport",
  
  "electron transport chain" =
    "Electron transport chain",
  
  "calcium import into the mitochondrion" =
    "Mitochondrial Ca²⁺ import",
  
  "branched-chain amino acid catabolic process" =
    "BCAA catabolism",
  
  "heme biosynthetic process" =
    "Heme biosynthesis",
  
  "apical protein localization" =
    "Apical protein localization",
  
  "NAD+ metabolic process" =
    "NAD⁺ metabolism",
  
  "glutathione metabolic process" =
    "Glutathione metabolism",
  
  "positive regulation of mitophagy" =
    "Mitophagy activation",
  
  "chaperone-mediated protein complex assembly" =
    "Chaperone complex assembly",
  
  "mitochondrial electron transport, NADH to ubiquinone" =
    "Mitochondrial ETC I",
  
  "membrane protein intracellular domain proteolysis" =
    "Membrane protein proteolysis",
  
  "GPI anchor biosynthetic process" =
    "GPI anchor biosynthesis",
  
  "establishment of planar polarity" =
    "Planar polarity",
  
  "response to reactive oxygen species" =
    "Response to ROS",
  
  "protein N-linked glycosylation via asparagine" =
    "N-linked glycosylation",
  
  "NADP+ metabolic process" =
    "NADP⁺ metabolism",
  
  "zinc ion transmembrane transport" =
    "Zn²⁺ transport",
  
  "negative regulation of ferroptosis" =
    "Ferroptosis inhibition",
  
  "intrinsic apoptotic signaling pathway in response to oxidative stress" =
    "ROS-induced intrinsic apoptosis",
  
  "inner mitochondrial membrane organization" =
    "Inner mitochondrial membrane",
  
  "postreplication repair" =
    "Postreplication repair",
  
  "Wnt signaling pathway, planar cell polarity pathway" =
    "Wnt/PCP signaling",
  
  "maturation of LSU-rRNA from tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA)" =
    "LSU-rRNA maturation"
  
)

go_names_DOWN <- c(
  
  "alpha-beta T cell differentiation" =
    "αβ T cell differentiation",
  
  "calcineurin-NFAT signaling cascade" =
    "Calcineurin–NFAT signaling",
  
  "centriole-centriole cohesion" =
    "Centriole cohesion",
  
  "cGMP-mediated signaling" =
    "cGMP signaling",
  
  "diacylglycerol metabolic process" =
    "Diacylglycerol metabolism",
  
  "glycerolipid metabolic process" =
    "Glycerolipid metabolism",
  
  "inflammatory response to antigenic stimulus" =
    "Antigen-driven inflammation",
  
  "leukocyte chemotaxis" =
    "Leukocyte chemotaxis",
  
  "lipid phosphorylation" =
    "Lipid phosphorylation",
  
  "lymph vessel development" =
    "Lymphatic vessel development",
  
  "lymphangiogenesis" =
    "Lymphangiogenesis",
  
  "natural killer cell activation" =
    "NK cell activation",
  
  "negative regulation of myeloid cell differentiation" =
    "Myeloid differentiation inhibition",
  
  "negative regulation of T cell apoptotic process" =
    "T-cell apoptosis inhibition",
  
  "nucleotide-binding oligomerization domain containing 2 signaling pathway" =
    "NOD2 signaling",
  
  "positive regulation of cholesterol efflux" =
    "Cholesterol efflux",
  
  "positive regulation of MHC class II biosynthetic process" =
    "MHC-II biosynthesis",
  
  "positive regulation of natural killer cell mediated cytotoxicity" =
    "NK cytotoxicity",
  
  "positive regulation of vascular endothelial growth factor receptor signaling pathway" =
    "VEGFR signaling",
  
  "regulation of calcium ion transport" =
    "Ca²⁺ transport regulation",
  
  "regulation of T cell differentiation" =
    "T-cell differentiation",
  
  "T cell migration" =
    "T-cell migration",
  
  "thrombin-activated receptor signaling pathway" =
    "Thrombin receptor signaling",
  
  "triglyceride catabolic process" =
    "Triglyceride catabolism",
  
  "DNA methylation-dependent constitutive heterochromatin formation" =
    "DNA methylation-dependent\nheterochromatin",
  
  "positive regulation of isotype switching to IgG isotypes" =
    "IgG class switching"
  
)
  

# Hallmark unification maps (old label -> canonical label)
hallmark_map_UP <- c(
  
  "Unlocking phenotypic plasticity"      = "Plasticity",
  "Stress adaptation"                    = "Stress Adaptation",
  "Sustaining proliferative signaling"   = "Proliferation",
  "Resisting cell death"                 = "Death Resistance",
  "Genome instability and mutation"      = "Genome Instability",
  "Deregulating cellular energetics"     = "Metabolism",
  "Avoiding immune destruction"          = "Immune Evasion",
  "Activating invasion and metastasis"   = "Invasion & Metastasis"
  
)
hallmark_map_DOWN <- c(
  "Tumor-promoting inflammation"       = "Inflammation",
  "Sustaining proliferative signaling" = "Proliferation",
  "Stress adaptation"                  = "Stress Adaptation",
  "Resisting cell death"               = "Death Resistance",
  "Genome instability and mutation"    = "Genome Instability",
  "Deregulating cellular energetics"   = "Metabolism",
  "Avoiding immune destruction"        = "Immune Evasion",
  "Activating invasion and metastasis" = "Invasion & Metastasis",
  "Unlocking phenotypic plasticity"    = "Plasticity",
  "Inducing angiogenesis" = "Angiogenesis"
)

hallmark_colors_UP <- c(
  "Proliferation"            = "#7B2D8B",
  "Death Resistance"     = "black",
  "Stress Adaptation"        = "#8C4F2D",
  "Genome Instability"       = "#4A4A4A",
  "Metabolism"               = "#C8A900",
  "Immune Evasion"           = "#2856A8",
  "Invasion & Metastasis"    = "#D4700A",
  "Plasticity"               = "#1B9E77",
  "Inflammation"             = "#B83040"
)

hallmark_order_UP <- c(
  "Proliferation",
    "Stress Adaptation",
  "Genome Instability",
  "Metabolism",
  "Immune Evasion",
  "Invasion & Metastasis",
  "Plasticity",
  "Inflammation",
  'Death Resistance'
)

hallmark_order_DOWN <- c(
  "Proliferation", "Inflammation", "Angiogenesis",
  "Genome Instability", "Metabolism", "Immune Evasion", "Invasion & Metastasis"
)


hallmark_colors_DOWN <- c(
  "Inflammation"         = "#B83040",
  "Angiogenesis"         = "#D4608A",
  "Genome Instability"   = "#4A4A4A",
  "Metabolism"           = "#C8A900",
  "Immune Evasion"       = "#2856A8",
  "Plasticity"           = "#1B9E77"
)

hallmark_order_DOWN <- c(
  "Inflammation",
  "Angiogenesis",
  "Genome Instability",
  "Metabolism",
  "Immune Evasion",
  "Plasticity"
)
# =============================================================================
# SECTION 3: MODULAR FUNCTIONS
# =============================================================================

# ---- prepare_heatmap_matrix() ----
# Loads an Excel file, applies name simplification, unifies hallmarks, orders rows
# Returns: list(matrix, hallmark_vector)
prepare_heatmap_matrix <- function(xlsx_path,
                                   go_names_dict,
                                   hallmark_map,
                                   hallmark_order) {
  df <- read_excel(xlsx_path) %>%
    filter(!is.na(HALLMARK)) %>%
    mutate(
      logFC_CS  = as.numeric(logFC_CS),
      logFC_CTM = as.numeric(logFC_CTM),
      # Unify hallmark labels using the map; keep original if not in map
      HALLMARK  = ifelse(HALLMARK %in% names(hallmark_map),
                         hallmark_map[HALLMARK],
                         HALLMARK)
    )
  
  # Build expression matrix (GO_ID as rownames)
  mat <- df %>%
    dplyr::select(GO_ID, logFC_CS, logFC_CTM) %>%
    column_to_rownames("GO_ID") %>%
    as.matrix()
  colnames(mat) <- c("CS", "CTM")
  
  # Apply simplified GO labels
  simplified <- go_names_dict[df$GO_Term]
  unmapped   <- is.na(simplified)
  if (any(unmapped)) {
    warning(sum(unmapped), " GO term(s) have no simplified name; original label used.")
    simplified[unmapped] <- df$GO_Term[unmapped]
  }
  rownames(mat) <- simplified
  
  # Order rows by hallmark group
  ord <- order(factor(df$HALLMARK, levels = hallmark_order))
  list(
    matrix           = mat[ord, , drop = FALSE],
    hallmark_vector  = factor(df$HALLMARK[ord], levels = hallmark_order)
  )
}

# ---- create_color_palette() ----
# Custom Nature-style gradient palettes using colorRamp2()
# direction = "up"  : white → light pink → magenta → dark crimson
# direction = "down": dark forest → teal → light mint → white
create_color_palette <- function(mat, direction = c("up", "down")) {
  direction <- match.arg(direction)
  
  if (direction == "up") {
    # All values should be positive for UP
    lo <- max(0, quantile(mat, 0.02, na.rm = TRUE))
    hi <- quantile(mat, 0.98, na.rm = TRUE)
    breaks <- seq(lo, hi, length.out = 7)
    colors <- c("#FFFFFF", "#FFE0EC", "#F9A8C9", "#F06292",
                "#D81B60", "#AD1457", "#6A0032")
  } else {
    # All values should be negative for DOWN
    lo <- quantile(mat, 0.02, na.rm = TRUE)
    hi <- min(0, quantile(mat, 0.98, na.rm = TRUE))
    breaks <- seq(lo, hi, length.out = 7)
    colors <- c("#00352A", "#00695C", "#26A69A", "#80CBC4",
                "#B2DFDB", "#E0F2F1", "#FFFFFF")
  }
  
  colorRamp2(breaks, colors)
}

# ---- create_row_annotation() ----
# Creates right-side hallmark annotation bar
create_row_annotation <- function(hallmark_vector, hallmark_colors) {
  rowAnnotation(
    Hallmark = hallmark_vector,
    col      = list(Hallmark = hallmark_colors),
    show_annotation_name = FALSE,
    show_legend          = TRUE,
    simple_anno_size     = unit(0.22, "cm"),   # narrow bar
    annotation_legend_param = list(
      title    = "Hallmark",
      title_gp = gpar(fontsize = 7.5, fontface = "bold", fontfamily = "Helvetica"),
      labels_gp = gpar(fontsize = 7, fontfamily = "Helvetica"),
      grid_height = unit(3, "mm"),
      grid_width  = unit(3, "mm")
    )
  )
}

# ---- make_heatmap() ----
# Assembles the Heatmap object with all publication-quality aesthetic settings
make_heatmap <- function(mat,
                         col_fun,
                         row_ha,
                         legend_name  = "logFC",
                         hallmark_vec,
                         cell_h_cm    = 0.37,    # row height in cm
                         n_legend_breaks = 5) {
  
  val_range <- range(mat, na.rm = TRUE)
  legend_at  <- seq(val_range[1], val_range[2], length.out = n_legend_breaks)
  
  Heatmap(
    mat,
    name = legend_name,
    col  = col_fun,
    
    # ---- Rows ----
    cluster_rows    = FALSE,
    row_split       = hallmark_vec,
    row_gap         = unit(1.0, "mm"),       # spacing between hallmark blocks
    row_title       = NULL,                   # hide split group titles
    row_title_gp    = gpar(fontsize = 0, col = "transparent"),
    row_names_gp    = gpar(fontsize = 8.5, fontfamily = "Helvetica"),
    row_names_side  = "right",
    row_names_max_width = unit(5.2, "cm"),   # compact labels
    
    # ---- Columns ----
    cluster_columns  = FALSE,
    column_names_gp  = gpar(fontsize = 9, fontface = "bold", fontfamily = "Helvetica"),
    column_names_rot = 0,
    column_names_centered = TRUE,
    
    # ---- Right annotation (Hallmark bar) ----
    right_annotation = row_ha,
    
    # ---- Cell aesthetics ----
    rect_gp = gpar(col = "#F0F0F0", lwd = 0.2),  # nearly invisible borders
    border  = FALSE,
    
    # ---- Dimensions ----
    # `width`  = width  of the MATRIX CELLS ONLY (2 columns of color).
    # `height` = height of the MATRIX CELLS ONLY (rows × row height).
    # ComplexHeatmap adds row names, annotations, and dendrograms outside
    # these values automatically. Using `width` (not `heatmap_width`) avoids
    # the "body width is negative" error that occurs when row_names_max_width
    # would exceed a too-small heatmap_width.
    width  = unit(2.2, "cm"),
    height = unit(nrow(mat) * cell_h_cm, "cm"),
    
    # ---- Legend ----
    heatmap_legend_param = list(
      title         = legend_name,
      title_gp      = gpar(fontsize = 7.5, fontface = "bold", fontfamily = "Helvetica"),
      labels_gp     = gpar(fontsize = 7, fontfamily = "Helvetica"),
      at            = round(legend_at, 2),
      border        = NA,
      grid_width    = unit(2.5, "mm"),
      legend_height = unit(2.8, "cm"),
      direction     = "vertical"
    )
  )
}

# ---- export_heatmap() ----
# Draws heatmap and saves as vectorial PDF via cairo_pdf
# width and height in cm; cairo_pdf takes inches
export_heatmap <- function(ht, filepath, width_cm, height_cm,
                           legend_side = "left", annot_side = "right") {
  cairo_pdf(
    filename = filepath,
    width    = width_cm  / 2.54,
    height   = height_cm / 2.54,
    family   = "Helvetica",
    onefile  = FALSE
  )
  draw(
    ht,
    heatmap_legend_side      = legend_side,
    annotation_legend_side   = annot_side,
    merge_legend             = FALSE,
    padding                  = unit(c(4, 4, 4, 4), "mm")   # top, right, bottom, left
  )
  dev.off()
  cat("Saved:", filepath, "\n")
}

# =============================================================================
# SECTION 4: BUILD UP HEATMAP
# =============================================================================

cat("--- Building UP pathway heatmap ---\n")

up_data <- prepare_heatmap_matrix(
  xlsx_path     = file.path(PATH_GSVA_TAB, "GO_selected_up_gsva_4.xlsx"),
  go_names_dict = go_names_UP,
  hallmark_map  = hallmark_map_UP,
  hallmark_order = hallmark_order_UP
)

col_fun_up <- create_color_palette(up_data$matrix, direction = "up")

row_ha_up  <- create_row_annotation(up_data$hallmark_vector, hallmark_colors_UP)

ht_up <- make_heatmap(
  mat          = up_data$matrix,
  col_fun      = col_fun_up,
  row_ha       = row_ha_up,
  legend_name  = "logFC",
  hallmark_vec = up_data$hallmark_vector
)

# Preview in RStudio viewer
draw(ht_up, heatmap_legend_side = "left", annotation_legend_side = "right",
     padding = unit(c(4, 4, 4, 4), "mm"))

cat("UP heatmap: ", nrow(up_data$matrix), "GO terms |",
    length(levels(up_data$hallmark_vector)), "hallmarks\n")

# =============================================================================
# SECTION 5: BUILD DOWN HEATMAP
# =============================================================================

cat("--- Building DOWN pathway heatmap ---\n")

down_data <- prepare_heatmap_matrix(
  xlsx_path      = file.path(PATH_GSVA_TAB, "GO_selected_down_gsva_4.xlsx"),
  go_names_dict  = go_names_DOWN,
  hallmark_map   = hallmark_map_DOWN,
  hallmark_order = hallmark_order_DOWN
)

col_fun_down <- create_color_palette(down_data$matrix, direction = "down")

row_ha_down  <- create_row_annotation(down_data$hallmark_vector, hallmark_colors_DOWN)

ht_down <- make_heatmap(
  mat          = down_data$matrix,
  col_fun      = col_fun_down,
  row_ha       = row_ha_down,
  legend_name  = "logFC",
  hallmark_vec = down_data$hallmark_vector
)

draw(ht_down, heatmap_legend_side = "left", annotation_legend_side = "right",
     padding = unit(c(4, 4, 4, 4), "mm"))

cat("DOWN heatmap:", nrow(down_data$matrix), "GO terms |",
    length(levels(down_data$hallmark_vector)), "hallmarks\n")

# =============================================================================
# SECTION 6: EXPORT VECTORIAL PDFs
# =============================================================================
# Page dimensions are computed from actual data so nothing is cropped.
# Both heatmaps use identical widths for panel alignment in Illustrator/Inkscape.
#
# Layout per PDF (left to right):
#   legend (~2 cm) | matrix (2.2 cm) | labels (~5.5 cm) | annot bar (~0.4 cm)
#   Total width ≈ 10.5 cm; adjust if your label set is longer
# =============================================================================

PAGE_WIDTH_CM  <- 10.8   # fixed for both panels to allow perfect alignment

# PDF page heights: cell area + gaps between hallmark groups + column names
# + top/bottom margins + legend space outside the heatmap body.
# Rule of thumb: ~0.18 cm per gap (1.8 mm) × number of hallmark groups,
# plus ~2.5 cm for column labels, padding, and legend.
n_gaps_up   <- nlevels(up_data$hallmark_vector)   - 1
n_gaps_down <- nlevels(down_data$hallmark_vector)  - 1

up_height_cm   <- nrow(up_data$matrix)   * 0.37 + n_gaps_up   * 0.18 + 3.5
down_height_cm <- nrow(down_data$matrix) * 0.37 + n_gaps_down * 0.18 + 3.5

export_heatmap(
  ht       = ht_up,
  filepath = file.path(PATH_GSVA_FIG, "heatmap_UP_pathways_btt.pdf"),
  width_cm  = PAGE_WIDTH_CM,
  height_cm = up_height_cm
)

export_heatmap(
  ht       = ht_down,
  filepath = file.path(PATH_GSVA_FIG, "heatmap_DOWN_pathways_btt.pdf"),
  width_cm  = PAGE_WIDTH_CM,
  height_cm = down_height_cm
)

# =============================================================================
# SECTION 7: SUMMARY
# =============================================================================

cat("\n=== EXPORT COMPLETE ===\n")
cat("UP heatmap  →", nrow(up_data$matrix),   "GO terms across",
    nlevels(up_data$hallmark_vector),   "hallmarks\n")
cat("DOWN heatmap→", nrow(down_data$matrix),  "GO terms across",
    nlevels(down_data$hallmark_vector),  "hallmarks\n")
cat("\nFiles saved to:", PATH_GSVA_FIG, "\n")
cat("\nDesign note:\n")
cat("  Both PDFs are identical width (", PAGE_WIDTH_CM, "cm) for panel alignment.\n")
cat("  Import both into Illustrator/Inkscape, stack vertically or side by side.\n")
cat("  All text is real vector text — fully editable without rasterization.\n")
