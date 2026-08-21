# =============================================================================
# Script:      25_tf_dorothea_networks.R
# Project:     Transcriptomic profiling in canine mammary cancer (microarray)
# Author:      Karen Rodriguez
# Advisor:     Geysson Javier Hernandez
# Description: Constructs TF-target regulatory networks for CS and CTM using
#              the DoRothEA curated regulon database (confidence levels A/B/C).
#              For each condition:
#                - Retains only edges where the TF AND the target are DEGs
#                - Detects and removes "contradictory" edges: interactions where
#                  the predicted regulatory effect (activator/repressor × TF
#                  direction) contradicts the observed target expression, BUT
#                  only if that target has exactly one regulating TF in the
#                  network (unambiguous contradiction)
#              Outputs are formatted for direct import into Cytoscape:
#                - edges CSV: source, target, mode of regulation, confidence,
#                  edge color, width, and arrow shape
#                - nodes CSV: node type, regulation direction, color, size, shape
# Inputs:      data/processed/TFs/tf_secreted_intersections.rds
#              data/processed/TFs/deg_lists.rds
# Outputs:     results/tables/TFs/network_{CS,CTM}_edges.csv
#              results/tables/TFs/network_{CS,CTM}_nodes.csv
#              results/tables/TFs/network_{CS,CTM}_verification.csv
# Dependencies: 00_setup.R, 24_tf_intersections_volcanos.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(dplyr)
library(dorothea)

PATH_TF_PROC <- file.path(PATH_PROCESSED, "TFs")
PATH_TF_TAB  <- file.path(PATH_TABLES,    "TFs")

tf_int <- readRDS(file.path(PATH_TF_PROC, "tf_intersections.rds"))
degs   <- readRDS(file.path(PATH_TF_PROC, "deg_lists.rds"))

# ---- 1. Load DoRothEA regulon ----
data(dorothea_hs, package = "dorothea")
cat("DoRothEA interactions loaded:", nrow(dorothea_hs), "\n")

# ---- 2. Network builder function ----
# Encapsulates the identical logic used for both CS and CTM,
# eliminating ~120 lines of duplicated code from the original script.
build_network <- function(tf_up, tf_down, genes_up, genes_down,
                          confidence_levels = c("A", "B", "C")) {

  # Clean strings to avoid whitespace mismatches
  tf_up    <- trimws(as.character(tf_up))
  tf_down  <- trimws(as.character(tf_down))
  genes_up <- trimws(as.character(genes_up))
  genes_down <- trimws(as.character(genes_down))

  # Filter regulon: keep only TFs and targets that are DEGs
  net <- dorothea_hs %>%
    filter(
      tf         %in% c(tf_up, tf_down),
      confidence %in% confidence_levels,
      target     %in% c(genes_up, genes_down)
    ) %>%
    mutate(
      tf               = trimws(as.character(tf)),
      target           = trimws(as.character(target)),
      regulacion_tf    = case_when(tf %in% tf_up   ~ "up",
                                   tf %in% tf_down ~ "down",
                                   TRUE            ~ "unknown"),
      regulacion_target = case_when(target %in% genes_up   ~ "up",
                                    target %in% genes_down ~ "down",
                                    TRUE                   ~ "unknown"),
      # An edge is contradictory when the predicted effect (mor × TF direction)
      # opposes the observed target direction
      contradiccion = case_when(
        regulacion_tf == "up"   & mor ==  1 & regulacion_target == "down" ~ TRUE,
        regulacion_tf == "down" & mor ==  1 & regulacion_target == "up"   ~ TRUE,
        regulacion_tf == "up"   & mor == -1 & regulacion_target == "up"   ~ TRUE,
        regulacion_tf == "down" & mor == -1 & regulacion_target == "down" ~ TRUE,
        TRUE ~ FALSE
      )
    )

  # Count unique regulators per target (to assess ambiguity)
  tf_count <- net %>%
    group_by(target) %>%
    summarise(n_tf = n(), .groups = "drop")

  net <- net %>%
    left_join(tf_count, by = "target") %>%
    # Only remove contradictions when target has a SINGLE regulator
    # (unambiguous: no alternative TF can explain the observed expression)
    filter(!(contradiccion & n_tf == 1))

  net
}

# ---- 3. Build networks ----
net_CS  <- build_network(
  tf_up    = tf_int$NCS$TF_up,
  tf_down  = tf_int$NCS$TF_down,
  genes_up = degs$NCS$up,
  genes_down = degs$NCS$down
)

net_CTM <- build_network(
  tf_up    = tf_int$NCTM$TF_up,
  tf_down  = tf_int$NCTM$TF_down,
  genes_up = degs$NCTM$up,
  genes_down = degs$NCTM$down
)

cat("Network CS  — edges:", nrow(net_CS),  "\n")
cat("Network CTM — edges:", nrow(net_CTM), "\n")

# ---- 4. Cytoscape export function ----
# Produces edges + nodes + verification files for one network
export_cytoscape <- function(net, tf_up, tf_down, genes_up, genes_down,
                             prefix, out_dir) {

  # ── Edges ────────────────────────────────────────────────────
  edges <- net %>%
    dplyr::select(source = tf, target, mor, confidence) %>%
    dplyr::mutate(
      source      = as.character(source),
      target      = as.character(target),
      edge_color  = case_when(mor ==  1 ~ "#8B0000",
                              mor == -1 ~ "#00008B",
                              TRUE      ~ "#808080"),
      arrow_shape = if_else(mor == 1, "Arrow", "None"),
      edge_width  = case_when(confidence == "A" ~ 5,
                              confidence == "B" ~ 4,
                              confidence == "C" ~ 3,
                              TRUE              ~ 2)
    )

  # ── Nodes ────────────────────────────────────────────────────
  all_nodes <- unique(c(net$tf, net$target))

  nodes <- data.frame(node = as.character(all_nodes),
                      stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      regulation = case_when(
        node %in% tf_up     ~ "TF_up",
        node %in% tf_down   ~ "TF_down",
        node %in% genes_up  ~ "target_up",
        node %in% genes_down~ "target_down",
        TRUE                ~ "other"
      ),
      color = case_when(
        regulation %in% c("TF_up",     "target_up")   ~ "#8B0000",
        regulation %in% c("TF_down",   "target_down") ~ "#00008B",
        TRUE                                           ~ "#D3D3D3"
      ),
      size  = case_when(regulation %in% c("TF_up","TF_down") ~ 70L, TRUE ~ 40L),
      shape = case_when(regulation %in% c("TF_up","TF_down") ~ "ELLIPSE",
                        TRUE                                  ~ "DIAMOND")
    )

  # ── Verification ─────────────────────────────────────────────
  edge_nodes <- sort(unique(c(edges$source, edges$target)))
  verif <- data.frame(node = edge_nodes) %>%
    left_join(nodes %>% dplyr::select(node, regulation), by = "node")

  match_ok <- length(edge_nodes) == nrow(nodes)
  cat(prefix, "— edges:", nrow(edges),
      "| nodes:", nrow(nodes),
      "| match:", match_ok, "\n")

  # ── Write ─────────────────────────────────────────────────────
  write.csv(edges, file.path(out_dir, paste0("network_", prefix, "_edges.csv")),
            row.names = FALSE, quote = FALSE)
  write.csv(nodes, file.path(out_dir, paste0("network_", prefix, "_nodes.csv")),
            row.names = FALSE, quote = FALSE)
  write.csv(verif, file.path(out_dir, paste0("network_", prefix, "_verification.csv")),
            row.names = FALSE, quote = FALSE)
}

# ---- 5. Export both conditions ----
export_cytoscape(
  net_CS,
  tf_up    = tf_int$NCS$TF_up,
  tf_down  = tf_int$NCS$TF_down,
  genes_up = degs$NCS$up,
  genes_down = degs$NCS$down,
  prefix   = "CS",
  out_dir  = PATH_TF_TAB
)

export_cytoscape(
  net_CTM,
  tf_up    = tf_int$NCTM$TF_up,
  tf_down  = tf_int$NCTM$TF_down,
  genes_up = degs$NCTM$up,
  genes_down = degs$NCTM$down,
  prefix   = "CTM",
  out_dir  = PATH_TF_TAB
)

cat("\nAll Cytoscape files saved to:", PATH_TF_TAB, "\n")

