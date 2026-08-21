# Transcriptomic Profiling of Canine Mammary Cancer via Microarray

R pipeline for differential gene expression analysis from microarray data
(Affymetrix CanGene-1.0-ST) in canine mammary cancer tissue, comparing
normal tissue (N) against two carcinoma histological subtypes: simple
carcinoma (CS) and carcinoma in mixed tumor (CTM).

**Author:** Karen Rodríguez
**Advisor:** Geysson Javier Hernández

##Data Availability Note

This repository is part of a manuscript in preparation. Raw data (`.CEL`
files), the sample metadata table, and study-specific results are **not
included** until the associated article is published. This repository
documents the **pipeline and methodology only**; the code is fully
functional and reproducible once the corresponding input files are placed
in `data/raw/`.

## Pipeline overview

```
Raw .CEL files
   │
   ▼
QC (pre-normalization: pseudoimages, boxplots, density plots, RLE via PLM)
   │
   ▼
RMA normalization (background correction + quantile normalization +
                    median polish summarization, transcript-cluster level)
   │
   ├──► Sample metadata (targets) loading
   │
   ▼
Exploratory PCA (top 5,000 most variable transcript clusters)
   │
   ▼
Probe annotation → gene symbol mapping
   │
   ▼
Linear model (limma) + contrasts (NCS, NCTM) + eBayes moderation
   │
   ▼
DEG extraction (|logFC| ≥ 1, adj.P.Val < 0.05; up / down / total per contrast)
   │
   ├──► PCA on DEGs (combined + per-contrast) → exported for GraphPad Prism
   │
   ├──► DEG hierarchical clustering heatmap (ComplexHeatmap)
   │
   ├──► Gene list export for volcano/waterfall visualization
   │
   ├──► Canine-to-human ortholog mapping (biomaRt, one-to-one orthologs only)
   │        │
   │        ▼
   │     GO Biological Process gene set construction (org.Hs.eg.db, GO.db)
   │        │
   │        ▼
   │     GSVA enrichment scoring (per-sample pathway activity)
   │        │
   │        ▼
   │     Differential pathway enrichment (limma on GSVA scores, NCS/NCTM)
   │        │
   │        ▼
   │     Pathway curation: full significant term tables (go_gene_table_up/down)
   │     → top-30 term selection + manual hallmark assignment (GO_selected_*)
   │        │
   │        ├──► Hub gene identification (genes shared across >1 GO term,
   │        │     intersected with DEGs) → alluvial diagrams (Gene → GO → Hallmark)
   │        │
   │        └──► Supplementary boxplot tables (GO-DEG intersection, manually
   │              curated gene sets per hallmark) → exported for GraphPad Prism
   │
   └──► Transcription factor (TF) analysis
            │
            ▼
        HPA download + TF reference list extraction
            │
            ▼
        TF–DEG intersection (up/down × NCS/NCTM) + TF-highlighted volcano plots
            │
            ├──► TF expression heatmaps (pheatmap, z-score, 4-cluster cutree)
            │      → exported for GraphPad Prism
            │
            ├──► ChEA 2022 regulon enrichment (via Enrichr) → heatmaps
            │
            └──► DoRothEA regulatory network construction (confidence A/B/C,
                  contradiction filtering) → Cytoscape 3.10.4 export
```

## Repository structure

```
├── R/                          Numbered scripts, in execution order
├── functions/                  Reusable helper functions
├── data/
│   ├── raw/                    Input data (not version-controlled, see .gitignore)
│   └── processed/               Intermediate .rds objects passed between scripts
├── results/
│   ├── tables/                  .csv / .xlsx outputs
│   └── figures/                 .pdf outputs
└── run_pipeline.R               Runs the full pipeline in order
```

## How to run it

1. Open the project in RStudio (`.Rproj`).
2. Place input files in `data/raw/` (see each script's header comment for
   the exact files it expects).
3. Run `R/00_setup.R` once to install/load dependencies.
4. Run the scripts in numeric order, or use `run_pipeline.R`.

Each script saves its intermediate results as `.rds` objects in
`data/processed/`, which the next script reloads with `readRDS()`. This
avoids re-running computationally expensive steps (such as RMA
normalization or GSVA) every time a different part of the analysis needs
to be revisited.

## Methodology summary

| Step | Method | Package(s) |
|---|---|---|
| Normalization | RMA (background correction + quantile normalization + median polish summarization, transcript-cluster level) | `oligo` |
| Annotation | Mapping `transcript_cluster_id` to `gene_symbol` via the official Affymetrix annotation file | `dplyr` |
| Differential expression | Linear model + contrasts (NCS, NCTM) + empirical Bayes moderation (eBayes) | `limma` |
| DEG criterion | \|logFC\| ≥ 1 and Benjamini–Hochberg adjusted p-value (FDR) < 0.05 | — |
| Dimensionality reduction | PCA on top variable transcript clusters and on DEGs | `prcomp` (base R) |
| Ortholog mapping | Canine-to-human one-to-one orthologs | `biomaRt` |
| Pathway enrichment | GO Biological Process gene sets + Gene Set Variation Analysis (GSVA) | `org.Hs.eg.db`, `GO.db`, `GSVA` |
| Hub gene / hallmark curation | GO-term/DEG intersection, manual hallmark and term curation | `dplyr`, `tidyr` |
| Alluvial visualization | Gene → GO Biological Process → Hallmark diagrams | `ggalluvial` |
| Transcription factor identification | Human Protein Atlas reference list, intersected with DEGs | `dplyr` |
| Regulatory network reconstruction | DoRothEA regulon (confidence A/B/C), contradiction filtering | `dorothea` |
| TF regulon enrichment | ChEA 2022 database via Enrichr | — |
| Visualization | Heatmaps with hierarchical clustering, PCA plots, volcano plots | `ComplexHeatmap`, `pheatmap`, `ggplot2`, `ggrepel` |
| Network visualization | Edge/node CSV export for external rendering | Cytoscape 3.10.4 |
| Final figure assembly | Panel alignment, font/axis correction | Inkscape 1.4.4 |

## Main dependencies

`limma`, `oligo`, `pd.cangene.1.0.st`, `dplyr`, `tidyr`, `stringr`,
`ggplot2`, `ggrepel`, `ggalluvial`, `matrixStats`, `ComplexHeatmap`,
`circlize`, `pheatmap`, `RColorBrewer`, `biomaRt`, `org.Hs.eg.db`, `GO.db`,
`AnnotationDbi`, `GSVA`, `Biobase`, `sva`, `dorothea`, `readxl`, `readr`,
`here`.