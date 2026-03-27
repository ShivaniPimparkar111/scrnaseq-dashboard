# scRNA-seq Analysis Dashboard

![R](https://img.shields.io/badge/R-%3E%3D4.2-276DC3?style=flat-square&logo=r&logoColor=white)
![Shiny](https://img.shields.io/badge/Shiny-1.7%2B-blue?style=flat-square&logo=rstudio&logoColor=white)
![Seurat](https://img.shields.io/badge/Seurat-v5-4e9af1?style=flat-square)
![License](https://img.shields.io/badge/license-MIT-green?style=flat-square)

Interactive R Shiny dashboard for single-cell RNA sequencing data exploration and analysis using the Seurat framework.

> **Note:** This dashboard is designed for initial exploration of scRNA-seq data. Results should be validated with domain expertise and peer-reviewed analytical pipelines before drawing biological conclusions.

---

## Overview

The scRNA-seq Analysis Dashboard provides a point-and-click graphical interface for running a complete Seurat preprocessing and analysis workflow — from raw count matrices to annotated cell type populations — without writing a single line of R code. It is ideal for biologists, clinical researchers, and computational scientists who want to interactively explore their single-cell data.

---

## Features

- **Data Loading** — Load the built-in PBMC3k example or upload 10X Genomics directories (ZIP), RDS Seurat objects, or CSV count matrices
- **QC & Filtering** — Compute mitochondrial percentage, visualize quality distributions, and apply cell-level filters
- **Normalization** — LogNormalize / CLR / RC normalization, variable feature selection (VST / MVP / DISP)
- **PCA** — Run PCA, visualize elbow plots, loadings, and heatmaps to choose dimensionality
- **Clustering & UMAP** — Louvain/SLM clustering with adjustable resolution, UMAP and tSNE projections
- **Marker Genes** — FindAllMarkers with violin, feature, dot, heatmap, and ridge plot visualizations
- **Cell Type Annotation** — Assign custom labels to each cluster, export annotated Seurat objects
- **Interactive Vignette** — Step-by-step walkthrough built into the app

---

## Installation

### Prerequisites

- R >= 4.2.0
- RStudio (recommended) or command-line R

### Clone and install dependencies

```bash
git clone https://github.com/ShivaniPimparkar111/scrnaseq-dashboard.git
cd scrnaseq-dashboard
Rscript install_packages.R
```

The installation script will handle all CRAN and Bioconductor dependencies automatically.

---

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/ShivaniPimparkar111/scrnaseq-dashboard.git

# 2. Install dependencies
Rscript install_packages.R

# 3. Launch the app
R -e "shiny::runApp('.')"
```

The app will open in your default web browser at `http://127.0.0.1:PORT`.

Alternatively, open `app.R` in RStudio and click **Run App**.

---

## Workflow

The dashboard guides you through 8 sequential analysis steps:

| Step | Tab | Description |
|------|-----|-------------|
| 1 | **Load Data** | Import scRNA-seq count data (PBMC3k example or custom upload) |
| 2 | **QC & Filtering** | Compute QC metrics, visualize distributions, filter low-quality cells |
| 3 | **Normalization** | Normalize counts, identify highly variable genes, scale data |
| 4 | **PCA & Dimensionality** | Run PCA, inspect elbow plot to determine significant PCs |
| 5 | **Clustering & UMAP** | Build kNN graph, cluster cells, project to UMAP or tSNE |
| 6 | **Marker Genes** | Identify cluster-specific marker genes with multiple plot types |
| 7 | **Cell Type Annotation** | Label clusters with biological cell type names, export results |
| 8 | **Vignette** | Built-in step-by-step documentation |

---

## Data Formats Supported

| Format | Extension | Description |
|--------|-----------|-------------|
| 10X Genomics directory | `.zip` | ZIP containing `matrix.mtx`, `barcodes.tsv`, `features.tsv` |
| Seurat RDS object | `.rds` | Pre-processed or raw Seurat object |
| CSV count matrix | `.csv` | Rows = genes, columns = cells |
| Built-in example | — | PBMC3k dataset downloaded automatically from 10x Genomics |

---

## Screenshots

> Screenshots will be added after first public release.

```
[Load Data Tab]    [QC Violin Plots]    [UMAP Clusters]
[Marker Heatmap]   [Annotated UMAP]     [Composition Bar]
```

---

## Tech Stack

| Package | Version | Role |
|---------|---------|------|
| [Shiny](https://shiny.posit.co/) | >= 1.7 | Web application framework |
| [shinydashboard](https://rstudio.github.io/shinydashboard/) | >= 0.7 | Dashboard layout |
| [Seurat](https://satijalab.org/seurat/) | >= 5.0 | scRNA-seq analysis core |
| [ggplot2](https://ggplot2.tidyverse.org/) | >= 3.4 | Visualization |
| [patchwork](https://patchwork.data-imaginist.com/) | >= 1.1 | Plot composition |
| [DT](https://rstudio.github.io/DT/) | >= 0.25 | Interactive data tables |
| [shinyWidgets](https://dreamrs.github.io/shinyWidgets/) | >= 0.7 | Enhanced UI widgets |
| [shinycssloaders](https://github.com/daattali/shinycssloaders) | >= 1.0 | Loading spinners |
| [ggrepel](https://ggrepel.slowkow.com/) | >= 0.9 | Non-overlapping labels |

---

## File Structure

```
scrnaseq-dashboard/
├── app.R               # Main Shiny application
├── VIGNETTE.md         # Step-by-step user guide (displayed in Tab 8)
├── install_packages.R  # Dependency installation script
├── README.md           # This file
├── .gitignore
└── www/
    └── custom.css      # Custom dashboard styling
```

---

## Author

**Shivani Pimparkar**

- GitHub: [github.com/shivani-pimparkar](https://github.com/shivani-pimparkar)

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Acknowledgements

- [Satija Lab](https://satijalab.org/) for the Seurat package and PBMC3k tutorial
- [10x Genomics](https://www.10xgenomics.com/) for the publicly available PBMC3k dataset
- [Posit / RStudio](https://posit.co/) for the Shiny framework
