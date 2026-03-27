# Vignette: scRNA-seq Dashboard Walkthrough

**Author:** Shivani Pimparkar
**Last updated:** 2026-03-27

---

## Introduction

This dashboard guides you through the standard **Seurat preprocessing workflow** for single-cell RNA sequencing (scRNA-seq) data. Each tab corresponds to a distinct analysis step, and steps should generally be performed in order from top to bottom in the sidebar.

The workflow implemented here closely follows the [Seurat PBMC3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) and is suitable for:

- 10x Genomics Chromium data
- Drop-seq, inDrop, and other droplet-based platforms
- Any data that can be represented as a genes × cells count matrix

> **Important:** This dashboard is designed for initial data exploration only. Always validate findings with domain-specific knowledge and peer-reviewed methodology.

---

## Step 1 — Load Data

### What this step does biologically
Raw scRNA-seq data consists of a count matrix where each row is a gene and each column is a cell barcode. The counts represent how many RNA molecules from each gene were detected in each cell. Loading this matrix into a Seurat object is the first step of any analysis.

### Options

**PBMC3k Example**
Click "Load PBMC3k Example" to automatically download a classic benchmark dataset: ~2,700 peripheral blood mononuclear cells from a healthy donor, sequenced on the 10x Chromium platform. This is an excellent dataset to learn the workflow.

**Upload Your Own Data**

| Format | When to use |
|--------|-------------|
| 10X directory (ZIP) | You have raw output from Cell Ranger (`matrix.mtx`, `barcodes.tsv`, `features.tsv`) |
| RDS (Seurat object) | You already have a Seurat object created externally in R |
| CSV count matrix | You have a genes × cells table exported from another tool |

### What to look for
After loading, the **Dataset Summary** cards will show:
- **Total Cells** — for PBMC3k, expect ~2,700
- **Total Genes** — for PBMC3k, expect ~13,700 after `min.cells=3` filtering
- **Median Genes / Cell** — a rough indicator of data quality; typical values are 500–3,000 for 10x data
- **Median UMIs / Cell** — unique molecular identifiers per cell; 1,000–10,000 is typical

The **Count Matrix Preview** shows raw integer counts. Many values will be 0 (scRNA-seq data is sparse).

---

## Step 2 — QC & Filtering

### What this step does biologically
Not all barcodes captured represent live, healthy cells. Empty droplets, dead cells, and doublets (two cells in one droplet) need to be removed. Three key metrics help identify low-quality cells:

1. **nFeature_RNA** — number of distinct genes detected. Too few suggests an empty droplet or dead cell; too many suggests a doublet.
2. **nCount_RNA** — total UMI counts. Correlated with nFeature, but extreme outliers warrant scrutiny.
3. **percent.mt** — percentage of reads from mitochondrial genes. Dying cells lose cytoplasmic RNA but retain mitochondrial RNA, so high percent.mt indicates cell stress or death.

### Steps

1. Enter the **mitochondrial gene pattern** (`^MT-` for human, `^mt-` for mouse).
2. Click **Calculate QC Metrics** — this adds `percent.mt` to the metadata.
3. Inspect the violin and distribution plots.
4. Adjust the filter sliders:
   - `min_features`: remove empty droplets (default 200)
   - `max_features`: remove likely doublets (default 2,500)
   - `max_mt`: remove dying cells (default 5%)
5. Click **Apply Filters**.

### Typical values for PBMC3k
- Min features: 200
- Max features: 2,500
- Max percent.mt: 5%
- After filtering: ~2,638 cells remain

### What to look for in plots
- **Violin plots:** Healthy distributions are roughly log-normal. Look for a distinct low tail (empty droplets) and a long upper tail that may need trimming.
- **Scatter plots:** nCount vs nFeature should show a tight linear relationship. Outlier cells above the main diagonal may be doublets.
- **Distribution histograms:** The red dashed lines show your current filter thresholds. Adjust sliders until the lines fall at the natural inflection points of the distributions.

---

## Step 3 — Normalization & Variable Features

### What this step does biologically
Raw counts are biased by sequencing depth: cells sequenced more deeply will appear to express all genes at higher levels. Normalization corrects for this, making cells comparable. Next, we identify **highly variable genes (HVGs)** — genes that vary substantially across cells and are most likely to capture biologically meaningful differences. Using only HVGs for downstream analysis reduces noise and computation.

### Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| Normalization method | LogNormalize | Divide each cell by total counts × scale factor, then log-transform. Recommended for most datasets. |
| Scale factor | 10,000 | Effectively normalizes to counts per 10,000 (CP10K). |
| Variable feature method | vst | Variance-stabilizing transformation. Best for 10x data. |
| Number of variable features | 2,000 | 1,500–3,000 is typical. |
| Scale all genes | FALSE | If TRUE, scales all genes (needed for heatmaps showing non-HVGs). Slower. |

### What to look for
- **Variable Feature Plot:** HVGs appear in red above the grey non-variable genes. The top 10 genes are labeled. You should see a clear separation between highly variable (red) and low-variability (grey) genes.
- **Top Variable Genes table:** Genes like ribosomal protein genes, mitochondrial genes, and cell-cycle genes often appear at the top. Consider checking whether these are driving downstream clustering if results seem non-biological.

---

## Step 4 — PCA & Dimensionality

### What this step does biologically
PCA (Principal Component Analysis) compresses the ~2,000 variable genes into a smaller set of orthogonal axes that capture maximum variance. The first few PCs typically capture the major axes of biological variation (cell types), while later PCs capture noise. Choosing the right number of PCs to carry forward is critical.

### Steps

1. Set **Number of PCs to compute** (default 50).
2. Click **Run PCA**.
3. Inspect the **Elbow Plot** to determine where the variance explained per PC "elbows" off — the point after which additional PCs capture mostly noise.
4. Set the **PC Cutoff** slider to mark your chosen cutoff on the elbow plot.

### Choosing the number of PCs

The elbow plot shows standard deviation explained by each PC. Look for the point where the curve flattens. Common strategies:

- **Visual elbow:** Find where the drop in standard deviation becomes small and consistent.
- **Conservative approach:** Choose a slightly higher number than the obvious elbow (e.g., if the elbow is at PC 7, use 10).
- **For PBMC3k:** The elbow falls around PC 7–10. Using 10 PCs is standard.

### Sub-tabs explained

| Tab | What it shows |
|-----|---------------|
| Elbow Plot | Standard deviation per PC; red dashed line at your chosen cutoff |
| PCA Plot | Cells projected onto PC1 vs PC2 |
| Dim Loadings | Genes with highest positive/negative loading on each PC |
| Dim Heatmap (PC1) | Cells and genes sorted by PC1; useful for validating PC signal |
| Multi-PC Heatmap | Multiple PCs simultaneously; helps identify where signal fades |
| PC Gene Table | Full numeric loading table for all PCs |

---

## Step 5 — Clustering & UMAP

### What this step does biologically
Cells are clustered using a graph-based approach: first, a kNN (k-nearest neighbor) graph is built in PCA space, then a community detection algorithm (Louvain or SLM) partitions the graph into clusters. UMAP (Uniform Manifold Approximation and Projection) provides a 2D visualization that preserves local neighborhood structure.

### Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| Dimensions | 10 | Number of PCs to use. Should match your PC cutoff from Step 4. |
| K nearest neighbors | 20 | More neighbors → coarser clusters |
| Algorithm | Louvain | Standard choice; Louvain multi-level may give more stable results |
| Resolution | 0.5 | Higher → more clusters. For PBMC3k: 0.5 gives ~9 clusters |

### Recommended order
1. Click **Find Clusters** first.
2. Then click **Run Reduction** (UMAP or tSNE).

### Choosing resolution
- **Too low (< 0.3):** Under-clusters; different cell types merged together.
- **Too high (> 1.5):** Over-clusters; one cell type split into many clusters.
- **For PBMC3k:** 0.5 is well-validated and gives biologically meaningful clusters.

### What to look for
- **UMAP:** Well-separated clusters indicate distinct cell populations. Continuous gradients may represent cell states or trajectories.
- **Cluster QC violin plots:** All clusters should have similar QC metric distributions. A cluster with notably low features or high MT% may be low quality.
- **Split UMAP:** Use to check for batch effects or sample-specific clusters.

---

## Step 6 — Marker Genes

### What this step does biologically
Marker genes are genes that are significantly more highly expressed in one cluster compared to all others. Identifying markers allows you to match clusters to known cell types using published gene signatures (e.g., CD3D for T cells, CD79A for B cells, LYZ for monocytes).

### Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| Min pct expressing | 0.25 | Gene must be expressed in at least 25% of cells in the cluster |
| Log FC threshold | 0.25 | Minimum fold change. Higher = more stringent |
| Statistical test | Wilcoxon | Non-parametric; robust and recommended for most datasets |
| Positive only | TRUE | Only report genes upregulated in the cluster |

Click **Find All Markers** — this may take several minutes for large datasets.

### Plot types

| Plot | Best for |
|------|----------|
| **Marker Table** | Scanning markers, sorting by fold change or p-value |
| **Violin Plot** | Comparing expression distribution of a gene across clusters |
| **Feature Plot** | Seeing where a gene is expressed on the UMAP |
| **Dot Plot** | Visualizing multiple genes × multiple clusters simultaneously |
| **Heatmap** | Overview of top N markers per cluster |
| **Ridge Plot** | Distribution of expression per cluster, particularly for bimodal genes |

### Known PBMC3k marker genes

| Cluster | Cell Type | Key Markers |
|---------|-----------|-------------|
| 0, 4 | CD4+ T cells | IL7R, CCR7, S100A4 |
| 1 | CD14+ Monocytes | CD14, LYZ, CST3 |
| 2 | CD4+ T cells | IL7R, S100A4 |
| 3 | B cells | MS4A1, CD79A |
| 5 | CD8+ T cells | CD8A, GZMK |
| 6 | FCGR3A+ Monocytes | FCGR3A, MS4A7 |
| 7 | NK cells | GNLY, NKG7, GZMB |
| 8 | Dendritic cells | FCER1A, CST3 |

---

## Step 7 — Cell Type Annotation

### What this step does biologically
Based on the marker genes identified in Step 6, each cluster can be assigned a biologically meaningful label. This step replaces the numeric cluster IDs (0, 1, 2, ...) with cell type names.

### Steps

1. For each cluster, enter the cell type name in the corresponding text box. The boxes are pre-populated with "Cluster X" as defaults.
2. Use the marker gene results and known biology to determine appropriate labels.
3. Click **Apply Cell Type Labels**.
4. The **Annotated UMAP** and **Cell Type Composition** bar chart will update.

### Exporting results

- **Download Seurat Object (.rds):** Save the complete annotated Seurat object for further analysis in R.
- **Download Marker Table (.csv):** Export the FindAllMarkers results for supplementary materials or further analysis.

### What to look for
- The **Annotated UMAP** should show clusters with biologically coherent labels.
- The **Cell Type Composition** bar chart shows relative cell proportions. For PBMC3k, T cells should be the most abundant population (~50–60%).

---

## Tips & Troubleshooting

### General tips
- Follow the tabs in order (1 → 7). Some steps require previous steps to be completed.
- The app stores state in memory — if you refresh the browser, you will need to re-load your data.
- For large datasets (>50,000 cells), consider pre-processing in R and uploading an RDS file.

### Common issues

| Issue | Likely cause | Solution |
|-------|-------------|---------|
| No QC plots appear | `percent.mt` not calculated | Click "Calculate QC Metrics" before plotting |
| UMAP empty | Clustering not run | Complete Step 5 in order: clusters first, then UMAP |
| Marker finding is very slow | Large dataset | Reduce dataset size, or use `only.pos = TRUE` and increase `logfc.threshold` |
| Error: "object of type closure is not subsettable" | PCA not run | Complete Step 4 before clustering |
| Mitochondrial % all zeros | Wrong pattern | Use `^mt-` for mouse, `^MT-` for human |

---

## References

1. Stuart T, et al. (2019). Comprehensive Integration of Single-Cell Data. *Cell*, 177(7):1888–1902. https://doi.org/10.1016/j.cell.2019.05.031
2. Hao Y, et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*, 184(13):3573–3587. https://doi.org/10.1016/j.cell.2021.04.048
3. Seurat Documentation: https://satijalab.org/seurat/
4. PBMC3k Tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
5. Luecken MD & Theis FJ (2019). Current best practices in single-cell RNA-seq analysis. *Molecular Systems Biology*, 15:e8746. https://doi.org/10.15252/msb.20188746

---

*Author: Shivani Pimparkar*
*This dashboard is designed for initial exploration of scRNA-seq data. Findings should be validated before drawing biological conclusions.*
