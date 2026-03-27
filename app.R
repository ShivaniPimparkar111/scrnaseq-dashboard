# ==============================================================================
# scRNA-seq Analysis Dashboard
# Author: Shivani Pimparkar
# Description: Interactive Shiny dashboard for scRNA-seq data exploration using
#              the Seurat framework.
# ==============================================================================

library(shiny)
library(shinydashboard)
library(Seurat)
library(ggplot2)
library(dplyr)
library(DT)
library(shinycssloaders)
library(shinyWidgets)
library(patchwork)
library(ggrepel)
library(scales)
options(shiny.maxRequestSize = 500 * 1024^2)

# ==============================================================================
# UI
# ==============================================================================

ui <- dashboardPage(
  skin = "blue",

  # ---- Header ----
  dashboardHeader(
    title = tags$span(
      tags$img(src = "https://img.icons8.com/fluency/24/dna-helix.png",
               style = "margin-right:8px; vertical-align:middle;"),
      "scRNA-seq Dashboard"
    ),
    titleWidth = 280
  ),

  # ---- Sidebar ----
  dashboardSidebar(
    width = 280,
    tags$head(
      tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"),
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    sidebarMenu(
      id = "sidebar",
      menuItem("Load Data",            tabName = "data",      icon = icon("database")),
      menuItem("QC & Filtering",       tabName = "qc",        icon = icon("filter")),
      menuItem("Normalization",         tabName = "norm",      icon = icon("chart-line")),
      menuItem("PCA & Dimensionality", tabName = "pca",       icon = icon("project-diagram")),
      menuItem("Clustering & UMAP",    tabName = "cluster",   icon = icon("circle-nodes")),
      menuItem("Marker Genes",         tabName = "markers",   icon = icon("dna")),
      menuItem("Cell Type Annotation", tabName = "celltypes", icon = icon("tags")),
      menuItem("How to Use",           tabName = "vignette",  icon = icon("book-open"))
    ),
    tags$div(
      class = "sidebar-footer",
      tags$hr(style = "border-color: rgba(255,255,255,0.15); margin: 10px 15px;"),
      tags$p("Made by Shivani Pimparkar",
             style = "color:rgba(255,255,255,0.7); font-size:12px; margin:4px 15px; font-weight:500;"),
      tags$p("For initial data exploration only",
             style = "color:rgba(255,255,255,0.45); font-size:11px; margin:4px 15px;")
    )
  ),

  # ---- Body ----
  dashboardBody(
    tabItems(

      # ==========================================================================
      # TAB 1: LOAD DATA
      # ==========================================================================
      tabItem(
        tabName = "data",
        fluidRow(
          column(12,
            tags$h2("Load Data", class = "tab-title"),
            tags$p("Load the PBMC3k example dataset or upload your own scRNA-seq data.",
                   class = "tab-subtitle")
          )
        ),
        fluidRow(
          # PBMC3k Example
          box(
            title = tags$span(icon("cloud-download-alt"), " Example Dataset"),
            width = 4, solidHeader = TRUE, status = "primary",
            tags$p("Load the classic PBMC3k dataset from 10x Genomics (~60 MB download)."),
            tags$p(tags$strong("2,700 peripheral blood mononuclear cells"), " from a healthy donor."),
            tags$br(),
            actionButton("load_pbmc3k", "Load PBMC3k Example",
                         icon = icon("play-circle"),
                         class = "btn-primary btn-lg btn-block")
          ),
          # Upload
          box(
            title = tags$span(icon("upload"), " Upload Your Data"),
            width = 8, solidHeader = TRUE, status = "info",
            radioButtons("upload_type", "Data Format:",
                         choices = c(
                           "10X directory (ZIP)" = "10x_zip",
                           "RDS (Seurat object)"  = "rds",
                           "CSV count matrix"     = "csv"
                         ),
                         inline = TRUE),
            conditionalPanel(
              condition = "input.upload_type == '10x_zip'",
              tags$p(tags$em("ZIP file should contain: barcodes.tsv, genes.tsv (or features.tsv), matrix.mtx"),
                     style = "font-size:12px; color:#666;")
            ),
            conditionalPanel(
              condition = "input.upload_type == 'csv'",
              tags$p(tags$em("CSV: rows = genes, columns = cells, first column = gene names"),
                     style = "font-size:12px; color:#666;")
            ),
            fileInput("user_file", "Choose File:",
                      accept = c(".rds", ".RDS", ".csv", ".CSV", ".zip", ".ZIP")),
            actionButton("load_uploaded", "Load Uploaded Data",
                         icon = icon("check-circle"),
                         class = "btn-success")
          )
        ),
        fluidRow(
          box(
            title = tags$span(icon("info-circle"), " Dataset Summary"),
            width = 12, solidHeader = TRUE, status = "primary",
            uiOutput("data_summary")
          )
        ),
        fluidRow(
          box(
            title = tags$span(icon("table"), " Count Matrix Preview (top 5 genes × 5 cells)"),
            width = 6, solidHeader = TRUE, status = "info",
            withSpinner(DTOutput("count_preview"), type = 6, color = "#2196F3")
          ),
          box(
            title = tags$span(icon("list"), " Cell Metadata"),
            width = 6, solidHeader = TRUE, status = "info",
            withSpinner(DTOutput("meta_preview"), type = 6, color = "#2196F3")
          )
        )
      ),

      # ==========================================================================
      # TAB 2: QC & FILTERING
      # ==========================================================================
      tabItem(
        tabName = "qc",
        fluidRow(
          column(12,
            tags$h2("QC & Filtering", class = "tab-title"),
            tags$p("Calculate quality control metrics and filter low-quality cells.", class = "tab-subtitle")
          )
        ),
        fluidRow(
          box(
            title = tags$span(icon("sliders-h"), " QC Parameters"),
            width = 4, solidHeader = TRUE, status = "warning",
            textInput("mt_pattern", "Mitochondrial Gene Pattern:", value = "^MT-"),
            tags$hr(),
            sliderInput("min_features", "Min Features (genes/cell):",
                        min = 0, max = 1000, value = 200, step = 50),
            sliderInput("max_features", "Max Features (genes/cell):",
                        min = 500, max = 10000, value = 2500, step = 100),
            sliderInput("max_mt", "Max Mitochondrial % :",
                        min = 1, max = 50, value = 5, step = 1),
            tags$hr(),
            actionButton("apply_qc", "Calculate QC Metrics",
                         icon = icon("calculator"), class = "btn-warning btn-block"),
            tags$br(),
            actionButton("apply_filter", "Apply Filters",
                         icon = icon("filter"), class = "btn-danger btn-block"),
            tags$br(),
            uiOutput("filter_summary")
          ),
          box(
            title = tags$span(icon("chart-area"), " QC Visualizations"),
            width = 8, solidHeader = TRUE, status = "primary",
            tabsetPanel(
              tabPanel("Violin Plots",
                       withSpinner(plotOutput("vln_qc", height = "400px"), type = 6, color = "#2196F3")),
              tabPanel("Scatter Plots",
                       withSpinner(plotOutput("scatter_qc", height = "400px"), type = 6, color = "#2196F3")),
              tabPanel("Distributions",
                       withSpinner(plotOutput("hist_qc", height = "400px"), type = 6, color = "#2196F3"))
            )
          )
        )
      ),

      # ==========================================================================
      # TAB 3: NORMALIZATION
      # ==========================================================================
      tabItem(
        tabName = "norm",
        fluidRow(
          column(12,
            tags$h2("Normalization & Variable Features", class = "tab-title"),
            tags$p("Normalize expression values and identify highly variable genes.", class = "tab-subtitle")
          )
        ),
        fluidRow(
          box(
            title = tags$span(icon("cog"), " Normalization Parameters"),
            width = 4, solidHeader = TRUE, status = "warning",
            selectInput("norm_method", "Normalization Method:",
                        choices = c("LogNormalize", "CLR", "RC"),
                        selected = "LogNormalize"),
            numericInput("scale_factor", "Scale Factor:", value = 10000, min = 1000, max = 100000),
            tags$hr(),
            selectInput("var_method", "Variable Feature Method:",
                        choices = c("vst", "mvp", "disp"),
                        selected = "vst"),
            sliderInput("n_features", "Number of Variable Features:",
                        min = 500, max = 5000, value = 2000, step = 100),
            tags$hr(),
            checkboxInput("scale_all_genes", "Scale All Genes (slower)", value = FALSE),
            tags$br(),
            actionButton("run_norm", "Run Normalization",
                         icon = icon("play"), class = "btn-success btn-block")
          ),
          box(
            title = tags$span(icon("chart-line"), " Normalization Results"),
            width = 8, solidHeader = TRUE, status = "primary",
            tabsetPanel(
              tabPanel("Variable Feature Plot",
                       withSpinner(plotOutput("var_feat_plot", height = "450px"), type = 6, color = "#2196F3")),
              tabPanel("Top Variable Genes",
                       withSpinner(DTOutput("var_genes_table"), type = 6, color = "#2196F3"))
            )
          )
        )
      ),

      # ==========================================================================
      # TAB 4: PCA & DIMENSIONALITY
      # ==========================================================================
      tabItem(
        tabName = "pca",
        fluidRow(
          column(12,
            tags$h2("PCA & Dimensionality Reduction", class = "tab-title"),
            tags$p("Run PCA and determine the number of significant principal components.", class = "tab-subtitle")
          )
        ),
        fluidRow(
          box(
            title = tags$span(icon("cog"), " PCA Parameters"),
            width = 3, solidHeader = TRUE, status = "warning",
            sliderInput("n_pcs", "Number of PCs to Compute:",
                        min = 10, max = 100, value = 50, step = 5),
            actionButton("run_pca", "Run PCA",
                         icon = icon("play"), class = "btn-success btn-block"),
            tags$hr(),
            tags$strong("Visualization Controls"),
            sliderInput("pc_cutoff", "PC Cutoff (Elbow):",
                        min = 2, max = 50, value = 10, step = 1),
            sliderInput("viz_dims", "Dims to Visualize (Loadings):",
                        min = 1, max = 20, value = c(1, 2), step = 1),
            sliderInput("heatmap_dims", "Dims for Multi-PC Heatmap:",
                        min = 1, max = 20, value = c(1, 15), step = 1),
            sliderInput("heatmap_cells", "Cells for Heatmap:",
                        min = 100, max = 1000, value = 500, step = 100)
          ),
          box(
            title = tags$span(icon("chart-bar"), " PCA Results"),
            width = 9, solidHeader = TRUE, status = "primary",
            tabsetPanel(
              tabPanel("Elbow Plot",
                       withSpinner(plotOutput("elbow_plot", height = "420px"), type = 6, color = "#2196F3")),
              tabPanel("PCA Plot",
                       withSpinner(plotOutput("pca_plot", height = "420px"), type = 6, color = "#2196F3")),
              tabPanel("Dim Loadings",
                       withSpinner(plotOutput("viz_dim_plot", height = "420px"), type = 6, color = "#2196F3")),
              tabPanel("Dim Heatmap (PC1)",
                       withSpinner(plotOutput("dim_heatmap_single", height = "420px"), type = 6, color = "#2196F3")),
              tabPanel("Multi-PC Heatmap",
                       withSpinner(plotOutput("dim_heatmap_multi", height = "600px"), type = 6, color = "#2196F3")),
              tabPanel("PC Gene Table",
                       withSpinner(DTOutput("pca_gene_table"), type = 6, color = "#2196F3"))
            )
          )
        )
      ),

      # ==========================================================================
      # TAB 5: CLUSTERING & UMAP
      # ==========================================================================
      tabItem(
        tabName = "cluster",
        fluidRow(
          column(12,
            tags$h2("Clustering & UMAP", class = "tab-title"),
            tags$p("Identify cell clusters and visualize in low-dimensional space.", class = "tab-subtitle")
          )
        ),
        fluidRow(
          box(
            title = tags$span(icon("cog"), " Clustering Parameters"),
            width = 3, solidHeader = TRUE, status = "warning",
            sliderInput("cluster_dims", "Dimensions for Clustering:",
                        min = 2, max = 50, value = 10, step = 1),
            numericInput("k_param", "K Nearest Neighbors:", value = 20, min = 5, max = 100),
            selectInput("cluster_algo", "Clustering Algorithm:",
                        choices = c("Louvain" = "1",
                                    "Louvain multi-level" = "2",
                                    "SLM" = "3"),
                        selected = "1"),
            sliderInput("resolution", "Resolution:",
                        min = 0.1, max = 2.0, value = 0.5, step = 0.1),
            actionButton("run_cluster", "Find Clusters",
                         icon = icon("circle-nodes"), class = "btn-danger btn-block"),
            tags$hr(),
            selectInput("reduction_type", "Reduction Method:",
                        choices = c("UMAP" = "umap", "tSNE" = "tsne")),
            actionButton("run_umap", "Run Reduction",
                         icon = icon("play"), class = "btn-success btn-block"),
            tags$hr(),
            checkboxInput("show_labels", "Show Cluster Labels", value = TRUE),
            uiOutput("color_by_ui")
          ),
          box(
            title = tags$span(icon("circle-nodes"), " Clustering Results"),
            width = 9, solidHeader = TRUE, status = "primary",
            tabsetPanel(
              tabPanel("UMAP / tSNE",
                       withSpinner(plotOutput("umap_plot", height = "500px"), type = 6, color = "#2196F3")),
              tabPanel("Cluster QC",
                       withSpinner(plotOutput("cluster_vln", height = "500px"), type = 6, color = "#2196F3")),
              tabPanel("Split UMAP",
                       uiOutput("split_by_ui"),
                       withSpinner(plotOutput("umap_split", height = "500px"), type = 6, color = "#2196F3"))
            )
          )
        )
      ),

      # ==========================================================================
      # TAB 6: MARKER GENES
      # ==========================================================================
      tabItem(
        tabName = "markers",
        fluidRow(
          column(12,
            tags$h2("Marker Genes", class = "tab-title"),
            tags$p("Identify differentially expressed marker genes for each cluster.", class = "tab-subtitle")
          )
        ),
        fluidRow(
          box(
            title = tags$span(icon("cog"), " Marker Detection Parameters"),
            width = 3, solidHeader = TRUE, status = "warning",
            sliderInput("min_pct", "Min Pct Expressing:",
                        min = 0.01, max = 0.5, value = 0.25, step = 0.01),
            sliderInput("logfc_threshold", "Log FC Threshold:",
                        min = 0.1, max = 2.0, value = 0.25, step = 0.05),
            selectInput("de_test", "Statistical Test:",
                        choices = c("Wilcoxon" = "wilcox",
                                    "Bimodal" = "bimod",
                                    "ROC" = "roc",
                                    "Student's t" = "t"),
                        selected = "wilcox"),
            checkboxInput("only_pos", "Positive Markers Only", value = TRUE),
            actionButton("run_markers", "Find All Markers",
                         icon = icon("search"), class = "btn-primary btn-block"),
            tags$hr(),
            tags$strong("Visualization"),
            uiOutput("gene_picker_ui"),
            numericInput("top_n_heatmap", "Top N Genes for Heatmap:", value = 5, min = 1, max = 20)
          ),
          box(
            title = tags$span(icon("dna"), " Marker Gene Results"),
            width = 9, solidHeader = TRUE, status = "primary",
            tabsetPanel(
              tabPanel("Marker Table",
                       withSpinner(DTOutput("marker_table"), type = 6, color = "#2196F3")),
              tabPanel("Violin Plot",
                       withSpinner(plotOutput("marker_vln", height = "480px"), type = 6, color = "#2196F3")),
              tabPanel("Feature Plot",
                       withSpinner(plotOutput("marker_feature", height = "480px"), type = 6, color = "#2196F3")),
              tabPanel("Dot Plot",
                       withSpinner(plotOutput("marker_dot", height = "480px"), type = 6, color = "#2196F3")),
              tabPanel("Heatmap",
                       withSpinner(plotOutput("marker_heatmap", height = "600px"), type = 6, color = "#2196F3")),
              tabPanel("Ridge Plot",
                       withSpinner(plotOutput("marker_ridge", height = "480px"), type = 6, color = "#2196F3"))
            )
          )
        )
      ),

      # ==========================================================================
      # TAB 7: CELL TYPE ANNOTATION
      # ==========================================================================
      tabItem(
        tabName = "celltypes",
        fluidRow(
          column(12,
            tags$h2("Cell Type Annotation", class = "tab-title"),
            tags$p("Assign biological cell type labels to each cluster.", class = "tab-subtitle")
          )
        ),
        fluidRow(
          box(
            title = tags$span(icon("tags"), " Assign Cell Type Labels"),
            width = 4, solidHeader = TRUE, status = "warning",
            tags$p("Enter a cell type label for each cluster. Pre-populated with cluster numbers."),
            uiOutput("celltype_inputs"),
            tags$br(),
            actionButton("apply_celltypes", "Apply Cell Type Labels",
                         icon = icon("check"), class = "btn-success btn-block"),
            tags$hr(),
            tags$strong("Export Results"),
            tags$br(), tags$br(),
            downloadButton("download_rds", "Download Seurat Object (.rds)",
                           class = "btn-primary btn-block"),
            tags$br(),
            downloadButton("download_markers", "Download Marker Table (.csv)",
                           class = "btn-info btn-block")
          ),
          box(
            title = tags$span(icon("palette"), " Annotated Visualizations"),
            width = 8, solidHeader = TRUE, status = "primary",
            tabsetPanel(
              tabPanel("Annotated UMAP",
                       withSpinner(plotOutput("annotated_umap", height = "500px"), type = 6, color = "#2196F3")),
              tabPanel("Cell Type Composition",
                       withSpinner(plotOutput("celltype_bar", height = "500px"), type = 6, color = "#2196F3"))
            )
          )
        )
      ),

      # ==========================================================================
      # TAB 8: VIGNETTE
      # ==========================================================================
      tabItem(
        tabName = "vignette",
        fluidRow(
          box(
            title = tags$span(icon("book-open"), " How to Use This Dashboard"),
            width = 12, solidHeader = TRUE, status = "info",
            includeMarkdown("VIGNETTE.md")
          )
        )
      )

    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage


# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {

  # Reactive values store
  rv <- reactiveValues(
    seurat    = NULL,
    filtered  = FALSE,
    markers   = NULL,
    step      = "none"   # tracks furthest completed step
  )

  # --------------------------------------------------------------------------
  # TAB 1: LOAD DATA
  # --------------------------------------------------------------------------

  # Load PBMC3k example
  observeEvent(input$load_pbmc3k, {
    withProgress(message = "Downloading PBMC3k data...", value = 0, {
      tryCatch({
        setProgress(0.1, detail = "Downloading from 10x Genomics...")
        tmp_dir  <- tempdir()
        tmp_file <- file.path(tmp_dir, "pbmc3k.tar.gz")
        url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
        download.file(url, tmp_file, mode = "wb", quiet = TRUE)

        setProgress(0.4, detail = "Extracting archive...")
        untar(tmp_file, exdir = tmp_dir)

        setProgress(0.6, detail = "Reading 10X data...")
        data_path <- file.path(tmp_dir, "filtered_gene_bc_matrices", "hg19")
        pbmc_data <- Read10X(data.dir = data_path)

        setProgress(0.8, detail = "Creating Seurat object...")
        rv$seurat   <- CreateSeuratObject(counts = pbmc_data,
                                          project = "PBMC3k",
                                          min.cells = 3,
                                          min.features = 200)
        rv$filtered <- FALSE
        rv$step     <- "loaded"

        setProgress(1.0, detail = "Done!")
        showNotification("PBMC3k dataset loaded successfully!", type = "message", duration = 5)
      }, error = function(e) {
        showNotification(paste("Error loading PBMC3k:", conditionMessage(e)),
                         type = "error", duration = 10)
      })
    })
  })

  # Load uploaded file
  observeEvent(input$load_uploaded, {
    req(input$user_file)
    withProgress(message = "Loading uploaded data...", value = 0, {
      tryCatch({
        path  <- input$user_file$datapath
        ftype <- input$upload_type

        if (ftype == "rds") {
          setProgress(0.5, detail = "Reading RDS file...")
          obj <- readRDS(path)
          if (!inherits(obj, "Seurat"))
            stop("RDS file does not contain a Seurat object.")
          rv$seurat <- obj

        } else if (ftype == "csv") {
          setProgress(0.3, detail = "Reading CSV file...")
          mat <- read.csv(path, row.names = 1, check.names = FALSE)
          setProgress(0.6, detail = "Creating Seurat object...")
          rv$seurat <- CreateSeuratObject(counts = as.matrix(mat),
                                          project = "Uploaded",
                                          min.cells = 3,
                                          min.features = 200)

        } else if (ftype == "10x_zip") {
          setProgress(0.2, detail = "Extracting ZIP...")
          tmp_dir <- tempfile()
          dir.create(tmp_dir)
          unzip(path, exdir = tmp_dir)
          # find the directory containing matrix.mtx
          mtx_file <- list.files(tmp_dir, pattern = "matrix.mtx", recursive = TRUE, full.names = TRUE)[1]
          if (is.na(mtx_file)) stop("No matrix.mtx found in ZIP.")
          data_path <- dirname(mtx_file)
          setProgress(0.6, detail = "Reading 10X data...")
          counts <- Read10X(data.dir = data_path)
          setProgress(0.8, detail = "Creating Seurat object...")
          rv$seurat <- CreateSeuratObject(counts = counts,
                                          project = "Uploaded10X",
                                          min.cells = 3,
                                          min.features = 200)
        }

        rv$filtered <- FALSE
        rv$step     <- "loaded"
        setProgress(1.0, detail = "Done!")
        showNotification("Data loaded successfully!", type = "message", duration = 5)

      }, error = function(e) {
        showNotification(paste("Error loading data:", conditionMessage(e)),
                         type = "error", duration = 10)
      })
    })
  })

  # Dataset summary
  output$data_summary <- renderUI({
    req(rv$seurat)
    s <- rv$seurat
    n_cells  <- ncol(s)
    n_genes  <- nrow(s)
    project  <- s@project.name
    med_feat <- median(s$nFeature_RNA)
    med_cnt  <- median(s$nCount_RNA)

    tags$div(
      fluidRow(
        valueBox(formatC(n_cells, big.mark = ","), "Total Cells",
                 icon = icon("circle"), color = "blue", width = 3),
        valueBox(formatC(n_genes, big.mark = ","), "Total Genes",
                 icon = icon("dna"), color = "green", width = 3),
        valueBox(formatC(med_feat, big.mark = ","), "Median Genes / Cell",
                 icon = icon("chart-bar"), color = "purple", width = 3),
        valueBox(formatC(med_cnt, big.mark = ","), "Median UMIs / Cell",
                 icon = icon("calculator"), color = "yellow", width = 3)
      ),
      tags$p(tags$strong("Project: "), project,
             tags$span(" | "), tags$strong("Assay: "), DefaultAssay(s),
             style = "margin-top:5px; color:#555; font-size:13px;")
    )
  })

  # Count matrix preview
  output$count_preview <- renderDT({
    req(rv$seurat)
    mat <- as.matrix(GetAssayData(rv$seurat, slot = "counts"))[1:5, 1:5]
    datatable(as.data.frame(mat),
              options = list(dom = "t", scrollX = TRUE, pageLength = 5),
              class = "compact stripe")
  })

  # Metadata preview
  output$meta_preview <- renderDT({
    req(rv$seurat)
    datatable(rv$seurat@meta.data,
              options = list(pageLength = 10, scrollX = TRUE),
              class = "compact stripe")
  })

  # --------------------------------------------------------------------------
  # TAB 2: QC & FILTERING
  # --------------------------------------------------------------------------

  observeEvent(input$apply_qc, {
    req(rv$seurat)
    tryCatch({
      rv$seurat[["percent.mt"]] <- PercentageFeatureSet(rv$seurat,
                                                        pattern = input$mt_pattern)
      showNotification("QC metrics calculated.", type = "message", duration = 3)
    }, error = function(e) {
      showNotification(paste("QC Error:", conditionMessage(e)), type = "error", duration = 8)
    })
  })

  observeEvent(input$apply_filter, {
    req(rv$seurat, "percent.mt" %in% colnames(rv$seurat@meta.data))
    tryCatch({
      cells_before <- ncol(rv$seurat)
      rv$seurat <- subset(rv$seurat,
                          subset = nFeature_RNA >= input$min_features &
                            nFeature_RNA <= input$max_features &
                            percent.mt  <= input$max_mt)
      cells_after  <- ncol(rv$seurat)
      rv$filtered  <- TRUE
      rv$step      <- "qc"
      showNotification(
        paste0("Filtered: ", cells_before - cells_after, " cells removed. ",
               cells_after, " cells remaining."),
        type = "message", duration = 5
      )
    }, error = function(e) {
      showNotification(paste("Filter Error:", conditionMessage(e)), type = "error", duration = 8)
    })
  })

  output$filter_summary <- renderUI({
    req(rv$seurat, rv$filtered)
    tags$div(
      class = "filter-summary-box",
      tags$p(tags$strong("Cells remaining: "),
             formatC(ncol(rv$seurat), big.mark = ","),
             style = "color:#27ae60; font-size:14px;")
    )
  })

  output$vln_qc <- renderPlot({
    req(rv$seurat, "percent.mt" %in% colnames(rv$seurat@meta.data))
    VlnPlot(rv$seurat,
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
            ncol = 3, pt.size = 0.1) &
      theme_minimal(base_size = 14)
  })

  output$scatter_qc <- renderPlot({
    req(rv$seurat, "percent.mt" %in% colnames(rv$seurat@meta.data))
    p1 <- FeatureScatter(rv$seurat, "nCount_RNA", "nFeature_RNA") +
      theme_minimal(base_size = 14)
    p2 <- FeatureScatter(rv$seurat, "nCount_RNA", "percent.mt") +
      theme_minimal(base_size = 14)
    p1 + p2
  })

  output$hist_qc <- renderPlot({
    req(rv$seurat, "percent.mt" %in% colnames(rv$seurat@meta.data))
    meta <- rv$seurat@meta.data
    p1 <- ggplot(meta, aes(x = nFeature_RNA)) +
      geom_histogram(bins = 60, fill = "#2196F3", color = "white", alpha = 0.8) +
      geom_vline(xintercept = input$min_features, color = "red", linetype = "dashed", linewidth = 1) +
      geom_vline(xintercept = input$max_features, color = "red", linetype = "dashed", linewidth = 1) +
      labs(title = "Genes per Cell", x = "nFeature_RNA", y = "Count") +
      theme_minimal(base_size = 14)
    p2 <- ggplot(meta, aes(x = nCount_RNA)) +
      geom_histogram(bins = 60, fill = "#4CAF50", color = "white", alpha = 0.8) +
      labs(title = "UMIs per Cell", x = "nCount_RNA", y = "Count") +
      theme_minimal(base_size = 14)
    p3 <- ggplot(meta, aes(x = percent.mt)) +
      geom_histogram(bins = 60, fill = "#FF5722", color = "white", alpha = 0.8) +
      geom_vline(xintercept = input$max_mt, color = "red", linetype = "dashed", linewidth = 1) +
      labs(title = "Mitochondrial %", x = "percent.mt", y = "Count") +
      theme_minimal(base_size = 14)
    p1 + p2 + p3 + plot_layout(ncol = 3)
  })

  # --------------------------------------------------------------------------
  # TAB 3: NORMALIZATION
  # --------------------------------------------------------------------------

  observeEvent(input$run_norm, {
    req(rv$seurat)
    withProgress(message = "Running normalization pipeline...", value = 0, {
      tryCatch({
        setProgress(0.2, detail = "Normalizing data...")
        rv$seurat <- NormalizeData(rv$seurat,
                                   normalization.method = input$norm_method,
                                   scale.factor = input$scale_factor)
        setProgress(0.4, detail = "Finding variable features...")
        rv$seurat <- FindVariableFeatures(rv$seurat,
                                          selection.method = input$var_method,
                                          nfeatures = input$n_features)
        setProgress(0.7, detail = "Scaling data...")
        genes_to_scale <- if (input$scale_all_genes) rownames(rv$seurat) else NULL
        rv$seurat <- ScaleData(rv$seurat, features = genes_to_scale)
        rv$step   <- "norm"
        setProgress(1.0, detail = "Done!")
        showNotification("Normalization complete!", type = "message", duration = 4)
      }, error = function(e) {
        showNotification(paste("Normalization Error:", conditionMessage(e)),
                         type = "error", duration = 8)
      })
    })
  })

  output$var_feat_plot <- renderPlot({
    req(rv$seurat, length(VariableFeatures(rv$seurat)) > 0)
    top_genes <- head(VariableFeatures(rv$seurat), 10)
    p <- VariableFeaturePlot(rv$seurat) + theme_minimal(base_size = 14)
    LabelPoints(plot = p, points = top_genes, repel = TRUE, xnudge = 0, ynudge = 0)
  })

  output$var_genes_table <- renderDT({
    req(rv$seurat, length(VariableFeatures(rv$seurat)) > 0)
    hvg <- HVFInfo(rv$seurat)
    hvg$gene   <- rownames(hvg)
    hvg$variable <- rownames(hvg) %in% VariableFeatures(rv$seurat)
    hvg <- hvg[order(-hvg$variable), ]
    datatable(hvg,
              options = list(pageLength = 15, scrollX = TRUE),
              class = "compact stripe")
  })

  # --------------------------------------------------------------------------
  # TAB 4: PCA
  # --------------------------------------------------------------------------

  observeEvent(input$run_pca, {
    req(rv$seurat)
    withProgress(message = "Running PCA...", value = 0.3, {
      tryCatch({
        rv$seurat <- RunPCA(rv$seurat, npcs = input$n_pcs)
        rv$step   <- "pca"
        setProgress(1.0)
        showNotification("PCA complete!", type = "message", duration = 3)
      }, error = function(e) {
        showNotification(paste("PCA Error:", conditionMessage(e)), type = "error", duration = 8)
      })
    })
  })

  output$elbow_plot <- renderPlot({
    req(rv$seurat, "pca" %in% names(rv$seurat@reductions))
    ElbowPlot(rv$seurat, ndims = rv$seurat@reductions$pca@cell.embeddings %>% ncol()) +
      geom_vline(xintercept = input$pc_cutoff, color = "red",
                 linetype = "dashed", linewidth = 1) +
      annotate("text", x = input$pc_cutoff + 0.5, y = Inf,
               label = paste0("PC ", input$pc_cutoff),
               color = "red", hjust = 0, vjust = 1.5, size = 4) +
      theme_minimal(base_size = 14) +
      labs(title = "Elbow Plot — PCA Standard Deviations")
  })

  output$pca_plot <- renderPlot({
    req(rv$seurat, "pca" %in% names(rv$seurat@reductions))
    DimPlot(rv$seurat, reduction = "pca") +
      theme_minimal(base_size = 14) +
      labs(title = "PCA Plot")
  })

  output$viz_dim_plot <- renderPlot({
    req(rv$seurat, "pca" %in% names(rv$seurat@reductions))
    dims <- seq(input$viz_dims[1], input$viz_dims[2])
    VizDimLoadings(rv$seurat, dims = dims, reduction = "pca") &
      theme_minimal(base_size = 12)
  })

  output$dim_heatmap_single <- renderPlot({
    req(rv$seurat, "pca" %in% names(rv$seurat@reductions))
    DimHeatmap(rv$seurat, dims = 1, cells = input$heatmap_cells, balanced = TRUE)
  })

  output$dim_heatmap_multi <- renderPlot({
    req(rv$seurat, "pca" %in% names(rv$seurat@reductions))
    dims <- seq(input$heatmap_dims[1], input$heatmap_dims[2])
    DimHeatmap(rv$seurat, dims = dims, cells = input$heatmap_cells, balanced = TRUE)
  })

  output$pca_gene_table <- renderDT({
    req(rv$seurat, "pca" %in% names(rv$seurat@reductions))
    loadings <- as.data.frame(rv$seurat@reductions$pca@feature.loadings)
    loadings$gene <- rownames(loadings)
    loadings <- loadings[, c("gene", colnames(loadings)[colnames(loadings) != "gene"])]
    datatable(loadings,
              options = list(pageLength = 15, scrollX = TRUE),
              class = "compact stripe") %>%
      formatRound(columns = 2:ncol(loadings), digits = 4)
  })

  # --------------------------------------------------------------------------
  # TAB 5: CLUSTERING & UMAP
  # --------------------------------------------------------------------------

  observeEvent(input$run_cluster, {
    req(rv$seurat, "pca" %in% names(rv$seurat@reductions))
    withProgress(message = "Finding clusters...", value = 0, {
      tryCatch({
        setProgress(0.3, detail = "Building neighbor graph...")
        rv$seurat <- FindNeighbors(rv$seurat,
                                   dims = 1:input$cluster_dims,
                                   k.param = input$k_param)
        setProgress(0.7, detail = "Finding clusters...")
        rv$seurat <- FindClusters(rv$seurat,
                                  resolution = input$resolution,
                                  algorithm = as.numeric(input$cluster_algo))
        rv$step <- "clustered"

        # Update color_by choices
        meta_cols <- colnames(rv$seurat@meta.data)
        updateSelectInput(session, "color_by",
                          choices = meta_cols,
                          selected = "seurat_clusters")

        setProgress(1.0)
        n_clusters <- length(levels(rv$seurat$seurat_clusters))
        showNotification(paste(n_clusters, "clusters found!"), type = "message", duration = 4)
      }, error = function(e) {
        showNotification(paste("Clustering Error:", conditionMessage(e)), type = "error", duration = 8)
      })
    })
  })

  observeEvent(input$run_umap, {
    req(rv$seurat, "pca" %in% names(rv$seurat@reductions))
    withProgress(message = paste("Running", toupper(input$reduction_type), "..."), value = 0.3, {
      tryCatch({
        if (input$reduction_type == "umap") {
          rv$seurat <- RunUMAP(rv$seurat, dims = 1:input$cluster_dims)
        } else {
          rv$seurat <- RunTSNE(rv$seurat, dims = 1:input$cluster_dims)
        }
        rv$step <- "umap"
        setProgress(1.0)
        showNotification(paste(toupper(input$reduction_type), "complete!"),
                         type = "message", duration = 3)
      }, error = function(e) {
        showNotification(paste("Reduction Error:", conditionMessage(e)), type = "error", duration = 8)
      })
    })
  })

  output$color_by_ui <- renderUI({
    req(rv$seurat)
    meta_cols <- colnames(rv$seurat@meta.data)
    selected  <- if ("seurat_clusters" %in% meta_cols) "seurat_clusters" else meta_cols[1]
    selectInput("color_by", "Color By:", choices = meta_cols, selected = selected)
  })

  output$umap_plot <- renderPlot({
    req(rv$seurat)
    red <- input$reduction_type
    req(red %in% names(rv$seurat@reductions))
    p <- DimPlot(rv$seurat,
                 reduction = red,
                 group.by  = input$color_by,
                 label     = input$show_labels,
                 repel     = TRUE) +
      theme_minimal(base_size = 14) +
      labs(title = paste(toupper(red), "—", input$color_by))
    p
  })

  output$cluster_vln <- renderPlot({
    req(rv$seurat, "seurat_clusters" %in% colnames(rv$seurat@meta.data),
        "percent.mt" %in% colnames(rv$seurat@meta.data))
    VlnPlot(rv$seurat,
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
            ncol = 3, pt.size = 0) &
      theme_minimal(base_size = 12)
  })

  output$split_by_ui <- renderUI({
    req(rv$seurat)
    meta_cols <- colnames(rv$seurat@meta.data)
    selectInput("split_by_var", "Split By:",
                choices = meta_cols, selected = meta_cols[1])
  })

  output$umap_split <- renderPlot({
    req(rv$seurat, input$split_by_var)
    red <- input$reduction_type
    req(red %in% names(rv$seurat@reductions))
    DimPlot(rv$seurat,
            reduction = red,
            split.by  = input$split_by_var,
            label     = input$show_labels) +
      theme_minimal(base_size = 13)
  })

  # --------------------------------------------------------------------------
  # TAB 6: MARKER GENES
  # --------------------------------------------------------------------------

  observeEvent(input$run_markers, {
    req(rv$seurat, "seurat_clusters" %in% colnames(rv$seurat@meta.data))
    withProgress(message = "Finding marker genes...", value = 0.2, {
      tryCatch({
        Idents(rv$seurat) <- "seurat_clusters"
        rv$markers <- FindAllMarkers(rv$seurat,
                                     only.pos      = input$only_pos,
                                     min.pct       = input$min_pct,
                                     logfc.threshold = input$logfc_threshold,
                                     test.use      = input$de_test)
        setProgress(0.9)
        # Populate gene picker
        top_genes <- rv$markers %>%
          group_by(cluster) %>%
          slice_max(order_by = avg_log2FC, n = 5) %>%
          pull(gene) %>%
          unique() %>%
          as.character()
        updatePickerInput(session, "selected_genes", choices = top_genes,
                          selected = top_genes[1:min(6, length(top_genes))])
        setProgress(1.0)
        showNotification(paste(nrow(rv$markers), "markers found!"),
                         type = "message", duration = 4)
      }, error = function(e) {
        showNotification(paste("Marker Error:", conditionMessage(e)), type = "error", duration = 8)
      })
    })
  })

  output$gene_picker_ui <- renderUI({
    req(rv$markers)
    all_genes <- unique(rv$markers$gene)
    pickerInput("selected_genes", "Select Genes to Plot:",
                choices = all_genes,
                selected = head(all_genes, 6),
                multiple = TRUE,
                options = list(
                  `live-search` = TRUE,
                  `actions-box` = TRUE,
                  `selected-text-format` = "count > 3"
                ))
  })

  output$marker_table <- renderDT({
    req(rv$markers)
    datatable(rv$markers,
              filter = "top",
              options = list(pageLength = 15, scrollX = TRUE),
              class = "compact stripe") %>%
      formatRound(columns = c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"),
                  digits  = 4)
  })

  output$marker_vln <- renderPlot({
    req(rv$seurat, input$selected_genes, length(input$selected_genes) > 0)
    VlnPlot(rv$seurat, features = input$selected_genes, pt.size = 0, ncol = 3) &
      theme_minimal(base_size = 12)
  })

  output$marker_feature <- renderPlot({
    req(rv$seurat, input$selected_genes, length(input$selected_genes) > 0)
    red <- if ("umap" %in% names(rv$seurat@reductions)) "umap"
    else if ("tsne" %in% names(rv$seurat@reductions)) "tsne"
    else "pca"
    FeaturePlot(rv$seurat,
                features  = input$selected_genes,
                reduction = red,
                ncol      = 3) &
      theme_minimal(base_size = 12)
  })

  output$marker_dot <- renderPlot({
    req(rv$seurat, input$selected_genes, length(input$selected_genes) > 0)
    DotPlot(rv$seurat, features = unique(input$selected_genes)) +
      coord_flip() +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  output$marker_heatmap <- renderPlot({
    req(rv$seurat, rv$markers)
    top_markers <- rv$markers %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = input$top_n_heatmap) %>%
      ungroup()
    DoHeatmap(rv$seurat, features = unique(top_markers$gene)) +
      theme(text = element_text(size = 11))
  })

  output$marker_ridge <- renderPlot({
    req(rv$seurat, input$selected_genes, length(input$selected_genes) > 0)
    RidgePlot(rv$seurat,
              features = input$selected_genes,
              ncol = 3) &
      theme_minimal(base_size = 12)
  })

  # --------------------------------------------------------------------------
  # TAB 7: CELL TYPE ANNOTATION
  # --------------------------------------------------------------------------

  output$celltype_inputs <- renderUI({
    req(rv$seurat, "seurat_clusters" %in% colnames(rv$seurat@meta.data))
    clusters <- levels(rv$seurat$seurat_clusters)
    inputs   <- lapply(clusters, function(cl) {
      textInput(
        inputId = paste0("celltype_", cl),
        label   = paste("Cluster", cl, ":"),
        value   = paste("Cluster", cl)
      )
    })
    do.call(tagList, inputs)
  })

  observeEvent(input$apply_celltypes, {
    req(rv$seurat, "seurat_clusters" %in% colnames(rv$seurat@meta.data))
    tryCatch({
      clusters <- levels(rv$seurat$seurat_clusters)
      new_ids  <- sapply(clusters, function(cl) {
        val <- input[[paste0("celltype_", cl)]]
        if (is.null(val) || trimws(val) == "") paste("Cluster", cl) else val
      })
      Idents(rv$seurat) <- "seurat_clusters"
      rv$seurat         <- RenameIdents(rv$seurat, new_ids)
      rv$seurat$cell_type <- as.character(Idents(rv$seurat))
      showNotification("Cell types applied!", type = "message", duration = 3)
    }, error = function(e) {
      showNotification(paste("Annotation Error:", conditionMessage(e)), type = "error", duration = 8)
    })
  })

  output$annotated_umap <- renderPlot({
    req(rv$seurat, "cell_type" %in% colnames(rv$seurat@meta.data))
    red <- if ("umap" %in% names(rv$seurat@reductions)) "umap"
    else if ("tsne" %in% names(rv$seurat@reductions)) "tsne"
    else "pca"
    DimPlot(rv$seurat,
            reduction = red,
            group.by  = "cell_type",
            label     = TRUE,
            repel     = TRUE,
            label.size = 4) +
      theme_minimal(base_size = 14) +
      labs(title = "Annotated Cell Types")
  })

  output$celltype_bar <- renderPlot({
    req(rv$seurat, "cell_type" %in% colnames(rv$seurat@meta.data))
    meta <- rv$seurat@meta.data
    meta$cell_type <- factor(meta$cell_type)
    counts_df <- meta %>%
      count(cell_type) %>%
      arrange(desc(n))
    ggplot(counts_df, aes(x = reorder(cell_type, n), y = n, fill = cell_type)) +
      geom_col(show.legend = FALSE, alpha = 0.9, width = 0.7) +
      geom_text(aes(label = n), hjust = -0.2, size = 4) +
      coord_flip() +
      scale_fill_brewer(palette = "Set2") +
      labs(title = "Cell Type Composition",
           x = "Cell Type", y = "Number of Cells") +
      theme_minimal(base_size = 14) +
      expand_limits(y = max(counts_df$n) * 1.15)
  })

  # Downloads
  output$download_rds <- downloadHandler(
    filename = function() paste0("seurat_annotated_", Sys.Date(), ".rds"),
    content  = function(file) {
      req(rv$seurat)
      saveRDS(rv$seurat, file)
    }
  )

  output$download_markers <- downloadHandler(
    filename = function() paste0("markers_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(rv$markers)
      write.csv(rv$markers, file, row.names = FALSE)
    }
  )

} # end server

# Run
shinyApp(ui = ui, server = server)
