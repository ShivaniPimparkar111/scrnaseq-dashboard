# ==============================================================================
# install_packages.R
# Installs all R packages required for the scRNA-seq Analysis Dashboard
# Author: Shivani Pimparkar
# ==============================================================================

cat("=== scRNA-seq Dashboard — Package Installation ===\n\n")

# CRAN packages required by the dashboard
pkgs <- c(
  "shiny",
  "shinydashboard",
  "Seurat",
  "ggplot2",
  "dplyr",
  "DT",
  "shinycssloaders",
  "shinyWidgets",
  "patchwork",
  "ggrepel",
  "scales",
  "BiocManager",
  "markdown"     # needed for includeMarkdown()
)

cat("Installing CRAN packages...\n")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  Installing: %s\n", p))
    install.packages(p, dependencies = TRUE, quiet = TRUE)
  } else {
    cat(sprintf("  Already installed: %s\n", p))
  }
}

# Bioconductor packages (recommended for Seurat v5)
cat("\nInstalling Bioconductor dependencies...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c("limma", "glmGamPoi", "BiocGenerics", "S4Vectors",
               "SummarizedExperiment", "SingleCellExperiment")

BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)

# Optional but recommended packages for enhanced Seurat functionality
cat("\nInstalling optional enhancement packages...\n")
optional_pkgs <- c("future", "future.apply", "MASS", "KernSmooth",
                   "fitdistrplus", "leiden", "uwot", "Rtsne", "igraph",
                   "lmtest", "RcppAnnoy", "RcppHNSW")

for (p in optional_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    tryCatch({
      install.packages(p, dependencies = TRUE, quiet = TRUE)
      cat(sprintf("  Installed: %s\n", p))
    }, error = function(e) {
      cat(sprintf("  Could not install %s (optional): %s\n", p, conditionMessage(e)))
    })
  } else {
    cat(sprintf("  Already installed: %s\n", p))
  }
}

cat("\n=== Installation complete! ===\n")
cat("Run the dashboard with: shiny::runApp('.')\n")
