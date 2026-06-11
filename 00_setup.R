# =============================================================================
# 00_setup.R
# Install required packages and set global configuration.
# Run this script ONCE before running the rest of the pipeline.
# =============================================================================

# ---- CRAN packages ----------------------------------------------------------
cran_pkgs <- c(
  "glmnet",      # Elastic Net logistic regression
  "pROC",        # ROC analysis and bootstrap CIs
  "WGCNA",       # weighted gene co-expression network analysis
  "umap",        # UMAP visualization
  "ggplot2",     # plotting
  "dplyr",
  "tibble",
  "readr",
  "Matrix"
)

cran_missing <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(cran_missing) > 0) {
  install.packages(cran_missing, repos = "https://cloud.r-project.org")
}

# ---- Bioconductor packages --------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c(
  "GEOquery",        # download datasets from GEO
  "limma",           # differential expression
  "sva",             # ComBat batch correction
  "GSVA",            # ssGSEA
  "clusterProfiler", # GO/KEGG enrichment
  "org.Mm.eg.db"     # mouse gene annotations
)

bioc_missing <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(bioc_missing) > 0) {
  BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
}

# ---- Global configuration ---------------------------------------------------
GLOBAL_SEED <- 12345
set.seed(GLOBAL_SEED)

# Project directory layout (relative to repository root)
DIR_DATA    <- "data"
DIR_RESULTS <- "results"
DIR_FIG     <- "figures"

dir.create(DIR_DATA,    showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG,     showWarnings = FALSE, recursive = TRUE)

# Save session information snapshot
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

cat("Setup complete. Seed =", GLOBAL_SEED, "\n")
