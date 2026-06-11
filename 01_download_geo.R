# =============================================================================
# 01_download_geo.R
# Download all five GEO datasets used in the manuscript.
# Datasets are stored under data/<accession>/.
# =============================================================================

source("00_setup.R")

suppressPackageStartupMessages({
  library(GEOquery)
})

GEO_ACCESSIONS <- c(
  "GSE131193",   # Discovery cohort (mouse)
  "GSE131712",   # Merged validation cohort 1 (mouse)
  "GSE202659",   # Merged validation cohort 2 (mouse)
  "GSE23163",    # Machine learning cohort (mouse)
  "GSE227317"    # Cross-species external evaluation (rat)
)

download_geo <- function(accession, destdir = DIR_DATA) {
  out_dir <- file.path(destdir, accession)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  cat(sprintf("[%s] downloading...\n", accession))

  # Series matrix file (contains expression matrix + phenoData)
  gset <- tryCatch(
    getGEO(accession, GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = out_dir),
    error = function(e) {
      cat(sprintf("[%s] ERROR: %s\n", accession, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(gset)) return(invisible(NULL))

  # Save series matrix as RDS for downstream scripts
  saveRDS(gset, file = file.path(out_dir, paste0(accession, "_series.rds")))

  # For RNA-seq datasets, also try to download supplementary count files
  if (accession %in% c("GSE131193", "GSE131712", "GSE202659", "GSE227317")) {
    tryCatch(
      getGEOSuppFiles(accession, baseDir = destdir, makeDirectory = FALSE),
      error = function(e) {
        cat(sprintf("[%s] supplementary file download skipped: %s\n",
                    accession, conditionMessage(e)))
      }
    )
  }

  cat(sprintf("[%s] done. Files in %s\n", accession, out_dir))
}

for (acc in GEO_ACCESSIONS) {
  download_geo(acc)
}

cat("\nAll downloads completed.\n")
