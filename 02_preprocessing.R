# =============================================================================
# 02_preprocessing.R
# - log2(TPM + 1) transformation for RNA-seq datasets (or log2 for microarray)
# - Row-wise median imputation for missing values
# - Cross-species symbol mapping (rat -> mouse) using NCBI HomoloGene Build 68
# Outputs are saved as RDS expression matrices keyed by GEO accession.
# =============================================================================

source("00_setup.R")

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(tibble)
})

# ----------------------------------------------------------------------------
# Helper: extract expression matrix and phenoData from GEOquery output
# ----------------------------------------------------------------------------
extract_expr <- function(gset_list) {
  # GEOquery returns a list with one ExpressionSet per platform
  if (is.null(gset_list) || length(gset_list) == 0) {
    stop("Empty GEO object.")
  }
  eset <- gset_list[[1]]
  expr <- as.matrix(Biobase::exprs(eset))
  pheno <- Biobase::pData(eset)
  list(expr = expr, pheno = pheno)
}

# ----------------------------------------------------------------------------
# Helper: log2(x + 1) transform (for raw TPM / count-derived matrices)
# ----------------------------------------------------------------------------
log2p1 <- function(mat) {
  log2(pmax(mat, 0) + 1)
}

# ----------------------------------------------------------------------------
# Helper: row-wise median imputation for missing values
# ----------------------------------------------------------------------------
impute_row_median <- function(mat) {
  if (!any(is.na(mat))) return(mat)
  for (i in seq_len(nrow(mat))) {
    row_vals <- mat[i, ]
    if (any(is.na(row_vals))) {
      med <- median(row_vals, na.rm = TRUE)
      row_vals[is.na(row_vals)] <- med
      mat[i, ] <- row_vals
    }
  }
  mat
}

# ----------------------------------------------------------------------------
# Helper: ensure expression values are on log2 scale
# Heuristic: if max(value) > 50, assume linear scale and apply log2(x+1).
# ----------------------------------------------------------------------------
ensure_log_scale <- function(mat) {
  if (max(mat, na.rm = TRUE) > 50) {
    mat <- log2p1(mat)
  }
  mat
}

# ----------------------------------------------------------------------------
# Helper: rat-to-mouse ortholog mapping using NCBI HomoloGene Build 68
# Build 68 is downloaded as a flat file from NCBI FTP.
# ----------------------------------------------------------------------------
build_homologene_map <- function() {
  cache_file <- file.path(DIR_DATA, "homologene_b68.rds")
  if (file.exists(cache_file)) return(readRDS(cache_file))

  url <- "https://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.data"
  tmp <- tempfile()
  utils::download.file(url, destfile = tmp, mode = "wb")

  hg <- read.table(
    tmp, sep = "\t", header = FALSE, quote = "",
    col.names = c("hid", "taxid", "gene_id", "symbol", "protein_gi", "protein_acc"),
    stringsAsFactors = FALSE
  )

  # taxid 10090 = Mus musculus, 10116 = Rattus norvegicus
  rat   <- hg[hg$taxid == 10116, c("hid", "symbol")]
  mouse <- hg[hg$taxid == 10090, c("hid", "symbol")]
  colnames(rat)   <- c("hid", "rat_symbol")
  colnames(mouse) <- c("hid", "mouse_symbol")

  map <- merge(rat, mouse, by = "hid")
  map <- map[!duplicated(map$rat_symbol), ]

  saveRDS(map, cache_file)
  map
}

map_rat_to_mouse <- function(expr_rat) {
  map <- build_homologene_map()
  rownames_rat <- rownames(expr_rat)
  matched <- map[match(rownames_rat, map$rat_symbol), ]
  keep <- !is.na(matched$mouse_symbol)

  expr_mapped <- expr_rat[keep, , drop = FALSE]
  rownames(expr_mapped) <- matched$mouse_symbol[keep]

  # collapse duplicate mouse symbols by mean
  if (any(duplicated(rownames(expr_mapped)))) {
    expr_mapped <- aggregate(expr_mapped,
                             by = list(symbol = rownames(expr_mapped)),
                             FUN = mean)
    rownames(expr_mapped) <- expr_mapped$symbol
    expr_mapped$symbol <- NULL
    expr_mapped <- as.matrix(expr_mapped)
  }
  expr_mapped
}

# ----------------------------------------------------------------------------
# Main loop: process each downloaded GEO dataset
# ----------------------------------------------------------------------------
GEO_ACCESSIONS <- c("GSE131193", "GSE131712", "GSE202659", "GSE23163", "GSE227317")

for (acc in GEO_ACCESSIONS) {
  rds_file <- file.path(DIR_DATA, acc, paste0(acc, "_series.rds"))
  if (!file.exists(rds_file)) {
    cat(sprintf("[%s] series RDS not found, skipping. Run 01_download_geo.R first.\n", acc))
    next
  }
  gset <- readRDS(rds_file)
  parsed <- extract_expr(gset)
  expr  <- parsed$expr
  pheno <- parsed$pheno

  expr <- impute_row_median(expr)
  expr <- ensure_log_scale(expr)

  # Cross-species mapping for the rat dataset
  if (acc == "GSE227317") {
    expr <- map_rat_to_mouse(expr)
  }

  out_expr  <- file.path(DIR_DATA, acc, paste0(acc, "_expr.rds"))
  out_pheno <- file.path(DIR_DATA, acc, paste0(acc, "_pheno.rds"))
  saveRDS(expr,  out_expr)
  saveRDS(pheno, out_pheno)
  cat(sprintf("[%s] preprocessed: %d genes x %d samples\n",
              acc, nrow(expr), ncol(expr)))
}

cat("\nPreprocessing complete.\n")
