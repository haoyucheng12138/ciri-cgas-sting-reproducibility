# =============================================================================
# 00_config.R  -  paths, seed, package loading. Sourced by every other script.
# Run from the repository root:  Rscript R/run_all.R
# =============================================================================
options(stringsAsFactors = FALSE, timeout = 3600)
SEED <- 12345
set.seed(SEED)

DIR_RAW <- "data/raw"          # GEO files downloaded by 01_download_geo.R
DIR_PRO <- "data/processed"    # processed matrices
DIR_RES <- "results"           # all tables reported in the manuscript
DIR_REF <- "reference"         # bundled reference files (gene lists, annotations)
for (d in c(DIR_RAW, DIR_PRO, DIR_RES)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

res_dir <- function(step) { p <- file.path(DIR_RES, step); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
save_session <- function(step) writeLines(capture.output(sessionInfo()), file.path(res_dir(step), "sessionInfo.txt"))
read_mat <- function(f) as.matrix(read.csv(f, row.names = 1, check.names = FALSE))

## Median of rows sharing the same identifier (single-row identifiers copied directly).
collapse_median <- function(m, ids) {
  ok <- !is.na(ids) & ids != ""
  m <- m[ok, , drop = FALSE]; ids <- ids[ok]
  tab <- table(ids); single <- names(tab)[tab == 1]; multi <- names(tab)[tab > 1]
  out <- m[match(single, ids), , drop = FALSE]; rownames(out) <- single
  if (length(multi)) {
    out <- rbind(out, t(vapply(multi, function(s) apply(m[ids == s, , drop = FALSE], 2, median), numeric(ncol(m)))))
  }
  out[order(rownames(out)), , drop = FALSE]
}

## Row-wise z-score then column mean (composite score).
zmean <- function(m) colMeans(t(scale(t(m))))
