# =============================================================================
# 00_config.R  -  paths, seed, package loading. Sourced by every other script.
# Run from the repository root:  Rscript R/run_all.R
# =============================================================================
options(stringsAsFactors = FALSE, timeout = 600)
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

## Display symbols for figures and supplementary tables: current official MGI symbol (org.Mm.eg.db). Official symbols are kept;
## a former symbol that is an alias of exactly one official symbol is replaced (e.g. Tmem173 -> Sting1, Mb21d1 -> Cgas);
## anything else (non-gene labels, ambiguous aliases) is returned unchanged. Analysis files keep the symbols of the GEO files.
official_symbol <- function(x) {
  x <- as.character(x); db <- org.Mm.eg.db::org.Mm.eg.db
  todo <- setdiff(unique(x[!is.na(x)]), AnnotationDbi::keys(db, "SYMBOL"))
  todo <- intersect(todo, AnnotationDbi::keys(db, "ALIAS"))
  if (!length(todo)) return(x)
  al <- suppressMessages(AnnotationDbi::select(db, todo, "SYMBOL", "ALIAS"))
  al <- unique(al[!is.na(al$SYMBOL), c("ALIAS", "SYMBOL")]); n <- table(al$ALIAS); al <- al[al$ALIAS %in% names(n)[n == 1], ]
  i <- match(x, al$ALIAS); x[!is.na(i)] <- al$SYMBOL[i[!is.na(i)]]
  x
}
