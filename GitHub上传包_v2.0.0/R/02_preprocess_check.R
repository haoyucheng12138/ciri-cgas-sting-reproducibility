# =============================================================================
# 02_preprocess_check.R  (optional; used by the authors, not needed by readers)
# Compares the matrices built by 02_preprocess.R with the matrices of the original LKFlow analysis.
# Set LKFLOW_DIR to the folder holding the original *_final.csv files.
# =============================================================================
source("R/00_config.R")
LKFLOW_DIR <- Sys.getenv("LKFLOW_DIR", "")
if (!nzchar(LKFLOW_DIR)) { cat("LKFLOW_DIR not set; skipping comparison with the original platform output.\n"); quit(save = "no") }
## LKFlow wrote the Excel-damaged symbols as "1-Mar" (GSE131712/GSE202659) or "MARCH1" (GSE131193); map both to "March1".
lkread <- function(f) { x <- read_mat(file.path(LKFLOW_DIR, f))
  rownames(x) <- sub("^MARCH(\\d+)$", "March\\1", sub("^(\\d+)-Mar$", "March\\1", rownames(x))); x }
same <- function(mine, lk, tag) {
  g <- intersect(rownames(mine), rownames(lk))
  cat(sprintf("%-10s genes mine %5d / LKFlow %5d / identical set %s | sample labels identical %s | max|diff| %.2e\n", tag, nrow(mine), nrow(lk),
              setequal(rownames(mine), rownames(lk)), identical(colnames(mine), colnames(lk)), max(abs(mine[g, colnames(lk)] - lk[g, ]))))
}
m <- read_mat(file.path(DIR_PRO, "merged_log2FPKM.csv"))
same(m[, 1:6], lkread("GSE131712_final.csv"), "GSE131712")
same(m[, 7:12], lkread("GSE202659_final.csv"), "GSE202659")
same(read_mat(file.path(DIR_PRO, "GSE23163_log2ratio.csv")), lkread("GSE23163_final.csv"), "GSE23163")
MERGED_LK <- Sys.getenv("LKFLOW_MERGED", "")
if (nzchar(MERGED_LK)) { x <- read_mat(MERGED_LK); rownames(x) <- sub("^(\\d+)-Mar$", "March\\1", rownames(x)); same(m, x[, colnames(m)], "merged") }

mine <- read_mat(file.path(DIR_PRO, "GSE131193_log2counts.csv")); lk <- lkread("GSE131193_final.csv")
g <- intersect(rownames(mine), rownames(lk))
cat(sprintf("GSE131193  genes mine %5d / LKFlow %5d / identical set %s\n", nrow(mine), nrow(lk), setequal(rownames(mine), rownames(lk))))
hold <- sapply(colnames(mine), function(s) { d <- apply(abs(lk[g, ] - mine[g, s]), 2, max); c(names(which.min(d)), min(d)) })
tab <- data.frame(GSM = colnames(mine), LKFlow_label_on_same_data = hold[1, ], max_abs_diff = signif(as.numeric(hold[2, ]), 3))
tab$label_correct <- tab$GSM == tab$LKFlow_label_on_same_data
print(tab, row.names = FALSE)
cat(sprintf("GSE131193: every column has an exact copy in the LKFlow matrix: %s | LKFlow label correct for %d of 24 samples\n",
            all(as.numeric(hold[2, ]) < 1e-8), sum(tab$label_correct)))
write.csv(tab, file.path(res_dir("02_preprocess"), "GSE131193_LKFlow_label_audit.csv"), row.names = FALSE)
