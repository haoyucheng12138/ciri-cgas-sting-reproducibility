# =============================================================================
# 06_ssgsea_sting_score.R
# - ssGSEA on the GSE131193 discovery cohort using a 27-cell-type immune
#   reference panel.
# - Composite STING activation score = mean z-score of five core axis genes:
#       Mb21d1 (cGAS), Tmem173 (STING), Irf3, Isg15, Cxcl10
#   (Trex1 was excluded because it was not detected in GSE131193.)
# - Mann-Whitney U tests with BH-FDR correction for differential immune
#   enrichment between Stroke and Sham.
# - Spearman correlation between the composite STING score and individual
#   immune-cell enrichment scores (BH-FDR corrected).
# =============================================================================

source("00_setup.R")

suppressPackageStartupMessages({
  library(GSVA)
  library(dplyr)
  library(tibble)
})

expr  <- readRDS(file.path(DIR_DATA, "GSE131193", "GSE131193_expr.rds"))
pheno <- readRDS(file.path(DIR_DATA, "GSE131193", "GSE131193_pheno.rds"))

group_field <- grep("characteristics_ch1|title|source_name",
                    colnames(pheno), value = TRUE, ignore.case = TRUE)
group_label <- as.character(pheno[, group_field[1]])
group <- factor(
  ifelse(grepl("(?i)ischemia|stroke|MCAO", group_label, perl = TRUE),
         "Stroke", "Sham"),
  levels = c("Sham", "Stroke")
)
names(group) <- rownames(pheno)
group <- group[colnames(expr)]

# ----------------------------------------------------------------------------
# 27-cell-type immune signature panel
# Replace the placeholder list below with the full 27-gene-set panel used in
# the manuscript (e.g. Bindea et al. 2013 / Charoentong et al. 2017). The
# structure is: list(cell_type_name = c("gene1", "gene2", ...))
# A complete reference panel can be downloaded from the original publication
# supplementary materials and saved as data/immune_signatures_27.rds.
# ----------------------------------------------------------------------------
sig_path <- file.path(DIR_DATA, "immune_signatures_27.rds")
if (!file.exists(sig_path)) {
  stop("Immune signature panel not found at ", sig_path,
       ". Please supply the 27-cell-type signature list as an RDS file ",
       "(named list of character vectors, gene symbols on mouse genome).")
}
immune_signatures <- readRDS(sig_path)

# ----------------------------------------------------------------------------
# ssGSEA
# ----------------------------------------------------------------------------
ssgsea_scores <- gsva(
  expr      = as.matrix(expr),
  gset.idx.list = immune_signatures,
  method    = "ssgsea",
  kcdf      = "Gaussian",
  ssgsea.norm = TRUE,
  verbose   = FALSE
)

write.csv(ssgsea_scores,
          file.path(DIR_RESULTS, "ssGSEA_immune_scores.csv"),
          row.names = TRUE)

# Mann-Whitney U tests (Stroke vs Sham) with BH-FDR
mw_test <- function(vec, grp) {
  wilcox.test(vec ~ grp, exact = FALSE)$p.value
}
pvals <- apply(ssgsea_scores, 1, mw_test, grp = group)
res_immune <- data.frame(
  cell_type = rownames(ssgsea_scores),
  p_value   = pvals,
  fdr       = p.adjust(pvals, method = "BH")
)
write.csv(res_immune,
          file.path(DIR_RESULTS, "ssGSEA_diff_test.csv"),
          row.names = FALSE)

# ----------------------------------------------------------------------------
# Composite STING activation score (z-mean of five core axis genes)
# ----------------------------------------------------------------------------
sting_genes <- c("Mb21d1", "Tmem173", "Irf3", "Isg15", "Cxcl10")
sting_present <- intersect(sting_genes, rownames(expr))

z_mat <- t(scale(t(expr[sting_present, , drop = FALSE])))
sting_score <- colMeans(z_mat, na.rm = TRUE)

# Group comparison
sting_test <- wilcox.test(sting_score ~ group, exact = FALSE)
cat(sprintf("STING composite score Sham vs Stroke: p = %.3e (Mann-Whitney)\n",
            sting_test$p.value))

write.csv(
  data.frame(sample_id = names(sting_score),
             group = as.character(group),
             sting_score = sting_score),
  file.path(DIR_RESULTS, "STING_composite_score.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------------
# Spearman correlation between STING score and each immune-cell signature
# ----------------------------------------------------------------------------
cor_results <- t(apply(ssgsea_scores, 1, function(x) {
  res <- suppressWarnings(cor.test(x, sting_score, method = "spearman"))
  c(rho = unname(res$estimate), p_value = res$p.value)
}))
cor_df <- data.frame(
  cell_type = rownames(cor_results),
  rho       = cor_results[, "rho"],
  p_value   = cor_results[, "p_value"]
)
cor_df$fdr <- p.adjust(cor_df$p_value, method = "BH")
write.csv(cor_df,
          file.path(DIR_RESULTS, "STING_immune_correlation.csv"),
          row.names = FALSE)

cat("\nssGSEA + STING composite score analysis complete.\n")
