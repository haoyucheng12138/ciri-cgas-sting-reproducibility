# =============================================================================
# 07_elasticnet_modeling.R
# Elastic Net logistic regression on the GSE23163 cohort.
#   - alpha = 0.4 (selected via comparative evaluation across alpha = 0.1 to 1.0
#     and alternative model types; see manuscript Methods)
#   - lambda chosen by 10-fold cross-validation, lambda.min criterion
#   - random seed = 12345
# Two biologically motivated stratified hold-out tasks are evaluated:
#   (1) Spatial generalization (Cortex -> Striatum)
#   (2) Genotypic generalization (WT -> Gpx1-/-)
# =============================================================================

source("00_setup.R")

suppressPackageStartupMessages({
  library(glmnet)
  library(dplyr)
  library(tibble)
})

# ----------------------------------------------------------------------------
# Curated candidate feature space (39 genes)
# Composition (per manuscript Methods):
#   - cGAS-STING / type I IFN axis genes
#   - MHC-I antigen processing and presentation components
#   - PPI-network hub genes
#   - Trpc6 as an a priori candidate of interest
# Genes not detectable in the GSE23163 expression matrix are dropped at runtime.
# ----------------------------------------------------------------------------
CURATED_FEATURES <- c(
  # MHC-I antigen processing and presentation
  "Tapbp", "Tap1", "Tap2", "Tapbpl",
  "H2-D1", "H2-K1", "H2-T23",
  # Translation / ribosomal hub genes from PPI
  "Eif4a1", "Eif6", "Rplp0", "Rpl18",
  # IFN response / ISG core
  "Ifi35", "Trim56", "Stat1", "Stat2", "Irf7", "Irf9", "Isg15",
  "Usp18", "Rsad2", "Mx1", "Oas1a", "Oas2", "Ifit3", "Ifih1", "Ddx58", "Nlrc5",
  # cGAS-STING associated chemokines
  "Cxcl9", "Cxcl10",
  # Immunoproteasome
  "Psmb8", "Psmb10",
  # GTPase / interferon-induced GTPases
  "Igtp", "Irgm1", "Irgm2",
  # Vasculature / context genes
  "Cdh5",
  # Chemokine
  "Ccl5",
  # Macrophage / immune checkpoint
  "Cd163", "Cd274",
  # A priori candidate
  "Trpc6"
)

# ----------------------------------------------------------------------------
# Load GSE23163 expression and phenoData
# ----------------------------------------------------------------------------
expr  <- readRDS(file.path(DIR_DATA, "GSE23163", "GSE23163_expr.rds"))
pheno <- readRDS(file.path(DIR_DATA, "GSE23163", "GSE23163_pheno.rds"))

# Build sample annotation (region: Cortex/Striatum; genotype: WT/KO; group: Sham/Stroke)
# Edit the regular expressions below if your phenoData uses different labels.
parse_field <- function(pheno, regex, label_yes, label_no) {
  txt <- apply(pheno, 1, function(r) paste(as.character(r), collapse = " | "))
  ifelse(grepl(regex, txt, ignore.case = TRUE, perl = TRUE),
         label_yes, label_no)
}
sample_meta <- tibble(
  sample_id = rownames(pheno),
  region   = parse_field(pheno, "(?i)cortex",       "Cortex",  "Striatum"),
  genotype = parse_field(pheno, "(?i)Gpx1|knockout|\\bKO\\b", "KO", "WT"),
  group    = parse_field(pheno, "(?i)ischemia|stroke|MCAO",   "Stroke", "Sham")
)

# Restrict to the 32 samples used in the manuscript (cortex + striatum,
# WT + Gpx1-/-). Adjust the filter if your phenoData annotations differ.
sample_meta <- sample_meta |>
  filter(region %in% c("Cortex", "Striatum"),
         genotype %in% c("WT", "KO"),
         group %in% c("Sham", "Stroke"))

stopifnot(nrow(sample_meta) > 0)

# Restrict expression matrix to curated features that are detectable
features_present <- intersect(CURATED_FEATURES, rownames(expr))
cat(sprintf("Curated features available in GSE23163: %d / %d\n",
            length(features_present), length(CURATED_FEATURES)))

X_full <- t(expr[features_present, sample_meta$sample_id, drop = FALSE])
y_full <- factor(sample_meta$group, levels = c("Sham", "Stroke"))

# ----------------------------------------------------------------------------
# Elastic Net training/evaluation function
# ----------------------------------------------------------------------------
fit_enet_holdout <- function(X_train, y_train, X_test, y_test,
                             alpha = 0.4, nfolds = 10, seed = GLOBAL_SEED) {
  set.seed(seed)
  cv_fit <- cv.glmnet(
    x       = as.matrix(X_train),
    y       = y_train,
    family  = "binomial",
    alpha   = alpha,
    nfolds  = nfolds,
    standardize = TRUE,
    type.measure = "deviance"
  )
  best_lambda <- cv_fit$lambda.min

  # Coefficients at lambda.min
  coef_vec <- as.matrix(coef(cv_fit, s = "lambda.min"))
  coef_df  <- data.frame(
    feature     = rownames(coef_vec),
    coefficient = as.numeric(coef_vec[, 1])
  )
  coef_df <- coef_df[coef_df$feature != "(Intercept)" &
                     coef_df$coefficient != 0, ]
  coef_df <- coef_df[order(-abs(coef_df$coefficient)), ]

  # Test predictions
  pred_prob <- predict(cv_fit, newx = as.matrix(X_test),
                       s = "lambda.min", type = "response")[, 1]
  list(
    cv_fit      = cv_fit,
    best_lambda = best_lambda,
    coefficients = coef_df,
    pred_prob   = pred_prob,
    y_test      = y_test
  )
}

# ----------------------------------------------------------------------------
# Task 1: spatial generalization (Cortex -> Striatum)
# ----------------------------------------------------------------------------
idx_cortex   <- which(sample_meta$region == "Cortex")
idx_striatum <- which(sample_meta$region == "Striatum")
res_spatial <- fit_enet_holdout(
  X_train = X_full[idx_cortex, , drop = FALSE],
  y_train = y_full[idx_cortex],
  X_test  = X_full[idx_striatum, , drop = FALSE],
  y_test  = y_full[idx_striatum]
)

# ----------------------------------------------------------------------------
# Task 2: genotypic generalization (WT -> KO)
# ----------------------------------------------------------------------------
idx_wt <- which(sample_meta$genotype == "WT")
idx_ko <- which(sample_meta$genotype == "KO")
res_genotype <- fit_enet_holdout(
  X_train = X_full[idx_wt, , drop = FALSE],
  y_train = y_full[idx_wt],
  X_test  = X_full[idx_ko, , drop = FALSE],
  y_test  = y_full[idx_ko]
)

# ----------------------------------------------------------------------------
# Save outputs (coefficients + predictions for downstream ROC/bootstrap script)
# ----------------------------------------------------------------------------
write.csv(res_spatial$coefficients,
          file.path(DIR_RESULTS, "ElasticNet_Cortex_to_Striatum_coefficients.csv"),
          row.names = FALSE)
write.csv(res_genotype$coefficients,
          file.path(DIR_RESULTS, "ElasticNet_WT_to_KO_coefficients.csv"),
          row.names = FALSE)

saveRDS(list(spatial = res_spatial, genotype = res_genotype),
        file.path(DIR_RESULTS, "ElasticNet_holdout_results.rds"))

cat(sprintf("\nSpatial task: lambda.min = %.4f, %d non-zero features\n",
            res_spatial$best_lambda, nrow(res_spatial$coefficients)))
cat(sprintf("Genotypic task: lambda.min = %.4f, %d non-zero features\n",
            res_genotype$best_lambda, nrow(res_genotype$coefficients)))
cat("\nElastic Net modeling complete. Run 08_roc_bootstrap.R for AUC + 95% CIs.\n")
