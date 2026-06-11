# =============================================================================
# 08_roc_bootstrap.R
# - Receiver Operating Characteristic (ROC) curves for the two Elastic Net
#   hold-out tasks.
# - Cross-species external evaluation in the GSE227317 rat cohort using a
#   composite z-mean panel score (Tap1, Tap2, Cxcl10, Psmb8, Stat1).
# - 95% confidence intervals computed by 2000-replicate bootstrap.
# =============================================================================

source("00_setup.R")

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
  library(dplyr)
})

BOOT_N <- 2000

# ----------------------------------------------------------------------------
# Helper: ROC + bootstrap 95% CI
# ----------------------------------------------------------------------------
auc_with_ci <- function(labels, scores, n_boot = BOOT_N, seed = GLOBAL_SEED) {
  set.seed(seed)
  roc_obj <- roc(response = labels,
                 predictor = scores,
                 levels = c("Sham", "Stroke"),
                 direction = "<",
                 quiet = TRUE)
  ci_obj <- ci.auc(roc_obj, method = "bootstrap",
                   boot.n = n_boot, boot.stratified = TRUE)
  list(
    roc = roc_obj,
    auc = as.numeric(roc_obj$auc),
    ci_low  = as.numeric(ci_obj[1]),
    ci_high = as.numeric(ci_obj[3])
  )
}

# ----------------------------------------------------------------------------
# Elastic Net hold-out tasks
# ----------------------------------------------------------------------------
en_results <- readRDS(file.path(DIR_RESULTS, "ElasticNet_holdout_results.rds"))

eval_holdout <- function(res, label) {
  out <- auc_with_ci(res$y_test, res$pred_prob)
  cat(sprintf("[%s] AUC = %.3f (95%% CI: %.3f - %.3f)\n",
              label, out$auc, out$ci_low, out$ci_high))
  out
}

eval_spatial  <- eval_holdout(en_results$spatial,  "Cortex -> Striatum")
eval_genotype <- eval_holdout(en_results$genotype, "WT -> Gpx1-/- (KO)")

# ----------------------------------------------------------------------------
# Cross-species external evaluation in GSE227317 (rat -> mouse symbol mapped)
# Composite panel score = z-mean of Tap1, Tap2, Cxcl10, Psmb8, Stat1
# ----------------------------------------------------------------------------
PANEL_GENES <- c("Tap1", "Tap2", "Cxcl10", "Psmb8", "Stat1")

expr_rat  <- readRDS(file.path(DIR_DATA, "GSE227317", "GSE227317_expr.rds"))
pheno_rat <- readRDS(file.path(DIR_DATA, "GSE227317", "GSE227317_pheno.rds"))

# Group label: 1 = Stroke, 0 = Sham; edit regex if your phenoData differs.
group_field <- grep("characteristics_ch1|title|source_name",
                    colnames(pheno_rat), value = TRUE, ignore.case = TRUE)
group_label <- as.character(pheno_rat[, group_field[1]])
group_rat <- factor(
  ifelse(grepl("(?i)ischemia|stroke|MCAO", group_label, perl = TRUE),
         "Stroke", "Sham"),
  levels = c("Sham", "Stroke")
)
names(group_rat) <- rownames(pheno_rat)
group_rat <- group_rat[colnames(expr_rat)]

panel_present <- intersect(PANEL_GENES, rownames(expr_rat))
cat(sprintf("Panel genes available in GSE227317: %d / %d\n",
            length(panel_present), length(PANEL_GENES)))

z_mat       <- t(scale(t(expr_rat[panel_present, , drop = FALSE])))
panel_score <- colMeans(z_mat, na.rm = TRUE)

eval_panel <- auc_with_ci(group_rat, panel_score)
cat(sprintf("[Rat panel score] AUC = %.3f (95%% CI: %.3f - %.3f)\n",
            eval_panel$auc, eval_panel$ci_low, eval_panel$ci_high))

# Per-gene AUC in the rat cohort
per_gene_auc <- lapply(panel_present, function(g) {
  out <- auc_with_ci(group_rat, expr_rat[g, ])
  data.frame(gene = g,
             AUC  = out$auc,
             CI_low  = out$ci_low,
             CI_high = out$ci_high)
}) |> do.call(rbind, args = _)

write.csv(per_gene_auc,
          file.path(DIR_RESULTS, "Rat_panel_per_gene_AUC.csv"),
          row.names = FALSE)

# ----------------------------------------------------------------------------
# Save AUC summary
# ----------------------------------------------------------------------------
auc_summary <- data.frame(
  task = c("Cortex_to_Striatum", "WT_to_KO", "Rat_panel"),
  AUC  = c(eval_spatial$auc, eval_genotype$auc, eval_panel$auc),
  CI_low  = c(eval_spatial$ci_low,  eval_genotype$ci_low,  eval_panel$ci_low),
  CI_high = c(eval_spatial$ci_high, eval_genotype$ci_high, eval_panel$ci_high)
)
write.csv(auc_summary,
          file.path(DIR_RESULTS, "AUC_summary_with_bootstrap_CI.csv"),
          row.names = FALSE)

# ----------------------------------------------------------------------------
# ROC plots
# ----------------------------------------------------------------------------
plot_roc <- function(roc_obj, title, file_path) {
  df <- data.frame(
    spec = rev(roc_obj$specificities),
    sens = rev(roc_obj$sensitivities)
  )
  p <- ggplot(df, aes(1 - spec, sens)) +
    geom_line(size = 1, color = "#E64B35") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey50") +
    labs(title = title, x = "1 - specificity", y = "sensitivity") +
    theme_minimal(base_size = 12) +
    coord_equal()
  ggsave(file_path, p, width = 5, height = 5)
}
plot_roc(eval_spatial$roc,
         sprintf("Cortex -> Striatum (AUC = %.3f)", eval_spatial$auc),
         file.path(DIR_FIG, "FigS2A_ROC_Cortex_to_Striatum.pdf"))
plot_roc(eval_genotype$roc,
         sprintf("WT -> KO (AUC = %.3f)", eval_genotype$auc),
         file.path(DIR_FIG, "FigS2B_ROC_WT_to_KO.pdf"))
plot_roc(eval_panel$roc,
         sprintf("Rat panel score (AUC = %.3f)", eval_panel$auc),
         file.path(DIR_FIG, "FigS3_ROC_rat_panel.pdf"))

cat("\nROC + bootstrap evaluation complete.\n")
cat("Outputs in results/ and figures/.\n")
