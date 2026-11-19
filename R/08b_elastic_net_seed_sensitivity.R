# =============================================================================
# 08b_elastic_net_seed_sensitivity.R  -  dependence of the Elastic Net results on the random CV fold assignment.
#  Same models as 08_elastic_net.R (alpha 0.4, 10-fold CV, lambda.min); only the seed of the fold assignment changes:
#  seeds 1-200 for the 39-gene and rebuilt panels, 1-30 for the full transcriptome (computational cost).
#  With 16 training samples, 10-fold CV leaves 1-2 samples per fold; glmnet warns that one class has < 8 observations
#  and enforces grouped = FALSE. These warnings are expected for this design and are recorded, not suppressed silently.
# =============================================================================
source("R/00_config.R")
suppressMessages({ library(glmnet); library(pROC) })
out <- res_dir("08_elastic_net")
x <- read_mat(file.path(DIR_PRO, "GSE23163_log2ratio.csv"))
meta <- read.csv(file.path(DIR_PRO, "GSE23163_sample_annotation.csv"))
meta <- meta[meta$group %in% c("Sham", "IR24h"), ]; meta$y <- as.integer(meta$group == "IR24h")
tasks <- list(Cortex_to_Striatum = list(train = meta$ID[meta$region == "Cortex"], test = meta$ID[meta$region == "Striatum"]),
              WT_to_Gpx1KO = list(train = meta$ID[meta$genotype == "WT"], test = meta$ID[meta$genotype == "Gpx1KO"]))
y_of <- function(ids) meta$y[match(ids, meta$ID)]
panel39 <- setdiff(read.csv(file.path(DIR_REF, "original_39gene_candidate_panel.csv"))$Gene, "Gene")
rebuilt <- readLines(file.path(out, "rebuilt_panel_definition.txt")); rebuilt <- unique(rebuilt[!grepl("^#", rebuilt)])
spaces <- list(A_original_39 = panel39, B_full_transcriptome = rownames(x), C_rebuilt_panel = rebuilt)
en <- read.csv(file.path(out, "elastic_net_summary.csv"))
warn_log <- list(); per_seed <- list(); summ <- list()
for (s in names(spaces)) for (t in names(tasks)) {
  tr <- tasks[[t]]$train; te <- tasks[[t]]$test
  g <- intersect(spaces[[s]], rownames(x)); g <- g[apply(x[g, tr, drop = FALSE], 1, sd) > 0]
  X <- t(x[g, tr]); y <- y_of(tr); Xt <- t(x[g, te]); yt <- y_of(te)
  seeds <- if (s == "B_full_transcriptome") 1:30 else 1:200
  r <- t(sapply(seeds, function(sd) {
    set.seed(sd)
    cv <- withCallingHandlers(cv.glmnet(X, y, family = "binomial", alpha = 0.4, nfolds = 10, type.measure = "deviance"),
                              warning = function(w) { warn_log[[length(warn_log) + 1]] <<- conditionMessage(w); invokeRestart("muffleWarning") })
    p <- as.numeric(predict(cv, newx = Xt, s = "lambda.min", type = "response"))
    c(seed = sd, auc = as.numeric(auc(roc(yt, p, levels = c(0, 1), direction = "<", quiet = TRUE))),
      nonzero = sum(as.matrix(coef(cv, s = "lambda.min"))[-1, 1] != 0), jerr = cv$glmnet.fit$jerr)
  }))
  per_seed[[paste(s, t)]] <- data.frame(feature_space = s, task = t, r)
  ref <- en[en$feature_space == s & en$task == t, ]
  summ[[paste(s, t)]] <- data.frame(feature_space = s, task = t, n_seeds = length(seeds), AUC_reported_seed12345 = ref$test_AUC,
    AUC_min = min(r[, "auc"]), AUC_median = median(r[, "auc"]), AUC_max = max(r[, "auc"]), fraction_AUC_ge_0.95 = mean(r[, "auc"] >= 0.95),
    nonzero_reported = ref$n_nonzero, nonzero_min = min(r[, "nonzero"]), nonzero_median = median(r[, "nonzero"]), nonzero_max = max(r[, "nonzero"]),
    fits_not_converged = sum(r[, "jerr"] != 0))
  cat(sprintf("%-22s %-20s seeds %3d | AUC reported %.4f, range %.4f-%.4f (median %.4f), %.1f%% >= 0.95 | non-zero %d-%d | non-converged %d\n",
              s, t, length(seeds), ref$test_AUC, min(r[, "auc"]), max(r[, "auc"]), median(r[, "auc"]), 100 * mean(r[, "auc"] >= 0.95),
              min(r[, "nonzero"]), max(r[, "nonzero"]), sum(r[, "jerr"] != 0)))
}
summ <- do.call(rbind, summ); write.csv(summ, file.path(out, "seed_sensitivity_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, per_seed), file.path(out, "seed_sensitivity_per_seed.csv"), row.names = FALSE)
wt <- as.data.frame(table(warning = unlist(warn_log))); write.csv(wt, file.path(out, "seed_sensitivity_warnings.csv"), row.names = FALSE)
cat("\nwarnings recorded during", sum(summ$n_seeds), "CV runs:\n"); print(wt, row.names = FALSE)
save_session("08b_elastic_net_seed_sensitivity")
