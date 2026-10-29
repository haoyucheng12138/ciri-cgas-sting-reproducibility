# =============================================================================
# 08_elastic_net.R  -  Elastic Net logistic regression on GSE23163 (glmnet).
#  Samples: Sham vs I/R 24 h (other time points excluded). Two hold-out tasks:
#    (1) train cortex (8 Sham + 8 I/R-24h, WT and Gpx1-/-) -> test striatum (16)
#    (2) train WT (cortex + striatum, 16) -> test Gpx1-/- (16)
#  Model: alpha = 0.4; lambda.min by 10-fold CV (binomial deviance); predictors standardised internally;
#  set.seed(12345) immediately before cv.glmnet. Test AUC with 2,000 stratified bootstrap 95% CI (pROC).
#  Feature spaces:
#    A. original 39-gene panel (reference/original_39gene_candidate_panel.csv; B2m not on the array) - reproduces
#       the original model; coefficients are reported from this fit (replace the earlier Table S3/S4 values).
#    B. full transcriptome (all 7,693 genes) - Reviewer 1.1.
#    C. rebuilt panel: the original mechanistic genes + Trpc6 + top-20 hubs of the corrected PPI (04).
# =============================================================================
source("R/00_config.R")
suppressMessages({ library(glmnet); library(pROC) })
out <- res_dir("08_elastic_net")
x <- read_mat(file.path(DIR_PRO, "GSE23163_log2ratio.csv"))
meta <- read.csv(file.path(DIR_PRO, "GSE23163_sample_annotation.csv"))
meta <- meta[meta$group %in% c("Sham", "IR24h"), ]; meta$y <- as.integer(meta$group == "IR24h")
tasks <- list(Cortex_to_Striatum = list(train = meta$ID[meta$region == "Cortex"], test = meta$ID[meta$region == "Striatum"]),
              WT_to_Gpx1KO = list(train = meta$ID[meta$genotype == "WT"], test = meta$ID[meta$genotype == "Gpx1KO"]))
for (t in names(tasks)) cat(sprintf("%s: train %d (%d I/R) | test %d (%d I/R)\n", t, length(tasks[[t]]$train), sum(meta$y[meta$ID %in% tasks[[t]]$train]),
                                    length(tasks[[t]]$test), sum(meta$y[meta$ID %in% tasks[[t]]$test])))
y_of <- function(ids) meta$y[match(ids, meta$ID)]

panel39 <- setdiff(read.csv(file.path(DIR_REF, "original_39gene_candidate_panel.csv"))$Gene, "Gene")
old_hub <- c("B2m", "Tapbp", "H2-K1", "Eif6", "H2-T23", "Tap1", "H2-D1", "Tap2", "Rpl7", "H2-Q4", "Eif4a1", "Ifi35", "Cdh5", "Rplp0", "Rpl18", "Il10", "Ccl5", "Cd163", "Lama5", "Cd274")
mech <- setdiff(panel39, c(old_hub, "Trpc6"))
new_hub <- head(read.csv(file.path(res_dir("04_enrichment_ppi"), "PPI_nodes_score700.csv"))$Gene, 20)
spaces <- list(A_original_39 = panel39, B_full_transcriptome = rownames(x), C_rebuilt_panel = unique(c(mech, new_hub, "Trpc6")))
writeLines(c("# mechanistic genes (unchanged):", mech, "# corrected PPI top-20 hubs:", new_hub, "# a priori:", "Trpc6"), file.path(out, "rebuilt_panel_definition.txt"))

fit_one <- function(genes, tr, te) {
  genes <- intersect(genes, rownames(x)); genes <- genes[apply(x[genes, tr, drop = FALSE], 1, sd) > 0]
  X <- t(x[genes, tr]); y <- y_of(tr)
  set.seed(SEED)
  cv <- cv.glmnet(X, y, family = "binomial", alpha = 0.4, nfolds = 10, type.measure = "deviance", standardize = TRUE)
  co <- as.matrix(coef(cv, s = "lambda.min"))[, 1]
  p <- as.numeric(predict(cv, newx = t(x[genes, te]), s = "lambda.min", type = "response"))
  r <- roc(y_of(te), p, levels = c(0, 1), direction = "<", quiet = TRUE)
  set.seed(SEED); ci <- ci.auc(r, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE, progress = "none")
  nz <- co[-1][co[-1] != 0]; nz <- nz[order(-abs(nz))]
  list(n_cand = length(genes), lambda = cv$lambda.min, intercept = co[1], nz = nz, auc = as.numeric(auc(r)), ci = c(ci[1], ci[3]),
       acc = mean((p >= 0.5) == y_of(te)), pred = data.frame(ID = te, y = y_of(te), prob = p))
}
summ <- list()
for (s in names(spaces)) for (t in names(tasks)) {
  f <- fit_one(spaces[[s]], tasks[[t]]$train, tasks[[t]]$test)
  write.csv(data.frame(rank = seq_along(f$nz), gene = names(f$nz), coefficient = round(f$nz, 6)), file.path(out, sprintf("coefficients_%s_%s.csv", s, t)), row.names = FALSE)
  write.csv(f$pred, file.path(out, sprintf("test_predictions_%s_%s.csv", s, t)), row.names = FALSE)
  ap <- c("Tap1", "Tap2", "Tapbp", "Tapbpl", "H2-K1", "H2-D1", "H2-T23", "Psmb8", "Psmb9", "Psmb10", "B2m")
  rk <- match(ap, names(f$nz))
  summ[[paste(s, t)]] <- data.frame(feature_space = s, task = t, n_candidates = f$n_cand, lambda_min = signif(f$lambda, 4), intercept = round(f$intercept, 4),
    n_nonzero = length(f$nz), test_AUC = round(f$auc, 4), CI_low = round(f$ci[1], 3), CI_high = round(f$ci[2], 3), accuracy = round(f$acc, 4),
    top10_by_abs_coef = paste(head(names(f$nz), 10), collapse = ";"),
    antigen_genes_rank = paste(sprintf("%s=%s", ap[!is.na(rk)], rk[!is.na(rk)]), collapse = ";"))
}
summ <- do.call(rbind, summ); write.csv(summ, file.path(out, "elastic_net_summary.csv"), row.names = FALSE)
print(summ[, c("feature_space", "task", "n_candidates", "n_nonzero", "test_AUC", "CI_low", "CI_high", "accuracy")], row.names = FALSE)
cat("\n"); for (i in seq_len(nrow(summ))) cat(sprintf("%s | %s\n  top10: %s\n  antigen-processing ranks: %s\n", summ$feature_space[i], summ$task[i], summ$top10_by_abs_coef[i], summ$antigen_genes_rank[i]))
save_session("08_elastic_net")
