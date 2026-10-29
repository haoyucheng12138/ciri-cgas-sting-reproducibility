# =============================================================================
# 03_combat_limma.R  -  merged validation cohort (6 Sham vs 6 tMCAO).
#  Primary (as in the original analysis): ComBat(batch = dataset, mod = ~group) -> limma ~group.
#  Sensitivity (Reviewer 2.5): (B) ComBat without the group covariate; (C) no ComBat, limma ~batch + group;
#  (D) no correction; (E) each dataset analysed alone (3 vs 3).
#  DEG: BH-FDR < 0.05 and |log2FC| > 0.585.
#  STING composite score: z-mean of Mb21d1, Tmem173, Irf3, Isg15, Cxcl10, Trex1 on the ComBat(mod = NULL) matrix,
#  Welch t-test (this is the definition behind the originally reported p = 4.79e-4).
# =============================================================================
source("R/00_config.R")
suppressMessages({ library(sva); library(limma) })
out <- res_dir("03_combat_limma")
m  <- read_mat(file.path(DIR_PRO, "merged_log2FPKM.csv"))
s2 <- read.csv(file.path(DIR_REF, "TableS2_merged_cohort_samples.csv"))
m <- m[, s2$Sample_ID]
group <- factor(s2$Group, levels = c("Sham", "Stroke")); batch <- factor(s2$Dataset)
d1 <- model.matrix(~ group)
cb_mod   <- ComBat(m, batch = batch, mod = d1, par.prior = TRUE)
cb_nomod <- ComBat(m, batch = batch, mod = NULL, par.prior = TRUE)
write.csv(cb_mod, file.path(DIR_PRO, "merged_ComBat_withGroup.csv"))
write.csv(cb_nomod, file.path(DIR_PRO, "merged_ComBat_noGroup.csv"))

fit <- function(x, design, coef) {
  tt <- topTable(eBayes(lmFit(x, design)), coef = coef, number = Inf, sort.by = "none")
  tt$DEG <- tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.585; tt
}
A <- fit(cb_mod, d1, 2)
B <- fit(cb_nomod, d1, 2)
C <- fit(m, model.matrix(~ batch + group), "groupStroke")
D <- fit(m, d1, 2)
E <- lapply(levels(batch), function(b) { i <- batch == b; fit(m[, i], model.matrix(~ droplevels(group[i])), 2) }); names(E) <- levels(batch)
write.csv(A, file.path(out, "limma_primary_ComBat_withGroup_all_genes.csv"))
write.csv(A[A$DEG, ], file.path(out, "TableS6_significant_DEGs.csv"))
for (nm in c("B", "C", "D")) write.csv(get(nm), file.path(out, paste0("limma_sensitivity_", nm, ".csv")))
for (b in names(E)) write.csv(E[[b]], file.path(out, paste0("limma_within_", b, ".csv")))
summ <- do.call(rbind, lapply(list(A = A, B = B, C = C, D = D, E1 = E[[1]], E2 = E[[2]]), function(x)
  data.frame(n_DEG = sum(x$DEG), up = sum(x$DEG & x$logFC > 0), down = sum(x$DEG & x$logFC < 0), shared_with_primary = sum(x$DEG & A$DEG))))
summ$strategy <- c("A ComBat(mod=~group) [primary]", "B ComBat(mod=NULL)", "C no ComBat, ~batch+group", "D no correction",
                   paste("E within", names(E)))
same_dir <- sign(E[[1]]$logFC) == sign(A$logFC) & sign(E[[2]]$logFC) == sign(A$logFC)
summ$note <- ""; summ$note[1] <- sprintf("%d of %d primary DEGs have the same sign in both datasets separately", sum(same_dir & A$DEG), sum(A$DEG))
write.csv(summ, file.path(out, "R2.5_batch_sensitivity_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)

key <- c("Cxcl10","Isg15","Tmem173","Mb21d1","Irf3","Trex1","Stat1","Il6","Trpc6","Tap1","Tap2","Tapbp","Tapbpl","B2m","H2-K1","H2-D1","Psmb8","Psmb9","Psmb10","Irgm1","Irgm2","Igtp","Trim56")
kt <- data.frame(gene = key)
for (nm in c("A", "B", "C")) { x <- get(nm); kt[[paste0(nm, "_logFC")]] <- x[key, "logFC"]; kt[[paste0(nm, "_FDR")]] <- x[key, "adj.P.Val"] }
for (b in names(E)) { kt[[paste0(b, "_logFC")]] <- E[[b]][key, "logFC"]; kt[[paste0(b, "_P")]] <- E[[b]][key, "P.Value"] }
write.csv(kt, file.path(out, "R2.5_key_genes_across_strategies.csv"), row.names = FALSE)

sting6 <- c("Mb21d1", "Tmem173", "Irf3", "Isg15", "Cxcl10", "Trex1")
sc <- zmean(cb_nomod[sting6, ]); st <- group == "Stroke"
tt <- t.test(sc[st], sc[!st]); wt <- wilcox.test(sc[st], sc[!st], exact = TRUE)
res <- data.frame(score = "STING6_zmean (ComBat mod=NULL)", diff = mean(sc[st]) - mean(sc[!st]), welch_p = tt$p.value, wilcoxon_exact_p = wt$p.value)
write.csv(res, file.path(out, "STING_score_merged_cohort.csv"), row.names = FALSE); print(res, row.names = FALSE)
save_session("03_combat_limma")
