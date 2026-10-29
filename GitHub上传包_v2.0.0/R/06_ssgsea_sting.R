# =============================================================================
# 06_ssgsea_sting.R  -  GSE131193 immune signatures and the composite STING score (correct GEO labels).
#  Gene sets: Charoentong et al., Cell Rep 2017;18:248-262, Table S6 (28 immune cell types), mapped to mouse
#  orthologs with babelgene (HCOP consensus, min_support = 3) -> reference/Charoentong28_mouse.gmt.
#  ssGSEA: GSVA::ssgseaParam(minSize = 5, alpha = 0.25, normalize = TRUE).
#  STING score: z-mean of Mb21d1, Tmem173, Irf3, Isg15, Cxcl10 (Trex1 is not in this matrix).
#  Comparisons: primary = ipsilateral tMCAO vs ipsilateral sham (6 vs 6); secondary = same at day 1 and at day 7
#  (3 vs 3; minimum attainable exact two-sided Wilcoxon p = 0.1); all tMCAO vs all sham (12 vs 12).
#  Group tests: exact Wilcoxon rank-sum (two-sided); Welch t-test reported alongside. BH over 28 cell types.
#  Correlations: Spearman over all 24 samples (label-independent) and within the primary 12 samples.
# =============================================================================
source("R/00_config.R")
suppressMessages(library(GSVA))
out <- res_dir("06_ssgsea_sting")
expr <- read_mat(file.path(DIR_PRO, "GSE131193_log2counts.csv"))
ann  <- read.csv(file.path(DIR_PRO, "GSE131193_sample_annotation.csv")); expr <- expr[, ann$ID]
gmt <- strsplit(readLines(file.path(DIR_REF, "Charoentong28_mouse.gmt")), "\t")
gs <- setNames(lapply(gmt, function(x) intersect(x[-(1:2)], rownames(expr))), gsub("_", " ", sapply(gmt, `[`, 1)))
cat("immune gene sets:", length(gs), "| sizes in matrix:", paste(range(lengths(gs)), collapse = "-"), "\n")
ss <- gsva(ssgseaParam(expr, gs, minSize = 5, alpha = 0.25, normalize = TRUE), verbose = FALSE)
write.csv(t(ss), file.path(out, "ssGSEA_scores.csv"))

sting5 <- c("Mb21d1", "Tmem173", "Irf3", "Isg15", "Cxcl10"); stopifnot(all(sting5 %in% rownames(expr)))
cmp <- list(
  primary_ipsi_tMCAO_vs_ipsi_sham = with(ann, list(ID[injury == "tMCAO" & side == "ipsilateral"], ID[injury == "Sham" & side == "ipsilateral"])),
  day1_ipsi_tMCAO_vs_ipsi_sham    = with(ann, list(ID[injury == "tMCAO" & side == "ipsilateral" & time == "1d"], ID[injury == "Sham" & side == "ipsilateral" & time == "1d"])),
  day7_ipsi_tMCAO_vs_ipsi_sham    = with(ann, list(ID[injury == "tMCAO" & side == "ipsilateral" & time == "7d"], ID[injury == "Sham" & side == "ipsilateral" & time == "7d"])),
  all_tMCAO_vs_all_sham           = with(ann, list(ID[injury == "tMCAO"], ID[injury == "Sham"])))
test <- function(x, a, b) c(diff = mean(x[a]) - mean(x[b]), wilcox_p = wilcox.test(x[a], x[b], exact = TRUE)$p.value, welch_p = t.test(x[a], x[b])$p.value)

genes_of_interest <- c(sting5, "Stat1", "Tap1", "Tap2", "Tapbp", "B2m", "H2-K1", "H2-D1", "Psmb8", "Psmb10", "Irgm1", "Igtp", "Irf7", "Hmox1", "Timp1")
rows <- list(); cells <- list()
for (nm in names(cmp)) {
  a <- cmp[[nm]][[1]]; b <- cmp[[nm]][[2]]; cols <- c(a, b)
  sc <- zmean(expr[sting5, cols])
  rows[[paste(nm, "score")]] <- data.frame(comparison = nm, n_case = length(a), n_ctrl = length(b), feature = "STING5_zmean", t(test(sc, a, b)))
  for (g in intersect(genes_of_interest, rownames(expr))) rows[[paste(nm, g)]] <- data.frame(comparison = nm, n_case = length(a), n_ctrl = length(b), feature = g, t(test(expr[g, ], a, b)))
  d <- t(sapply(rownames(ss), function(k) test(ss[k, ], a, b)))
  cells[[nm]] <- data.frame(comparison = nm, cell = rownames(d), d, FDR_wilcox = p.adjust(d[, "wilcox_p"], "BH"))
}
gt <- do.call(rbind, rows); ct <- do.call(rbind, cells)
write.csv(gt, file.path(out, "genes_and_STING_score_group_tests.csv"), row.names = FALSE)
write.csv(ct, file.path(out, "ssGSEA_group_tests.csv"), row.names = FALSE)
cat("\n== STING score and genes ==\n"); print(transform(gt, diff = round(diff, 3), wilcox_p = signif(wilcox_p, 3), welch_p = signif(welch_p, 3)), row.names = FALSE)
for (nm in names(cmp)) { x <- ct[ct$comparison == nm, ]; x <- x[order(x$wilcox_p), ][1:6, ]
  cat("\n== ssGSEA top 6:", nm, "==\n"); print(transform(x[, c("cell", "diff", "wilcox_p", "FDR_wilcox")], diff = round(diff, 4), wilcox_p = signif(wilcox_p, 3), FDR_wilcox = signif(FDR_wilcox, 3)), row.names = FALSE) }

variants <- list(STING5 = sting5, STING3_noISG = c("Mb21d1", "Tmem173", "Irf3"), STING4_noIsg15 = c("Mb21d1", "Tmem173", "Irf3", "Cxcl10"), STING4_noCxcl10 = c("Mb21d1", "Tmem173", "Irf3", "Isg15"))
cr <- list()
for (set in c("all24", "primary12")) {
  cols <- if (set == "all24") ann$ID else unlist(cmp$primary_ipsi_tMCAO_vs_ipsi_sham)
  for (v in names(variants)) {
    sc <- zmean(expr[variants[[v]], cols])
    r <- t(sapply(rownames(ss), function(k) { z <- cor.test(sc, ss[k, cols], method = "spearman", exact = FALSE); c(rho = unname(z$estimate), p = z$p.value) }))
    cr[[paste(set, v)]] <- data.frame(samples = set, score = v, cell = rownames(r), rho = r[, "rho"], p = r[, "p"], FDR = p.adjust(r[, "p"], "BH"))
  }
}
cr <- do.call(rbind, cr); write.csv(cr, file.path(out, "STING_score_vs_ssGSEA_spearman.csv"), row.names = FALSE)
cat("\n== STING5 vs ssGSEA, top 6 by p ==\n")
for (set in c("all24", "primary12")) { x <- cr[cr$samples == set & cr$score == "STING5", ]; x <- x[order(x$p), ][1:6, ]
  cat("--", set, "\n"); print(transform(x[, c("cell", "rho", "p", "FDR")], rho = round(rho, 3), p = signif(p, 3), FDR = signif(FDR, 3)), row.names = FALSE) }
cat("\n== Activated CD8 T cell, all score variants ==\n")
x <- cr[cr$cell == "Activated CD8 T cell", ]; print(transform(x[, c("samples", "score", "rho", "FDR")], rho = round(rho, 3), FDR = signif(FDR, 3)), row.names = FALSE)
cat("\nSTING-score genes contained in any Charoentong set:", paste(intersect(unique(unlist(gs)), c(sting5, "Stat1")), collapse = ",") , "(empty = none)\n")
save_session("06_ssgsea_sting")
