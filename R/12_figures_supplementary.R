# =============================================================================
# 12_figures_supplementary.R  -  Supplementary Figures S1-S3, drawn only from files in data/processed and results/.
#  S1: batch-correction sensitivity (merged cohort): DEG overlap across strategies and per-dataset concordance.
#  S2: Elastic Net (GSE23163) test-set ROC curves for the three feature spaces and both hold-out tasks.
#  S3: GSE227317 rat cortex (5 MCAO vs 5 Sham): panel genes and panel score.
# =============================================================================
source("R/00_config.R")
suppressMessages({ library(ggplot2); library(patchwork); library(pROC) })
FIG <- "figures"; dir.create(FIG, showWarnings = FALSE)
th <- theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 8.5, face = "bold"),
                                     legend.key.size = unit(3, "mm"), strip.background = element_rect(fill = "grey92"))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, units = "in", device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".tiff")), p, width = w, height = h, units = "in", dpi = 300, device = ragg::agg_tiff, compression = "lzw")
  cat("saved", name, "\n")
}
fmtp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))
rp <- function(step, f) file.path(DIR_RES, step, f)

## ---------------- Figure S1: batch-correction sensitivity ----------------
A <- read.csv(rp("03_combat_limma", "limma_primary_ComBat_withGroup_all_genes.csv"), row.names = 1)
B <- read.csv(rp("03_combat_limma", "limma_sensitivity_B.csv"), row.names = 1)
C <- read.csv(rp("03_combat_limma", "limma_sensitivity_C.csv"), row.names = 1)
E1 <- read.csv(rp("03_combat_limma", "limma_within_GSE131712.csv"), row.names = 1)
E2 <- read.csv(rp("03_combat_limma", "limma_within_GSE202659.csv"), row.names = 1)
g <- rownames(A)
d1 <- data.frame(A = A[g, "logFC"], B = B[g, "logFC"], C = C[g, "logFC"], deg = A[g, "DEG"])
r_ab <- cor(d1$A, d1$B); r_ac <- cor(d1$A, d1$C)
s1a <- ggplot(d1, aes(A, B)) + geom_point(aes(colour = deg), size = 0.3, alpha = 0.5) + geom_abline(linetype = 2, linewidth = 0.3) +
  scale_colour_manual(values = c(`TRUE` = "#E45756", `FALSE` = "grey70"), labels = c(`TRUE` = "primary DEG", `FALSE` = "other"), name = NULL) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 2.3, label = sprintf("Pearson r = %.4f (%d genes)", r_ab, length(g))) +
  th + labs(x = "log2FC, ComBat with group covariate", y = "log2FC, ComBat without covariate", title = "Primary vs ComBat without group covariate")
s1b <- ggplot(d1, aes(A, C)) + geom_point(aes(colour = deg), size = 0.3, alpha = 0.5) + geom_abline(linetype = 2, linewidth = 0.3) +
  scale_colour_manual(values = c(`TRUE` = "#E45756", `FALSE` = "grey70"), guide = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 2.3, label = sprintf("Pearson r = %.4f (%d genes)", r_ac, length(g))) +
  th + labs(x = "log2FC, ComBat with group covariate", y = "log2FC, no ComBat (batch in limma model)", title = "Primary vs batch as model covariate")
dd <- A$DEG
d3 <- data.frame(E1 = E1[g, "logFC"][dd], E2 = E2[g, "logFC"][dd], up = A[g, "logFC"][dd] > 0)
same <- sum(sign(d3$E1) == sign(A[g, "logFC"][dd]) & sign(d3$E2) == sign(A[g, "logFC"][dd]))
s1c <- ggplot(d3, aes(E1, E2, colour = up)) + geom_hline(yintercept = 0, linewidth = 0.3) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(size = 0.5, alpha = 0.6) + scale_colour_manual(values = c(`TRUE` = "#E45756", `FALSE` = "#4C78A8"), labels = c(`TRUE` = "up in primary", `FALSE` = "down in primary"), name = NULL) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 2.3, label = sprintf("%d of %d DEGs: same sign\nin both datasets", same, sum(dd))) +
  th + labs(x = "log2FC within GSE131712 (3 vs 3)", y = "log2FC within GSE202659 (3 vs 3)", title = "Per-dataset direction of the primary DEGs")
ss <- read.csv(rp("03_combat_limma", "R2.5_batch_sensitivity_summary.csv"))
ss$lab <- c("ComBat, group covariate\n(primary)", "ComBat, no covariate", "No ComBat,\nbatch in model", "No correction", "GSE131712 alone", "GSE202659 alone")
sl <- rbind(data.frame(lab = ss$lab, what = "DEGs", n = ss$n_DEG), data.frame(lab = ss$lab, what = "shared with primary", n = ss$shared_with_primary))
sl$lab <- factor(sl$lab, levels = rev(ss$lab))
s1d <- ggplot(sl, aes(n, lab, fill = what)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(0.75), hjust = -0.15, size = 2.1) + scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c(DEGs = "grey55", `shared with primary` = "#E45756"), name = NULL) + th +
  labs(x = "Genes (FDR < 0.05, |log2FC| > 0.585)", y = NULL, title = "DEG counts across batch-handling strategies")
save_fig((s1a | s1b) / (s1c | s1d) + plot_annotation(tag_levels = "A"), "Figure_S1", 7.2, 6.4)

## ---------------- Figure S2: Elastic Net ROC curves ----------------
en <- read.csv(rp("08_elastic_net", "elastic_net_summary.csv"))
lab_fs <- c(A_original_39 = "Original 39-gene panel", B_full_transcriptome = "Full transcriptome (7,693 genes)", C_rebuilt_panel = "Rebuilt panel (corrected PPI hubs)")
lab_t <- c(Cortex_to_Striatum = "Train cortex \u2192 test striatum", WT_to_Gpx1KO = "Train WT \u2192 test Gpx1-KO")
rc <- list(); ann <- list()
for (i in seq_len(nrow(en))) {
  pr <- read.csv(rp("08_elastic_net", sprintf("test_predictions_%s_%s.csv", en$feature_space[i], en$task[i])))
  r <- roc(pr$y, pr$prob, levels = c(0, 1), direction = "<", quiet = TRUE)
  stopifnot(abs(as.numeric(auc(r)) - en$test_AUC[i]) < 1e-6)
  rc[[i]] <- data.frame(fs = unname(lab_fs[en$feature_space[i]]), task = unname(lab_t[en$task[i]]), fpr = 1 - r$specificities, tpr = r$sensitivities)
  ann[[i]] <- data.frame(fs = unname(lab_fs[en$feature_space[i]]), task = unname(lab_t[en$task[i]]),
                         lab = sprintf("AUC %.3f (95%% CI %.3f\u2013%.3f)\n%d of %s genes retained", en$test_AUC[i], en$CI_low[i], en$CI_high[i], en$n_nonzero[i], format(en$n_candidates[i], big.mark = ",")))
}
rc <- do.call(rbind, rc); ann <- do.call(rbind, ann)
rc$fs <- factor(rc$fs, levels = lab_fs); ann$fs <- factor(ann$fs, levels = lab_fs)
rc <- rc[order(rc$fs, rc$task, rc$fpr, rc$tpr), ]
s2 <- ggplot(rc, aes(fpr, tpr)) + geom_abline(linetype = 2, linewidth = 0.3, colour = "grey50") + geom_step(colour = "#E45756", linewidth = 0.6, direction = "vh") +
  geom_text(data = ann, aes(x = 1, y = 0.05, label = lab), hjust = 1, vjust = 0, size = 2.2, lineheight = 0.9) +
  facet_grid(task ~ fs) + coord_equal() + th + labs(x = "1 \u2212 specificity", y = "Sensitivity",
  title = "Elastic Net (alpha 0.4, lambda.min, 10-fold CV) in GSE23163: held-out test sets (8 Sham vs 8 I/R 24 h each)")
save_fig(s2, "Figure_S2", 7.2, 5.2)

## ---------------- Figure S3: rat GSE227317 ----------------
fp <- read_mat(rp("09_rat_roc", "gene_FPKM.csv")); rr <- read.csv(rp("09_rat_roc", "rat_ROC.csv"))
soft <- readLines(gzfile(file.path(DIR_RAW, "GSE227317_family.soft.gz")), warn = FALSE)
ttl <- sub("^!Sample_title = ", "", grep("^!Sample_title", soft, value = TRUE))
lab <- data.frame(sample = sub(".*\\[(T\\d+)\\]$", "\\1", ttl), Group = ifelse(grepl("MCAO", ttl), "MCAO", "Sham"))
lg <- log2(fp[, lab$sample] + 1)
genes <- c("Tap1", "Tap2", "Cxcl10", "Psmb8", "Stat1", "Psmb10", "Trpc6")
pl <- do.call(rbind, lapply(genes, function(g) data.frame(feature = g, sample = lab$sample, Group = lab$Group, value = lg[g, ])))
panel <- zmean(lg[c("Tap1", "Tap2", "Cxcl10", "Psmb8", "Stat1"), ])
pl <- rbind(data.frame(feature = "Panel score\n(z-mean of 5 genes)", sample = lab$sample, Group = lab$Group, value = panel[lab$sample]), pl)
pl$feature <- factor(pl$feature, levels = unique(pl$feature)); pl$Group <- factor(pl$Group, levels = c("Sham", "MCAO"))
rr$feature <- ifelse(rr$feature == "Panel5_zmean", "Panel score\n(z-mean of 5 genes)", rr$feature)
lb <- data.frame(feature = factor(rr$feature, levels = levels(pl$feature)), lab = sprintf("AUC %.2f\np = %s", rr$AUC, fmtp(rr$wilcox_exact_p)))
s3 <- ggplot(pl, aes(Group, value, fill = Group)) + geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, alpha = 0.6) +
  geom_point(position = position_jitter(width = 0.08, seed = 1), size = 1.2) +
  geom_text(data = lb, aes(x = 1.5, y = Inf, label = lab), vjust = 1.2, size = 2.1, lineheight = 0.9, inherit.aes = FALSE) +
  facet_wrap(~feature, scales = "free_y", nrow = 2) + scale_fill_manual(values = c(Sham = "#4C78A8", MCAO = "#E45756"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) + th + theme(strip.text = element_text(face = "italic")) +
  labs(x = NULL, y = "log2(FPKM + 1) (panel: z-mean)", title = "GSE227317 rat cortex, MCAO 2 h/R 24 h vs Sham (n = 5 vs 5); exploratory directional check",
       caption = "AUC from the expression value of each gene; with 5 vs 5 samples an AUC of 1 has a saturated bootstrap CI. p: exact Wilcoxon rank-sum test.")
save_fig(s3, "Figure_S3", 7.2, 4.4)
