# =============================================================================
# 11_figures_main.R  -  Figures 2-6 of the revised manuscript, drawn only from files in data/processed and results/.
# Every number printed on a figure is read from those files (no typed values). Output: figures/*.pdf and 300-dpi *.tiff.
# =============================================================================
source("R/00_config.R")
suppressMessages({ library(ggplot2); library(patchwork); library(ggrepel); library(igraph); library(scales) })
FIG <- "figures"; dir.create(FIG, showWarnings = FALSE)
th <- theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 8.5, face = "bold"),
                                     legend.key.size = unit(3, "mm"), strip.background = element_rect(fill = "grey92"))
COL_GRP <- c(Sham = "#4C78A8", Stroke = "#E45756", tMCAO = "#E45756")
COL_BAT <- c(GSE131712 = "#54A24B", GSE202659 = "#B279A2")
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(FIG, paste0(name, ".tiff")), p, width = w, height = h, units = "in", dpi = 300, device = ragg::agg_tiff, compression = "lzw", bg = "white")
  cat("saved", name, "\n")
}
fmtp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))
rp <- function(step, f) file.path(DIR_RES, step, f)

## ---------------- Figure 2: merged cohort, batch correction ----------------
s2 <- read.csv(file.path(DIR_REF, "TableS2_merged_cohort_samples.csv"))
raw <- read_mat(file.path(DIR_PRO, "merged_log2FPKM.csv"))[, s2$Sample_ID]
cb  <- read_mat(file.path(DIR_PRO, "merged_ComBat_withGroup.csv"))[, s2$Sample_ID]
pca_df <- function(m, lab) {
  m <- m[apply(m, 1, sd) > 0, ]; p <- prcomp(t(m), center = TRUE, scale. = FALSE); ve <- 100 * p$sdev^2 / sum(p$sdev^2)
  data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2], Group = s2$Group, Dataset = s2$Dataset, panel = sprintf("%s (PC1 %.1f%%, PC2 %.1f%%)", lab, ve[1], ve[2]))
}
pd <- rbind(pca_df(raw, "Before ComBat"), pca_df(cb, "After ComBat"))
pd$panel <- factor(pd$panel, levels = unique(pd$panel))
f2a <- ggplot(pd, aes(PC1, PC2, colour = Dataset, shape = Group)) + geom_point(size = 2.2) + facet_wrap(~panel, scales = "free") +
  scale_colour_manual(values = COL_BAT) + scale_shape_manual(values = c(Sham = 16, Stroke = 17)) + th + labs(title = "PCA of all 20,078 genes")
bx <- rbind(data.frame(stage = "Before ComBat", sample = rep(colnames(raw), each = nrow(raw)), value = as.vector(raw)),
            data.frame(stage = "After ComBat", sample = rep(colnames(cb), each = nrow(cb)), value = as.vector(cb)))
bx$Dataset <- s2$Dataset[match(bx$sample, s2$Sample_ID)]; bx$stage <- factor(bx$stage, levels = c("Before ComBat", "After ComBat"))
bx$sample <- factor(bx$sample, levels = s2$Sample_ID)
f2b <- ggplot(bx, aes(sample, value, fill = Dataset)) + geom_boxplot(outlier.size = 0.1, linewidth = 0.25) + facet_wrap(~stage, ncol = 1) +
  scale_fill_manual(values = COL_BAT) + th + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  labs(x = NULL, y = "log2(FPKM + 1)", title = "Sample expression distributions")
save_fig((f2a / f2b) + plot_layout(heights = c(1, 1.3)) + plot_annotation(tag_levels = "A"), "Figure_2", 7.2, 7)

## ---------------- Figure 3: merged cohort, DEGs ----------------
A <- read.csv(rp("03_combat_limma", "limma_primary_ComBat_withGroup_all_genes.csv"), row.names = 1)
A$cls <- ifelse(A$DEG & A$logFC > 0, "Up", ifelse(A$DEG & A$logFC < 0, "Down", "NS"))
nU <- sum(A$cls == "Up"); nD <- sum(A$cls == "Down")
lab_genes <- c("Cxcl10", "Isg15", "Tmem173", "Stat1", "Tap1", "Tap2", "H2-K1", "H2-D1", "Psmb8", "Irgm1", "Igtp", "Mb21d1", "Irf3")
A$lab <- ifelse(rownames(A) %in% lab_genes, rownames(A), NA)
f3a <- ggplot(A, aes(logFC, -log10(adj.P.Val), colour = cls)) + geom_point(size = 0.4, alpha = 0.6) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = 2, linewidth = 0.3) + geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
  geom_text_repel(aes(label = lab), size = 2.3, colour = "black", min.segment.length = 0, max.overlaps = 50, na.rm = TRUE) +
  scale_colour_manual(values = c(Up = "#E45756", Down = "#4C78A8", NS = "grey75"), labels = c(Up = paste0("Up (", nU, ")"), Down = paste0("Down (", nD, ")"), NS = "NS")) +
  th + labs(x = "log2 fold change (tMCAO vs Sham)", y = "-log10 FDR", colour = NULL, title = sprintf("Merged cohort: %d DEGs (FDR < 0.05, |log2FC| > 0.585)", nU + nD))
kt <- read.csv(rp("03_combat_limma", "R2.5_key_genes_across_strategies.csv"))
show <- c("Cxcl10", "Isg15", "Tmem173", "Mb21d1", "Irf3", "Stat1", "Tap1", "Tap2", "Tapbp", "B2m", "H2-K1", "H2-D1", "Psmb8", "Psmb10")
kl <- do.call(rbind, lapply(c("A", "B", "C"), function(s) data.frame(gene = kt$gene, strategy = s, logFC = kt[[paste0(s, "_logFC")]], FDR = kt[[paste0(s, "_FDR")]])))
kl <- kl[kl$gene %in% show, ]; kl$gene <- factor(kl$gene, levels = rev(show))
kl$strategy <- factor(kl$strategy, levels = c("A", "B", "C"), labels = c("ComBat\n(group kept)", "ComBat\n(no covariate)", "No ComBat\n(batch in model)"))
kl$star <- ifelse(kl$FDR < 0.001, "***", ifelse(kl$FDR < 0.01, "**", ifelse(kl$FDR < 0.05, "*", "")))
f3b <- ggplot(kl, aes(strategy, gene, fill = logFC)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f%s", logFC, star)), size = 2.1) +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#E45756", midpoint = 0) + th + labs(x = NULL, y = NULL, title = "Key genes: log2FC under three batch strategies") +
  theme(axis.text.y = element_text(face = "italic"))
ss <- read.csv(rp("03_combat_limma", "R2.5_batch_sensitivity_summary.csv"))
ss <- ss[1:3, ]; ss$lab <- c("ComBat\n(group kept)", "ComBat\n(no covariate)", "No ComBat\n(batch in model)")
sl <- rbind(data.frame(lab = ss$lab, dir = "Up", n = ss$up), data.frame(lab = ss$lab, dir = "Down", n = ss$down)); sl$lab <- factor(sl$lab, levels = ss$lab)
f3c <- ggplot(sl, aes(lab, n, fill = dir)) + geom_col(width = 0.6) + geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 2.3) +
  scale_fill_manual(values = c(Up = "#E45756", Down = "#4C78A8")) + th + labs(x = NULL, y = "DEGs", fill = NULL, title = "DEG counts by batch strategy")
save_fig(f3a + (f3b / f3c + plot_layout(heights = c(1.6, 1))) + plot_layout(widths = c(1.2, 1)) + plot_annotation(tag_levels = "A"), "Figure_3", 7.2, 5.6)

## ---------------- Figure 4: WGCNA, GSE131193 ----------------
sf <- read.csv(rp("05_wgcna", "soft_threshold.csv")); pw <- sf$Power[which(sf$signedR2 >= 0.85)[1]]
f4a <- ggplot(sf, aes(Power, signedR2)) + geom_hline(yintercept = 0.85, linetype = 2, colour = "red", linewidth = 0.3) + geom_line(linewidth = 0.3) +
  geom_point(aes(colour = Power == pw), size = 1.4, show.legend = FALSE) + scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = "black")) +
  annotate("text", x = pw, y = sf$signedR2[sf$Power == pw] - 0.08, label = sprintf("power %d\nR\u00b2 = %.3f", pw, sf$signedR2[sf$Power == pw]), size = 2.3, colour = "red") +
  th + labs(x = "Soft-thresholding power", y = "Scale-free fit (signed R\u00b2)", title = "Soft threshold (unsigned network, n = 24)")
mt <- read.csv(rp("05_wgcna", "module_trait_r.csv"), check.names = FALSE); mp <- read.csv(rp("05_wgcna", "module_trait_p.csv"), check.names = FALSE)
traits <- c(ipsi_tMCAO_vs_ipsi_sham = "tMCAO vs sham\n(ipsilateral)", tMCAO_all = "tMCAO vs sham\n(all 24)", ipsi_vs_contra_in_tMCAO = "Ipsi vs contra\n(tMCAO)", day7 = "Day 7 vs day 1")
ml <- do.call(rbind, lapply(names(traits), function(t) data.frame(module = mt$module, size = mt$size, trait = traits[[t]], r = mt[[t]], p = mp[[t]])))
ml <- ml[ml$module != "MEgrey", ]
ord <- mt$module[mt$module != "MEgrey"][order(mt$ipsi_tMCAO_vs_ipsi_sham[mt$module != "MEgrey"])]
ml$module <- factor(sprintf("%s (%d)", sub("^ME", "", ml$module), ml$size), levels = sprintf("%s (%d)", sub("^ME", "", ord), mt$size[match(ord, mt$module)]))
ml$trait <- factor(ml$trait, levels = traits)
ml$lab <- sprintf("%.2f\n(%s)", ml$r, fmtp(ml$p))
f4b <- ggplot(ml, aes(trait, module, fill = r)) + geom_tile(colour = "white") + geom_text(aes(label = lab), size = 1.7, lineheight = 0.85) +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#E45756", midpoint = 0, limits = c(-1, 1)) + th +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + labs(x = NULL, y = "Module (genes)", fill = "r", title = "Module eigengene\u2013trait correlation\n(Pearson r, p)")
fg <- read.csv(rp("05_wgcna", "focal_genes.csv"))
fg$gene <- factor(fg$gene, levels = fg$gene[order(fg$module_r_primary, fg$kME)])
fg$assoc <- ifelse(fg$module_p_primary < 0.05, "module associated (p < 0.05)", "module not associated")
## kME can be negative in an unsigned network (gene anti-correlated with its module eigengene, e.g. Mb21d1, Tapbpl):
## point size = |kME|, label shows the signed kME so no gene is dropped.
fg$lab <- sprintf("%s, kME %.2f", fg$module, fg$kME)
f4c <- ggplot(fg, aes(module_r_primary, gene, colour = assoc)) + geom_vline(xintercept = 0, linewidth = 0.3) + geom_point(aes(size = abs(kME))) +
  geom_text(aes(label = lab), hjust = -0.18, size = 1.9, colour = "grey30") + scale_colour_manual(values = c("module associated (p < 0.05)" = "#E45756", "module not associated" = "grey55")) +
  scale_size_continuous(range = c(0.8, 2.6), limits = c(0, 1)) + scale_x_continuous(limits = c(min(fg$module_r_primary) - 0.1, 1.45)) + th +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom", legend.box = "vertical") +
  labs(x = "Module r with ipsilateral tMCAO vs sham", y = NULL, colour = NULL, size = "|kME|", title = "Module of each axis gene\n(label: module, signed kME)")
save_fig(((f4a / f4c + plot_layout(heights = c(1, 2.6))) | f4b) + plot_layout(widths = c(1.1, 1)) + plot_annotation(tag_levels = "A"), "Figure_4", 7.2, 7.8)

## ---------------- Figure 5: GSE131193 STING axis and immune signatures ----------------
ann <- read.csv(file.path(DIR_PRO, "GSE131193_sample_annotation.csv"))
e131 <- read_mat(file.path(DIR_PRO, "GSE131193_log2counts.csv"))[, ann$ID]
ipsi <- ann[ann$side == "ipsilateral", ]; ipsi$Group <- factor(ifelse(ipsi$injury == "tMCAO", "tMCAO", "Sham"), levels = c("Sham", "tMCAO"))
gt <- read.csv(rp("06_ssgsea_sting", "genes_and_STING_score_group_tests.csv"))
gt <- gt[gt$comparison == "primary_ipsi_tMCAO_vs_ipsi_sham", ]
sting5 <- c("Mb21d1", "Tmem173", "Irf3", "Isg15", "Cxcl10")
sc <- zmean(e131[sting5, ipsi$ID])
pl <- rbind(data.frame(feature = "IFN/STING score", ID = ipsi$ID, value = sc),
            do.call(rbind, lapply(c(sting5, "Tap1", "H2-K1", "Psmb8"), function(g) data.frame(feature = g, ID = ipsi$ID, value = e131[g, ipsi$ID]))))
pl$Group <- ipsi$Group[match(pl$ID, ipsi$ID)]; pl$Time <- ipsi$time[match(pl$ID, ipsi$ID)]
key <- c("IFN/STING score", sting5, "Tap1", "H2-K1", "Psmb8"); pl$feature <- factor(pl$feature, levels = key)
pv <- gt[match(c("STING5_zmean", key[-1]), gt$feature), ]
pl_lab <- data.frame(feature = factor(key, levels = key), lab = sprintf("p = %s", fmtp(pv$wilcox_p)))
f5a <- ggplot(pl, aes(Group, value)) + geom_boxplot(aes(fill = Group), outlier.shape = NA, width = 0.55, linewidth = 0.3, alpha = 0.6) +
  geom_point(aes(shape = Time), position = position_jitter(width = 0.08, seed = 1), size = 1.2) +
  geom_text(data = pl_lab, aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, size = 2.2, inherit.aes = FALSE) +
  facet_wrap(~feature, scales = "free_y", nrow = 2) + scale_fill_manual(values = COL_GRP) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  th + theme(strip.text = element_text(face = "italic")) + labs(x = NULL, y = "log2 normalised count (score: z-mean)", shape = "Time",
  title = "GSE131193 brain microvessels, ipsilateral hemisphere: tMCAO vs sham (n = 6 vs 6; exact Wilcoxon)")
ct <- read.csv(rp("06_ssgsea_sting", "ssGSEA_group_tests.csv")); ct <- ct[ct$comparison == "primary_ipsi_tMCAO_vs_ipsi_sham", ]
ct <- ct[order(ct$diff), ]; ct$cell <- factor(ct$cell, levels = ct$cell); ct$sig <- ifelse(ct$FDR_wilcox < 0.05, "FDR < 0.05", "FDR \u2265 0.05")
f5b <- ggplot(ct, aes(diff, cell, fill = sig)) + geom_col(width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("FDR < 0.05" = "#E45756", "FDR \u2265 0.05" = "grey70")) + th + theme(axis.text.y = element_text(size = 6)) +
  labs(x = "ssGSEA score difference (tMCAO \u2212 sham)", y = NULL, fill = NULL, title = sprintf("28 immune signatures (Charoentong 2017):\n%d with FDR < 0.05", sum(ct$FDR_wilcox < 0.05)))
ss <- read_mat(rp("06_ssgsea_sting", "ssGSEA_scores.csv"))
cr <- read.csv(rp("06_ssgsea_sting", "STING_score_vs_ssGSEA_spearman.csv"))
r12 <- cr[cr$samples == "primary12" & cr$score == "STING5" & cr$cell == "Activated CD8 T cell", ]
d5c <- data.frame(score = sc, cd8 = ss[ipsi$ID, "Activated CD8 T cell"], Group = ipsi$Group)
f5c <- ggplot(d5c, aes(score, cd8, colour = Group)) + geom_point(size = 1.8) + scale_colour_manual(values = COL_GRP) + th +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 2.3, label = sprintf("Spearman rho = %.3f\nFDR = %s (n = 12)", r12$rho, fmtp(r12$FDR))) +
  labs(x = "IFN/STING score (z-mean of 5 genes)", y = "Activated CD8 T cell (ssGSEA)", title = "IFN/STING score vs activated CD8 T cell")
vv <- cr[cr$cell == "Activated CD8 T cell", ]
vv$score <- factor(vv$score, levels = c("STING5", "STING4_noIsg15", "STING4_noCxcl10", "STING3_noISG"),
                   labels = c("5 genes", "without Isg15", "without Cxcl10", "without Isg15 & Cxcl10"))
vv$samples <- factor(vv$samples, levels = c("primary12", "all24"), labels = c("ipsilateral (n = 12)", "all samples (n = 24)"))
f5d <- ggplot(vv, aes(score, rho, fill = samples)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", rho)), position = position_dodge(0.75), vjust = -0.3, size = 2) +
  scale_fill_manual(values = c("#9ecae1", "#3182bd")) + scale_y_continuous(limits = c(0, 1.08)) + th + theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(x = "IFN/STING score version", y = "Spearman rho", fill = NULL, title = "IFN/STING score vs activated CD8 T cell:\nscore versions without ISGs") + theme(legend.position = "bottom")
save_fig(f5a / (f5b | (f5c / f5d)) + plot_layout(heights = c(1, 1.45)) + plot_annotation(tag_levels = "A"), "Figure_5", 7.2, 8.6)

## ---------------- Figure 6: enrichment and PPI ----------------
kg <- read.csv(rp("04_enrichment_ppi", "KEGG.csv")); bp <- read.csv(rp("04_enrichment_ppi", "GO_BP.csv"))
cc <- read.csv(rp("04_enrichment_ppi", "GO_CC.csv")); mf <- read.csv(rp("04_enrichment_ppi", "GO_MF.csv"))
top <- function(d, src, n = 8) { d <- head(d[order(d$p.adjust), ], n); data.frame(src = src, term = d$Description, Count = d$Count, FDR = d$p.adjust) }
kg_sel <- kg[kg$ID %in% c("mmu04612", "mmu04623", "mmu04621", "mmu04620"), ]
en <- rbind(top(bp, "GO BP"), top(cc, "GO CC", 5), top(mf, "GO MF", 5),
            unique(rbind(top(kg, "KEGG", 8), data.frame(src = "KEGG", term = kg_sel$Description, Count = kg_sel$Count, FDR = kg_sel$p.adjust))))
en <- en[!duplicated(paste(en$src, en$term)), ]
en$term <- factor(paste0(en$term, "  "), levels = rev(unique(paste0(en$term, "  "))))
en$hl <- grepl("Antigen processing|Cytosolic DNA|MHC class I|interferon", en$term, ignore.case = TRUE)
f6a <- ggplot(en, aes(-log10(FDR), term)) + geom_col(aes(fill = src), width = 0.7) + geom_text(aes(label = Count), hjust = -0.2, size = 2) +
  facet_grid(src ~ ., scales = "free_y", space = "free_y") + scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
  scale_fill_brewer(palette = "Set2", guide = "none") + th + theme(axis.text.y = element_text(size = 6)) +
  labs(x = "-log10 FDR (number = DEGs in term)", y = NULL, title = sprintf("GO / KEGG enrichment of\nthe %d DEGs (top terms)", nrow(read.csv(rp("03_combat_limma", "TableS6_significant_DEGs.csv")))))
## (B) PPI, concentric layout grouped by gene class. Nodes shown: the 40 highest-degree nodes of the full STRING network
## plus the selected cGAS-STING/IFN and MHC-I genes. Each class occupies its own ring (inner: cGAS-STING/IFN; middle:
## MHC-I antigen processing; outer: other high-degree nodes); within a ring nodes are ordered by degree.
## Edges: all STRING edges (score >= 0.7) among the shown nodes. Degree and rank are those of the full network.
nodes <- read.csv(rp("04_enrichment_ppi", "PPI_nodes_score700.csv")); edges <- read.csv(rp("04_enrichment_ppi", "PPI_edges_score700.csv"))
ifn_g <- c("Stat1", "Cxcl10", "Isg15", "Irf7", "Tmem173"); mhc_g <- c("Tap1", "Tap2", "Tapbp", "B2m", "H2-K1", "H2-D1", "Psmb8")
sub <- unique(c(head(nodes$Gene, 40), intersect(c(ifn_g, mhc_g), nodes$Gene)))
vd <- data.frame(name = sub, Degree = nodes$Degree[match(sub, nodes$Gene)], rank = match(sub, nodes$Gene))
vd$class <- factor(ifelse(vd$name %in% mhc_g, "MHC-I antigen processing", ifelse(vd$name %in% ifn_g, "cGAS-STING / IFN", "other high-degree node")),
                   levels = c("cGAS-STING / IFN", "MHC-I antigen processing", "other high-degree node"))
vd$ring <- as.integer(vd$class)
Rr <- c(1.0, 2.35, 4.1)
vd <- vd[order(vd$ring, -vd$Degree), ]
off <- c(pi / 2, pi / 2 + pi / 7, pi / 2)
vd$theta <- NA_real_
for (k in 1:3) { i <- which(vd$ring == k); n <- length(i); vd$theta[i] <- off[k] - (seq_len(n) - 1) / n * 2 * pi }
vd$x <- Rr[vd$ring] * cos(vd$theta); vd$y <- Rr[vd$ring] * sin(vd$theta)
es <- edges[edges[, 1] %in% sub & edges[, 2] %in% sub, ]
ed <- data.frame(x = vd$x[match(es[, 1], vd$name)], y = vd$y[match(es[, 1], vd$name)], xend = vd$x[match(es[, 2], vd$name)], yend = vd$y[match(es[, 2], vd$name)])
ca <- vd$class[match(es[, 1], vd$name)]; cb <- vd$class[match(es[, 2], vd$name)]
ed$type <- factor(ifelse(ca == "MHC-I antigen processing" | cb == "MHC-I antigen processing", "edge with an MHC-I gene",
                         ifelse(ca == "cGAS-STING / IFN" | cb == "cGAS-STING / IFN", "edge with a cGAS-STING / IFN gene", "other edge")),
                  levels = c("edge with an MHC-I gene", "edge with a cGAS-STING / IFN gene", "other edge"))
ed <- ed[order(ed$type, decreasing = TRUE), ]
## labels: radial outside the node on the outer ring; inner rings get horizontal labels offset outward
vd$lr <- Rr[vd$ring] + 0.2 + 0.03 * sqrt(vd$Degree)
vd$lx <- vd$lr * cos(vd$theta); vd$ly <- vd$lr * sin(vd$theta)
deg <- vd$theta * 180 / pi; left <- cos(vd$theta) < -1e-9
vd$angle <- ifelse(vd$ring == 3, ifelse(left, deg + 180, deg), 0)
vd$hjust <- ifelse(vd$ring == 3, ifelse(left, 1, 0), 0.5 - 0.5 * cos(vd$theta))
vd$vjust <- ifelse(vd$ring == 3, 0.5, 0.5 - 0.5 * sin(vd$theta))
rings <- do.call(rbind, lapply(1:3, function(k) data.frame(ring = k, t = seq(0, 2 * pi, length.out = 361), r = Rr[k])))
nc <- table(vd$class)
f6b <- ggplot() +
  geom_path(data = rings, aes(r * cos(t), r * sin(t), group = ring), colour = "grey88", linetype = 2, linewidth = 0.25) +
  geom_segment(data = ed, aes(x, y, xend = xend, yend = yend, colour = type, alpha = type), linewidth = 0.22) +
  scale_colour_manual(values = c("edge with an MHC-I gene" = "#E45756", "edge with a cGAS-STING / IFN gene" = "#E0A100", "other edge" = "grey60"), name = NULL, drop = FALSE) +
  scale_alpha_manual(values = c("edge with an MHC-I gene" = 0.8, "edge with a cGAS-STING / IFN gene" = 0.55, "other edge" = 0.16), guide = "none") +
  geom_point(data = vd, aes(x, y, size = Degree, fill = class), shape = 21, stroke = 0.25, colour = "grey20") +
  geom_text(data = vd, aes(lx, ly, label = name, angle = angle, hjust = hjust, vjust = vjust), size = 1.9, fontface = "italic") +
  scale_fill_manual(values = c("cGAS-STING / IFN" = "#F2B701", "MHC-I antigen processing" = "#E45756", "other high-degree node" = "#9ecae1"),
                    labels = c(sprintf("cGAS-STING / IFN (inner ring, n = %d)", nc[1]), sprintf("MHC-I antigen processing (middle ring, n = %d)", nc[2]),
                               sprintf("other top-40 node (outer ring, n = %d)", nc[3])), name = NULL) +
  scale_size_continuous(range = c(1.3, 5.5), name = "Degree\n(full network)") +
  coord_equal(xlim = c(-5.3, 5.3), ylim = c(-5.0, 5.0), clip = "off") + theme_void(base_size = 8) +
  theme(plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5), plot.title.position = "plot",
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 6.5)) +
  guides(fill = guide_legend(override.aes = list(size = 3), order = 1), colour = guide_legend(order = 2), size = guide_legend(order = 3)) +
  labs(title = sprintf("STRING v11.5 PPI of the DEGs (score \u2265 0.7): full network %d nodes, %d edges\nShown: top-40 nodes by degree plus selected genes (%d nodes, %d edges); one ring per gene class",
                       nrow(nodes), nrow(edges), nrow(vd), nrow(es)))
write.csv(vd[, c("name", "class", "Degree", "rank", "ring")], file.path(FIG, "Figure_6C_nodes_shown.csv"), row.names = FALSE)
## (C) mmu04612: every DEG of the pathway (KEGG REST gene list) with its log2FC and FDR
pg <- strsplit(kg$geneID[kg$ID == "mmu04612"], "/")[[1]]
mapE <- read.csv(rp("04_enrichment_ppi", "DEG_to_Entrez_map.csv")); inp <- mapE$input[match(pg, mapE$official)]
d6c <- data.frame(gene = pg, logFC = A[inp, "logFC"], FDR = A[inp, "adj.P.Val"])
d6c <- d6c[order(d6c$logFC), ]; d6c$gene <- factor(d6c$gene, levels = d6c$gene)
f6c <- ggplot(d6c, aes(logFC, gene)) + geom_col(fill = "#E45756", width = 0.7) + geom_text(aes(label = fmtp(FDR)), hjust = -0.15, size = 1.9) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.62))) + coord_cartesian(clip = "off") + th + theme(axis.text.y = element_text(face = "italic", size = 6), plot.margin = margin(5.5, 16, 5.5, 5.5)) +
  labs(x = "log2FC (label: FDR)", y = NULL, title = sprintf("KEGG mmu04612 (antigen processing\nand presentation): %d DEGs, FDR %s", nrow(d6c), fmtp(kg$p.adjust[kg$ID == "mmu04612"])))
save_fig(((f6a | f6c) + plot_layout(widths = c(1.35, 1))) / f6b + plot_layout(heights = c(1.15, 1)) + plot_annotation(tag_levels = "A"), "Figure_6", 7.2, 9.6)
