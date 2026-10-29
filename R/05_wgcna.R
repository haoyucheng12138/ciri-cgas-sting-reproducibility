# =============================================================================
# 05_wgcna.R  -  WGCNA on GSE131193 (brain microvessel fragments, 24 samples).
#  Samples: all 24 (no sample has standardized connectivity Z.k < -2.5).
#  Network: unsigned, Pearson; soft-thresholding power = lowest power with signed scale-free R^2 >= 0.85
#  (powers 1-10, 12-30); blockwiseModules(minModuleSize = 30, deepSplit = 2, mergeCutHeight = 0.3,
#  pamRespectsDendro = FALSE, randomSeed = 12345, single block). Single-threaded for reproducibility on Windows.
#  Traits (from GEO sample titles): primary = ipsilateral tMCAO vs ipsilateral sham (6 vs 6; other samples NA).
#  Secondary traits: tMCAO (all 24), day 7, ipsilateral, ipsilateral vs contralateral within tMCAO.
# =============================================================================
source("R/00_config.R")
suppressMessages(library(WGCNA)); disableWGCNAThreads()
out <- res_dir("05_wgcna")
expr <- read_mat(file.path(DIR_PRO, "GSE131193_log2counts.csv"))
ann  <- read.csv(file.path(DIR_PRO, "GSE131193_sample_annotation.csv"))
dat  <- t(expr[, ann$ID])
gsg <- goodSamplesGenes(dat, verbose = 0); stopifnot(gsg$allOK)

Zk <- as.numeric(scale(colSums(adjacency(t(dat), type = "distance")))); names(Zk) <- rownames(dat)
write.csv(data.frame(sample = names(Zk), Zk = Zk), file.path(out, "sample_connectivity_Zk.csv"), row.names = FALSE)
cat("samples with Z.k < -2.5:", sum(Zk < -2.5), "| min Z.k:", round(min(Zk), 2), "\n")
stopifnot(sum(Zk < -2.5) == 0)

powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(dat, powerVector = powers, networkType = "unsigned", RsquaredCut = 0.85, verbose = 0)
fi <- sft$fitIndices; fi$signedR2 <- -sign(fi$slope) * fi$SFT.R.sq
write.csv(fi, file.path(out, "soft_threshold.csv"), row.names = FALSE)
power <- fi$Power[which(fi$signedR2 >= 0.85)[1]]
cat(sprintf("soft power = %d (signed R^2 = %.3f, mean k = %.1f)\n", power, fi$signedR2[fi$Power == power], fi$mean.k.[fi$Power == power]))

net <- blockwiseModules(dat, power = power, networkType = "unsigned", TOMType = "unsigned", minModuleSize = 30, deepSplit = 2,
                        pamRespectsDendro = FALSE, mergeCutHeight = 0.3, numericLabels = TRUE, saveTOMs = FALSE,
                        maxBlockSize = ncol(dat), randomSeed = 12345, verbose = 0)
col <- labels2colors(net$colors); names(col) <- colnames(dat)
MEs <- orderMEs(moduleEigengenes(dat, col)$eigengenes)
saveRDS(list(colors = col, MEs = MEs, power = power, net_colors = net$colors, dendro = net$dendrograms), file.path(out, "wgcna_network.rds"))

a <- ann[match(rownames(dat), ann$ID), ]
traits <- data.frame(row.names = rownames(dat),
  ipsi_tMCAO_vs_ipsi_sham  = ifelse(a$side == "ipsilateral", as.numeric(a$injury == "tMCAO"), NA),
  tMCAO_all                = as.numeric(a$injury == "tMCAO"),
  day7                     = as.numeric(a$time == "7d"),
  ipsilateral              = as.numeric(a$side == "ipsilateral"),
  ipsi_vs_contra_in_tMCAO  = ifelse(a$injury == "tMCAO", as.numeric(a$side == "ipsilateral"), NA),
  ipsi_tMCAO_vs_ipsi_sham_1d = ifelse(a$side == "ipsilateral" & a$time == "1d", as.numeric(a$injury == "tMCAO"), NA))
write.csv(traits, file.path(out, "traits.csv"))
R <- sapply(traits, function(t) { ok <- !is.na(t); sapply(MEs, function(me) cor(me[ok], t[ok])) })
P <- sapply(traits, function(t) { ok <- !is.na(t); sapply(MEs, function(me) cor.test(me[ok], t[ok])$p.value) })
size <- as.integer(table(col)[sub("^ME", "", rownames(R))])
mt <- data.frame(module = rownames(R), size = size, R, check.names = FALSE); mp <- data.frame(module = rownames(P), P, check.names = FALSE)
write.csv(mt, file.path(out, "module_trait_r.csv"), row.names = FALSE); write.csv(mp, file.path(out, "module_trait_p.csv"), row.names = FALSE)
prim <- data.frame(module = rownames(R), size = size, r = round(R[, 1], 4), p = signif(P[, 1], 3), FDR = signif(p.adjust(P[, 1], "BH"), 3))
prim <- prim[order(prim$p), ]; write.csv(prim, file.path(out, "module_trait_primary_sorted.csv"), row.names = FALSE)
cat("\nPrimary trait (ipsilateral tMCAO vs ipsilateral sham, n = 12):\n"); print(prim, row.names = FALSE)

kME <- signedKME(dat, MEs); write.csv(kME, file.path(out, "kME_all_genes.csv"))
write.csv(data.frame(gene = names(col), module = col), file.path(out, "gene_module_assignment.csv"), row.names = FALSE)
focal <- c("Mb21d1","Tmem173","Irf3","Isg15","Cxcl10","Stat1","Il6","Trpc6","Tap1","Tap2","Tapbp","Tapbpl","B2m","H2-K1","H2-D1",
           "Psmb8","Psmb9","Psmb10","Irgm1","Irgm2","Igtp","Trim56","Irf7","Ifit3")
focal <- intersect(focal, names(col))
fg <- data.frame(gene = focal, module = col[focal], kME = round(sapply(focal, function(g) kME[g, paste0("kME", col[g])]), 3))
fg <- merge(fg, prim[, c("module", "r", "p")], by.x = "module", by.y = "module", all.x = TRUE, sort = FALSE)
fg$module <- sub("^ME", "", fg$module)
fg <- data.frame(gene = focal, module = col[focal], kME = round(sapply(focal, function(g) kME[g, paste0("kME", col[g])]), 3),
                 module_r_primary = round(R[paste0("ME", col[focal]), 1], 3), module_p_primary = signif(P[paste0("ME", col[focal]), 1], 3))
write.csv(fg, file.path(out, "focal_genes.csv"), row.names = FALSE)
cat("\nFocal genes:\n"); print(fg, row.names = FALSE)
save_session("05_wgcna")
