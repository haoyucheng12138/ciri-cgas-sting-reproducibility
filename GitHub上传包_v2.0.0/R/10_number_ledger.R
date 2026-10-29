# =============================================================================
# 10_number_ledger.R  -  collect every number used in the revised manuscript from the result files.
# Output: results/number_ledger.csv (item, value, source file). Nothing in this table is typed by hand.
# =============================================================================
source("R/00_config.R")
L <- list(); add <- function(section, item, value, file) L[[length(L) + 1]] <<- data.frame(section = section, item = item, value = as.character(value), source = file)
f <- function(x, d = 3) formatC(x, digits = d, format = "g")
rp <- function(step, file) file.path(DIR_RES, step, file)

## 02 preprocessing
m <- read_mat(file.path(DIR_PRO, "merged_log2FPKM.csv")); g1 <- read_mat(file.path(DIR_PRO, "GSE131193_log2counts.csv")); g2 <- read_mat(file.path(DIR_PRO, "GSE23163_log2ratio.csv"))
add("Preprocessing", "merged cohort genes x samples", sprintf("%d x %d", nrow(m), ncol(m)), "data/processed/merged_log2FPKM.csv")
add("Preprocessing", "GSE131193 genes x samples", sprintf("%d x %d", nrow(g1), ncol(g1)), "data/processed/GSE131193_log2counts.csv")
add("Preprocessing", "GSE23163 genes x samples", sprintf("%d x %d", nrow(g2), ncol(g2)), "data/processed/GSE23163_log2ratio.csv")
cm <- read.csv(file.path(DIR_PRO, "sample_column_map.csv")); a <- cm[cm$dataset == "GSE131193", ]
add("Preprocessing", "GSE131193 file columns not in accession order", sum(a$GSM != sort(a$GSM)), "data/processed/sample_column_map.csv")

## 03 merged cohort
s <- read.csv(rp("03_combat_limma", "R2.5_batch_sensitivity_summary.csv"))
for (i in seq_len(nrow(s))) add("Merged cohort DEG", s$strategy[i], sprintf("%d DEG (%d up / %d down); shared with primary %d", s$n_DEG[i], s$up[i], s$down[i], s$shared_with_primary[i]), "results/03_combat_limma/R2.5_batch_sensitivity_summary.csv")
add("Merged cohort DEG", "direction concordance", s$note[1], "results/03_combat_limma/R2.5_batch_sensitivity_summary.csv")
k <- read.csv(rp("03_combat_limma", "R2.5_key_genes_across_strategies.csv"))
for (i in seq_len(nrow(k))) add("Merged cohort key genes", k$gene[i], sprintf("log2FC %.3f; FDR A %s, B %s, C %s; P within GSE131712 %s, GSE202659 %s", k$A_logFC[i], f(k$A_FDR[i]), f(k$B_FDR[i]), f(k$C_FDR[i]), f(k$GSE131712_P[i]), f(k$GSE202659_P[i])), "results/03_combat_limma/R2.5_key_genes_across_strategies.csv")
st <- read.csv(rp("03_combat_limma", "STING_score_merged_cohort.csv"))
add("Merged cohort STING score", st$score, sprintf("diff %.3f; Welch p %s; exact Wilcoxon p %s", st$diff, f(st$welch_p), f(st$wilcoxon_exact_p)), "results/03_combat_limma/STING_score_merged_cohort.csv")

## 04 enrichment / PPI
kg <- read.csv(rp("04_enrichment_ppi", "KEGG.csv"))
add("Enrichment", "KEGG significant pathways", nrow(kg), "results/04_enrichment_ppi/KEGG.csv")
for (id in c("mmu04621", "mmu04612", "mmu04623", "mmu04620")) { r <- kg[kg$ID == id, ]; if (nrow(r)) add("Enrichment", paste(id, r$Description), sprintf("%d genes; FDR %s", r$Count, f(r$p.adjust)), "results/04_enrichment_ppi/KEGG.csv") }
for (o in c("BP", "CC", "MF")) { x <- read.csv(rp("04_enrichment_ppi", paste0("GO_", o, ".csv"))); add("Enrichment", paste("GO", o, "terms; top 3"), sprintf("%d; %s", nrow(x), paste(head(x$Description, 3), collapse = " | ")), paste0("results/04_enrichment_ppi/GO_", o, ".csv")) }
pn <- read.csv(rp("04_enrichment_ppi", "PPI_nodes_score700.csv")); pe <- read.csv(rp("04_enrichment_ppi", "PPI_edges_score700.csv"))
add("PPI", "nodes / edges (STRING 11.5, >= 0.7)", sprintf("%d / %d", nrow(pn), nrow(pe)), "results/04_enrichment_ppi/PPI_nodes_score700.csv")
add("PPI", "top 20 hubs (degree)", paste(sprintf("%s(%d)", head(pn$Gene, 20), head(pn$Degree, 20)), collapse = " "), "results/04_enrichment_ppi/PPI_nodes_score700.csv")
ap <- c("Tap1", "Tap2", "Tapbp", "B2m", "H2-K1", "H2-D1", "Psmb8", "Stat1", "Cxcl10", "Isg15", "Tmem173", "Irf7")
add("PPI", "degree rank of selected genes", paste(sprintf("%s=%s", ap, match(ap, pn$Gene)), collapse = " "), "results/04_enrichment_ppi/PPI_nodes_score700.csv")

## 05 WGCNA
sf <- read.csv(rp("05_wgcna", "soft_threshold.csv")); zk <- read.csv(rp("05_wgcna", "sample_connectivity_Zk.csv"))
pw <- sf$Power[which(sf$signedR2 >= 0.85)[1]]
add("WGCNA", "samples / outliers (Z.k < -2.5) / min Z.k", sprintf("%d / %d / %.2f", nrow(zk), sum(zk$Zk < -2.5), min(zk$Zk)), "results/05_wgcna/sample_connectivity_Zk.csv")
add("WGCNA", "soft power / signed R2 / mean k", sprintf("%d / %.3f / %.1f", pw, sf$signedR2[sf$Power == pw], sf$mean.k.[sf$Power == pw]), "results/05_wgcna/soft_threshold.csv")
pr <- read.csv(rp("05_wgcna", "module_trait_primary_sorted.csv"))
add("WGCNA", "modules (excluding grey)", sum(pr$module != "MEgrey"), "results/05_wgcna/module_trait_primary_sorted.csv")
for (i in which(pr$p < 0.05)) add("WGCNA primary trait", pr$module[i], sprintf("size %d; r %.3f; p %s; FDR %s", pr$size[i], pr$r[i], f(pr$p[i]), f(pr$FDR[i])), "results/05_wgcna/module_trait_primary_sorted.csv")
fg <- read.csv(rp("05_wgcna", "focal_genes.csv"))
for (i in seq_len(nrow(fg))) add("WGCNA focal genes", fg$gene[i], sprintf("%s; kME %.3f; module r %.3f, p %s", fg$module[i], fg$kME[i], fg$module_r_primary[i], f(fg$module_p_primary[i])), "results/05_wgcna/focal_genes.csv")

## 06 ssGSEA / STING score
gt <- read.csv(rp("06_ssgsea_sting", "genes_and_STING_score_group_tests.csv"))
for (i in seq_len(nrow(gt))) add(paste("GSE131193", gt$comparison[i]), gt$feature[i], sprintf("n %d vs %d; diff %.3f; exact Wilcoxon p %s; Welch p %s", gt$n_case[i], gt$n_ctrl[i], gt$diff[i], f(gt$wilcox_p[i]), f(gt$welch_p[i])), "results/06_ssgsea_sting/genes_and_STING_score_group_tests.csv")
ct <- read.csv(rp("06_ssgsea_sting", "ssGSEA_group_tests.csv"))
for (cmp in unique(ct$comparison)) { x <- ct[ct$comparison == cmp, ]; x <- x[order(x$wilcox_p), ]
  add("GSE131193 ssGSEA", paste(cmp, "cell types FDR < 0.05"), sum(x$FDR_wilcox < 0.05), "results/06_ssgsea_sting/ssGSEA_group_tests.csv")
  for (i in which(x$FDR_wilcox < 0.05)) add("GSE131193 ssGSEA", paste(cmp, x$cell[i]), sprintf("diff %.4f; p %s; FDR %s", x$diff[i], f(x$wilcox_p[i]), f(x$FDR_wilcox[i])), "results/06_ssgsea_sting/ssGSEA_group_tests.csv") }
cr <- read.csv(rp("06_ssgsea_sting", "STING_score_vs_ssGSEA_spearman.csv"))
for (smp in c("all24", "primary12")) for (v in unique(cr$score)) { x <- cr[cr$samples == smp & cr$score == v, ]; x <- x[order(x$p), ]
  add("STING score vs ssGSEA (Spearman)", paste(smp, v, "top 3"), paste(sprintf("%s rho %.3f FDR %s", x$cell[1:3], x$rho[1:3], f(x$FDR[1:3])), collapse = "; "), "results/06_ssgsea_sting/STING_score_vs_ssGSEA_spearman.csv")
  y <- x[x$cell == "Activated CD8 T cell", ]; add("STING score vs ssGSEA (Spearman)", paste(smp, v, "Activated CD8 T cell"), sprintf("rho %.3f; FDR %s", y$rho, f(y$FDR)), "results/06_ssgsea_sting/STING_score_vs_ssGSEA_spearman.csv") }

## 07 intersection
ic <- read.csv(rp("07_intersection", "intersection_counts.csv"))
for (i in seq_len(nrow(ic))) add("Intersection", ic$pair[i], ic$n[i], "results/07_intersection/intersection_counts.csv")
for (nm in c("DEG", "WGCNA_core", "cGAS_STING")) add("Intersection", paste("set size", nm), length(readLines(rp("07_intersection", paste0("set_", nm, ".txt")))), paste0("results/07_intersection/set_", nm, ".txt"))
ev <- read.csv(rp("07_intersection", "intersection_evidence.csv"))
add("Intersection", "three-set genes", paste(ev$gene[ev$in_all_three], collapse = ", "), "results/07_intersection/intersection_evidence.csv")
for (i in seq_len(nrow(ev))) add("Intersection evidence", ev$gene[i], sprintf("log2FC %.3f; FDR %s; module %s (assoc %s); kME %.3f; all three %s", ev$logFC[i], f(ev$FDR[i]), ev$module[i], ev$module_assoc_primary[i], ev$kME_own[i], ev$in_all_three[i]), "results/07_intersection/intersection_evidence.csv")

## 08 Elastic Net
en <- read.csv(rp("08_elastic_net", "elastic_net_summary.csv"))
for (i in seq_len(nrow(en))) add("Elastic Net", paste(en$feature_space[i], en$task[i]),
  sprintf("candidates %d; non-zero %d; lambda.min %s; test AUC %.3f (95%% CI %.3f-%.3f); accuracy %.3f; top10 %s; antigen ranks %s", en$n_candidates[i], en$n_nonzero[i], f(en$lambda_min[i]),
          en$test_AUC[i], en$CI_low[i], en$CI_high[i], en$accuracy[i], en$top10_by_abs_coef[i], ifelse(nzchar(en$antigen_genes_rank[i]), en$antigen_genes_rank[i], "none")), "results/08_elastic_net/elastic_net_summary.csv")

## 09 rat
rr <- read.csv(rp("09_rat_roc", "rat_ROC.csv"))
for (i in seq_len(nrow(rr))) add("Rat GSE227317 (5 vs 5)", rr$feature[i], sprintf("AUC %.3f (95%% CI %.2f-%.2f); mean log2 MCAO %.3f vs Sham %.3f; exact Wilcoxon p %s", rr$AUC[i], rr$CI_low[i], rr$CI_high[i], rr$mean_MCAO[i], rr$mean_Sham[i], f(rr$wilcox_exact_p[i])), "results/09_rat_roc/rat_ROC.csv")

led <- do.call(rbind, L); write.csv(led, file.path(DIR_RES, "number_ledger.csv"), row.names = FALSE)
cat("ledger rows:", nrow(led), "\n")
