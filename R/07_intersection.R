# =============================================================================
# 07_intersection.R  -  three-set intersection.
#  Set 1: 1,421 merged-cohort DEGs (03).
#  Set 2: WGCNA core genes = genes in modules whose eigengene correlates with the primary trait
#         (ipsilateral tMCAO vs ipsilateral sham; |r| >= 0.4 and p < 0.05; grey excluded) and |kME| >= 0.5 (05).
#  Set 3: curated cGAS-STING pathway genes (reference/cGAS_STING_curated_genes.csv; 29 valid symbols).
# =============================================================================
source("R/00_config.R")
out <- res_dir("07_intersection")
deg <- read.csv(file.path(res_dir("03_combat_limma"), "TableS6_significant_DEGs.csv"), row.names = 1)
mt  <- read.csv(file.path(res_dir("05_wgcna"), "module_trait_r.csv"), check.names = FALSE)
mp  <- read.csv(file.path(res_dir("05_wgcna"), "module_trait_p.csv"), check.names = FALSE)
kME <- read_mat(file.path(res_dir("05_wgcna"), "kME_all_genes.csv"))
asg <- read.csv(file.path(res_dir("05_wgcna"), "gene_module_assignment.csv"))
tr  <- "ipsi_tMCAO_vs_ipsi_sham"
sig <- sub("^ME", "", mt$module[abs(mt[[tr]]) >= 0.4 & mp[[tr]] < 0.05 & mt$module != "MEgrey"])
cat("modules associated with the primary trait:", paste(sig, collapse = ", "), "\n")
core <- unlist(lapply(sig, function(m) { g <- asg$gene[asg$module == m]; g[abs(kME[g, paste0("kME", m)]) >= 0.5] }))
cg <- read.csv(file.path(DIR_REF, "cGAS_STING_curated_genes.csv"))$Gene
cg <- intersect(cg, asg$gene)
sets <- list(DEG = rownames(deg), WGCNA_core = core, cGAS_STING = cg)
for (nm in names(sets)) writeLines(sort(sets[[nm]]), file.path(out, paste0("set_", nm, ".txt")))
pair <- data.frame(pair = c("DEG & WGCNA_core", "DEG & cGAS_STING", "WGCNA_core & cGAS_STING", "all three"),
                   n = c(length(intersect(sets$DEG, core)), length(intersect(sets$DEG, cg)), length(intersect(core, cg)), length(Reduce(intersect, sets))))
cat(sprintf("set sizes: DEG %d | WGCNA core %d | cGAS-STING %d\n", length(sets$DEG), length(core), length(cg))); print(pair, row.names = FALSE)
write.csv(pair, file.path(out, "intersection_counts.csv"), row.names = FALSE)
tri <- Reduce(intersect, sets); cat("three-set intersection:", paste(tri, collapse = ", "), "\n")
ev <- data.frame(gene = sort(union(intersect(sets$DEG, cg), tri)))
ev$logFC <- round(deg[ev$gene, "logFC"], 3); ev$FDR <- signif(deg[ev$gene, "adj.P.Val"], 3)
ev$module <- asg$module[match(ev$gene, asg$gene)]
ev$kME_own <- round(sapply(seq_len(nrow(ev)), function(i) kME[ev$gene[i], paste0("kME", ev$module[i])]), 3)
ev$module_assoc_primary <- ev$module %in% sig; ev$in_WGCNA_core <- ev$gene %in% core; ev$in_all_three <- ev$gene %in% tri
print(ev, row.names = FALSE); write.csv(ev, file.path(out, "intersection_evidence.csv"), row.names = FALSE)
save_session("07_intersection")
