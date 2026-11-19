# =============================================================================
# 04_enrichment_check.R  -  independent verification of the GO / KEGG results of 04_enrichment_ppi.R.
#  (1) every reported p-value is recomputed from the table's own counts with the hypergeometric test
#      (phyper(k-1, M, N-M, n, lower.tail = FALSE)) and every FDR with p.adjust(method = "BH");
#  (2) every gene listed for a term must be one of the DEGs (official symbol after alias mapping);
#  (3) for key KEGG pathways the gene list is re-downloaded from KEGG REST (link/mmu/<pathway>) and the reported
#      genes must equal DEG genes in that official list.
# Writes results/04_enrichment_ppi/enrichment_verification.csv. Stops with an error if any check fails.
# =============================================================================
source("R/00_config.R")
suppressMessages({ library(org.Mm.eg.db); library(httr) })
out <- res_dir("04_enrichment_ppi")
map <- read.csv(file.path(out, "DEG_to_Entrez_map.csv"))
deg_sym <- unique(map$official); deg_ent <- unique(as.character(map$ENTREZID))
frac <- function(x) { z <- do.call(rbind, strsplit(x, "/")); matrix(as.numeric(z), ncol = 2) }
res <- list(); fail <- 0
for (tab in c("KEGG", "GO_BP", "GO_CC", "GO_MF")) {
  d <- read.csv(file.path(out, paste0(tab, ".csv")))
  gr <- frac(d$GeneRatio); br <- frac(d$BgRatio)
  k <- gr[, 1]; n <- gr[, 2]; M <- br[, 1]; N <- br[, 2]
  p_re <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
  cnt_ok <- all(d$Count == k)
  genes_ok <- vapply(strsplit(d$geneID, "/"), function(g) all(g %in% deg_sym), TRUE)
  pdiff <- max(abs(log10(p_re) - log10(d$pvalue)))
  ## BH FDR recomputed over ALL tested terms (saved by 04_enrichment_ppi.R), then compared term by term.
  allt <- read.csv(file.path(out, paste0(tab, "_all_tested.csv")))
  fdr_re <- p.adjust(allt$pvalue, method = "BH"); names(fdr_re) <- allt$ID
  fdr_diff <- max(abs(fdr_re[d$ID] - d$p.adjust))
  same_set <- setequal(d$ID, allt$ID[fdr_re < 0.05 & allt$qvalue < 0.2])
  ok <- cnt_ok && all(genes_ok) && pdiff < 1e-6 && fdr_diff < 1e-12 && same_set
  if (!ok) fail <- fail + 1
  res[[tab]] <- data.frame(table = tab, n_tested = nrow(allt), n_significant = nrow(d), Count_equals_k = cnt_ok, all_genes_are_DEGs = all(genes_ok),
                           max_abs_log10_p_difference = signif(pdiff, 3), max_abs_FDR_difference = signif(fdr_diff, 3),
                           significant_set_reproduced = same_set, n_DEGs_annotated = paste(unique(n), collapse = ";"),
                           background_N = paste(unique(N), collapse = ";"), PASS = ok)
}
v <- do.call(rbind, res); print(v, row.names = FALSE)

## (3) KEGG REST cross-check of key pathways
kd <- read.csv(file.path(out, "KEGG.csv"))
ky <- list()
for (pid in c("mmu04612", "mmu04623", "mmu04621", "mmu04620")) {
  r <- GET(paste0("https://rest.kegg.jp/link/mmu/", pid), timeout(120)); stop_for_status(r)
  lk <- read.delim(text = content(r, "text", encoding = "UTF-8"), header = FALSE)
  ent <- sub("^mmu:", "", lk$V2)
  writeLines(ent, file.path(out, paste0("KEGG_REST_", pid, "_entrez.txt")))
  expected <- sort(unique(suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, intersect(ent, deg_ent), "SYMBOL", "ENTREZID"))))
  reported <- sort(strsplit(kd$geneID[kd$ID == pid], "/")[[1]])
  same <- identical(expected, reported)
  if (!same) fail <- fail + 1
  ky[[pid]] <- data.frame(pathway = pid, KEGG_pathway_genes = length(unique(ent)), DEGs_in_pathway_by_KEGG_REST = length(expected),
                          reported_Count = length(reported), identical_gene_sets = same,
                          only_in_KEGG_REST = paste(setdiff(expected, reported), collapse = ","), only_in_report = paste(setdiff(reported, expected), collapse = ","))
}
ky <- do.call(rbind, ky); cat("\nKEGG REST cross-check:\n"); print(ky, row.names = FALSE)
info <- tryCatch(readLines(url("https://rest.kegg.jp/info/kegg"), warn = FALSE)[2], error = function(e) NA)
cat("KEGG release line:", info, "\n")
write.csv(v, file.path(out, "enrichment_verification.csv"), row.names = FALSE)
write.csv(ky, file.path(out, "KEGG_REST_crosscheck.csv"), row.names = FALSE)
writeLines(paste("KEGG REST info (queried", format(Sys.Date()), "):", info), file.path(out, "KEGG_release.txt"))
if (fail > 0) stop(fail, " verification check(s) failed") else cat("\nALL ENRICHMENT CHECKS PASSED\n")
