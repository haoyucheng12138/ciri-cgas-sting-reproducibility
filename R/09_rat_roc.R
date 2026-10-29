# =============================================================================
# 09_rat_roc.R  -  exploratory cross-species consistency check in GSE227317 (rat cortex, 5 MCAO 2h/R24h vs 5 Sham).
#  Labels from the GEO SOFT file (T01-T05 = MCAO, T06-T10 = Sham).
#  Transcript-level FPKM -> gene symbol via Ensembl BioMart (Rattus norvegicus); gene FPKM = sum of transcripts;
#  log2(FPKM + 1). Panel score = z-mean of Tap1, Tap2, Cxcl10, Psmb8, Stat1 (rat genes with the same symbol).
#  AUC with 2,000 stratified bootstrap 95% CI (pROC); exact Wilcoxon test. With 5 vs 5, AUC = 1 gives a
#  saturated, uninformative CI; this analysis is a directional check only.
# =============================================================================
source("R/00_config.R")
suppressMessages({ library(pROC); library(httr) })
out <- res_dir("09_rat_roc")
soft <- readLines(gzfile(file.path(DIR_RAW, "GSE227317_family.soft.gz")), warn = FALSE)
ttl <- sub("^!Sample_title = ", "", grep("^!Sample_title", soft, value = TRUE))
lab <- data.frame(sample = sub(".*\\[(T\\d+)\\]$", "\\1", ttl), MCAO = as.integer(grepl("MCAO", ttl)))
print(lab)
fp <- read.delim(gzfile(file.path(DIR_RAW, "GSE227317_mRNA.FPKM.txt.gz")), check.names = FALSE); names(fp)[1] <- "tx"
targets <- c("Tap1", "Tap2", "Cxcl10", "Psmb8", "Stat1", "Psmb10", "Trpc6")
## Transcript -> gene map. The BioMart response is stored in reference/ so the analysis does not depend on the
## server being available; it is (re)downloaded only if missing, and an HTML error page is rejected.
map_file <- file.path(DIR_REF, "GSE227317_Ensembl_transcript_map.tsv")
if (!file.exists(map_file)) {
  q <- paste0('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1">',
              '<Dataset name="rnorvegicus_gene_ensembl" interface="default"><Filter name="external_gene_name" value="', paste(targets, collapse = ","), '"/>',
              '<Attribute name="ensembl_transcript_id"/><Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/></Dataset></Query>')
  ok <- FALSE
  for (h in c("https://www.ensembl.org", "https://useast.ensembl.org", "https://asia.ensembl.org")) {
    r <- tryCatch(GET(paste0(h, "/biomart/martservice"), query = list(query = q), timeout(180)), error = function(e) NULL)
    if (is.null(r) || status_code(r) != 200) next
    txt <- content(r, "text", encoding = "UTF-8")
    if (!grepl("^Transcript stable ID", txt)) next
    writeLines(txt, map_file); ok <- TRUE; break
  }
  stopifnot(ok)
  rel <- tryCatch(content(GET("https://rest.ensembl.org/info/data/?content-type=application/json", timeout(60)))$releases[[1]], error = function(e) NA)
  writeLines(paste("Ensembl release:", rel), file.path(DIR_REF, "GSE227317_Ensembl_release.txt"))
}
map <- read.delim(map_file); stopifnot(ncol(map) == 3); names(map) <- c("tx", "gene_id", "symbol")
hit <- merge(map, fp, by = "tx"); smp <- lab$sample
print(table(hit$symbol))
gene <- t(sapply(targets, function(g) colSums(hit[hit$symbol == g, smp, drop = FALSE])))
write.csv(gene, file.path(out, "gene_FPKM.csv"))
lg <- log2(gene + 1); y <- setNames(lab$MCAO, smp)
panel <- zmean(lg[c("Tap1", "Tap2", "Cxcl10", "Psmb8", "Stat1"), ])
res <- do.call(rbind, lapply(c(targets, "Panel5_zmean"), function(g) {
  v <- if (g == "Panel5_zmean") panel else lg[g, ]
  r <- roc(y, v, levels = c(0, 1), direction = "<", quiet = TRUE)
  set.seed(SEED); ci <- ci.auc(r, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE, progress = "none")
  data.frame(feature = g, AUC = round(as.numeric(auc(r)), 3), CI_low = round(ci[1], 3), CI_high = round(ci[3], 3),
             mean_MCAO = round(mean(v[y == 1]), 3), mean_Sham = round(mean(v[y == 0]), 3), wilcox_exact_p = signif(wilcox.test(v[y == 1], v[y == 0], exact = TRUE)$p.value, 3))
}))
print(res, row.names = FALSE); write.csv(res, file.path(out, "rat_ROC.csv"), row.names = FALSE)
save_session("09_rat_roc")
