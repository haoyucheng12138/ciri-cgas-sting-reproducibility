# =============================================================================
# 13_supplementary_tables.R  -  Supplementary Tables S1-S11 (one Excel workbook, one sheet per table),
# built only from GEO metadata (data/raw), processed data and results/. Nothing is typed by hand except the
# descriptive dataset fields of Table S1, which are copied from the GEO series records. Table S11 (qPCR primers of the
# experiments added in the revision) is read from reference/qPCR_primers_TableS11.csv, copied from the experimental protocol.
# =============================================================================
source("R/00_config.R")
suppressMessages(library(openxlsx))
rp <- function(step, f) file.path(DIR_RES, step, f)
soft_series <- function(acc) {
  s <- readLines(gzfile(file.path(DIR_RAW, paste0(acc, "_family.soft.gz"))), warn = FALSE)
  get <- function(k) paste(unique(sub(paste0("^!", k, " = "), "", grep(paste0("^!", k, " = "), s, value = TRUE))), collapse = "; ")
  list(title = get("Series_title"), platform = get("Series_platform_id"), pubmed = get("Series_pubmed_id"), n = length(grep("^\\^SAMPLE", s)))
}
sm <- readLines(gzfile(file.path(DIR_RAW, "GSE23163_series_matrix.txt.gz")))
gse23 <- list(title = gsub('"', "", sub("^!Series_title\t", "", grep("^!Series_title", sm, value = TRUE))),
              platform = gsub('"', "", sub("^!Series_platform_id\t", "", grep("^!Series_platform_id", sm, value = TRUE))),
              pubmed = gsub('"', "", sub("^!Series_pubmed_id\t", "", grep("^!Series_pubmed_id", sm, value = TRUE)[1])),
              n = length(strsplit(sub("^!Sample_geo_accession\t", "", grep("^!Sample_geo_accession", sm, value = TRUE)), "\t")[[1]]))
meta <- list(GSE131712 = soft_series("GSE131712"), GSE202659 = soft_series("GSE202659"), GSE131193 = soft_series("GSE131193"),
             GSE23163 = gse23, GSE227317 = soft_series("GSE227317"))
ann <- read.csv(file.path(DIR_PRO, "GSE131193_sample_annotation.csv"))
m23 <- read.csv(file.path(DIR_PRO, "GSE23163_sample_annotation.csv"))
S1 <- data.frame(Dataset = names(meta), GEO_title = sapply(meta, `[[`, "title"), Platform = sapply(meta, `[[`, "platform"),
                 PubMed_ID = sapply(meta, `[[`, "pubmed"), Samples_in_GEO = sapply(meta, `[[`, "n"),
                 Species = c("Mus musculus", "Mus musculus", "Mus musculus", "Mus musculus", "Rattus norvegicus"),
                 Tissue_per_GEO = c("ipsilateral cortex", "hemisphere", "brain microvessel fragments", "cortex and striatum", "cerebral cortex"),
                 Samples_used = c("3 Sham + 3 tMCAO 24 h; C57BL/6; all mice received the control vector AAV-SYN-shRNA-Con (titles Sham_C1-3, R24h_C1-3)",
                                  "3 Sham + 3 tMCAO; C57BL/6, male, 8-10 weeks; ischemic hemisphere (titles Sham_1-3, tMCAO_1-3); celastrol group excluded",
                                  sprintf("all %d; C57BL/6J, male, 10-12 weeks; tMCAO 60 min or sham; 24 h (1 d) or 7 d; ipsilateral and contralateral hemisphere", nrow(ann)),
                                  sprintf("%d of %d (Sham and I/R 24 h; WT and Gpx1-/-; cortex and striatum); SuperSeries, 2 h and 8 h samples excluded", sum(m23$group %in% c("Sham", "IR24h")), nrow(m23)),
                                  "5 MCAO 2 h/R 24 h + 5 Sham; Sprague-Dawley, male"),
                 Data_file = c("GSE131712_1_genes_fpkm_expression.txt.gz", "GSE202659_genes_fpkm_expression.txt.gz", "GSE131193_counts.matrix.norm_anno.txt.gz",
                               "GSE23163_series_matrix.txt.gz (+ GPL6885.annot.gz)", "GSE227317_mRNA.FPKM.txt.gz"),
                 Role = c("merged cohort (DEG, enrichment, PPI)", "merged cohort (DEG, enrichment, PPI)", "WGCNA, ssGSEA, IFN/STING score",
                          "Elastic Net hold-out tasks", "exploratory cross-species check"),
                 check.names = FALSE)
## Primary publication of each dataset. GEO lists a PubMed ID only for GSE23163; for the other series the publication was
## identified from the paper's own data-availability statement / text naming the accession (reference/dataset_primary_publications.csv).
pubs <- read.csv(file.path(DIR_REF, "dataset_primary_publications.csv"), check.names = FALSE)
k <- match(S1$Dataset, pubs$Dataset)
S1 <- data.frame(S1, Primary_publication_PMID = pubs$PMID[k], Primary_publication_DOI = pubs$DOI[k],
                 Primary_publication = pubs$Citation[k], How_linked = pubs$Verification[k], check.names = FALSE)
names(S1)[names(S1) == "PubMed_ID"] <- "PubMed_ID_in_GEO_record"

cm <- read.csv(file.path(DIR_PRO, "sample_column_map.csv"))
s2 <- read.csv(file.path(DIR_REF, "TableS2_merged_cohort_samples.csv"))
## bundled in reference/: produced once by the optional R/02_preprocess_check.R, which needs the original platform output
lab_audit <- read.csv(file.path(DIR_REF, "GSE131193_LKFlow_label_audit.csv"))
S2a <- merge(cm[cm$dataset != "GSE131193", ], s2, by.x = "GSM", by.y = "Sample_ID")[, c("dataset", "GSM", "geo_title", "file_column", "Group")]
names(S2a) <- c("Dataset", "GSM", "GEO_sample_title", "Column_in_GEO_file", "Group")
g131 <- merge(ann, cm[cm$dataset == "GSE131193", c("GSM", "file_column")], by.x = "ID", by.y = "GSM")
g131 <- merge(g131, lab_audit[, c("GSM", "LKFlow_label_on_same_data", "label_correct")], by.x = "ID", by.y = "GSM")
S2b <- data.frame(GSM = g131$ID, GEO_sample_title = g131$title, Injury = g131$injury, Time = g131$time, Side = g131$side,
                  Column_in_GEO_file = g131$file_column, In_primary_comparison = g131$side == "ipsilateral",
                  Original_analysis_label_on_these_data = g131$LKFlow_label_on_same_data, Original_label_correct = g131$label_correct)
S2b <- S2b[order(S2b$GSM), ]
S2c <- m23[m23$group %in% c("Sham", "IR24h"), c("ID", "title", "region", "genotype", "group")]; names(S2c) <- c("GSM", "GEO_sample_title", "Region", "Genotype", "Group")
soft <- readLines(gzfile(file.path(DIR_RAW, "GSE227317_family.soft.gz")), warn = FALSE)
S2d <- data.frame(GSM = sub("^\\^SAMPLE = ", "", grep("^\\^SAMPLE", soft, value = TRUE)), GEO_sample_title = sub("^!Sample_title = ", "", grep("^!Sample_title", soft, value = TRUE)))
S2d$Column_in_GEO_file <- sub(".*\\[(T\\d+)\\]$", "\\1", S2d$GEO_sample_title); S2d$Group <- ifelse(grepl("MCAO", S2d$GEO_sample_title), "MCAO", "Sham")

A <- read.csv(rp("03_combat_limma", "limma_primary_ComBat_withGroup_all_genes.csv"), row.names = 1)
B <- read.csv(rp("03_combat_limma", "limma_sensitivity_B.csv"), row.names = 1); C <- read.csv(rp("03_combat_limma", "limma_sensitivity_C.csv"), row.names = 1)
E1 <- read.csv(rp("03_combat_limma", "limma_within_GSE131712.csv"), row.names = 1); E2 <- read.csv(rp("03_combat_limma", "limma_within_GSE202659.csv"), row.names = 1)
dg <- rownames(A)[A$DEG]
S3 <- data.frame(Gene = dg, log2FC = A[dg, "logFC"], AveExpr = A[dg, "AveExpr"], t = A[dg, "t"], P = A[dg, "P.Value"], FDR = A[dg, "adj.P.Val"],
                 Direction = ifelse(A[dg, "logFC"] > 0, "up", "down"),
                 FDR_ComBat_no_covariate = B[dg, "adj.P.Val"], DEG_ComBat_no_covariate = B[dg, "DEG"],
                 FDR_batch_in_model = C[dg, "adj.P.Val"], DEG_batch_in_model = C[dg, "DEG"],
                 log2FC_GSE131712_alone = E1[dg, "logFC"], P_GSE131712_alone = E1[dg, "P.Value"],
                 log2FC_GSE202659_alone = E2[dg, "logFC"], P_GSE202659_alone = E2[dg, "P.Value"])
S3 <- S3[order(S3$FDR), ]
S4 <- read.csv(rp("03_combat_limma", "R2.5_batch_sensitivity_summary.csv"))[, c("strategy", "n_DEG", "up", "down", "shared_with_primary")]
st <- read.csv(rp("03_combat_limma", "STING_score_merged_cohort.csv"))

kg <- read.csv(rp("04_enrichment_ppi", "KEGG.csv")); cols <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count", "geneID")
S5 <- rbind(cbind(Source = "KEGG", kg[, cols]),
            do.call(rbind, lapply(c("BP", "CC", "MF"), function(o) cbind(Source = paste("GO", o), read.csv(rp("04_enrichment_ppi", paste0("GO_", o, ".csv")))[, cols]))))
names(S5)[names(S5) == "p.adjust"] <- "FDR_BH"
S5b <- read.csv(rp("04_enrichment_ppi", "PPI_nodes_score700.csv"))
S5c <- read.csv(rp("04_enrichment_ppi", "DEG_to_Entrez_map.csv"))

mt <- read.csv(rp("05_wgcna", "module_trait_r.csv"), check.names = FALSE); mp <- read.csv(rp("05_wgcna", "module_trait_p.csv"), check.names = FALSE)
S6a <- merge(mt, mp, by = "module", suffixes = c("_r", "_p"))
S6a$ipsi_tMCAO_vs_ipsi_sham_FDR <- p.adjust(S6a$ipsi_tMCAO_vs_ipsi_sham_p, "BH")   # BH across the 23 module eigengenes (incl. grey), as in 05_wgcna.R
S6a <- S6a[order(S6a$ipsi_tMCAO_vs_ipsi_sham_p), ]
S6b <- read.csv(rp("05_wgcna", "focal_genes.csv"))
S6c <- read.csv(rp("05_wgcna", "soft_threshold.csv"))

gt <- read.csv(rp("06_ssgsea_sting", "genes_and_STING_score_group_tests.csv"))
ct <- read.csv(rp("06_ssgsea_sting", "ssGSEA_group_tests.csv")); names(ct)[names(ct) == "FDR_wilcox"] <- "FDR_BH"
cr <- read.csv(rp("06_ssgsea_sting", "STING_score_vs_ssGSEA_spearman.csv"))
gmt <- strsplit(readLines(file.path(DIR_REF, "Charoentong28_mouse.gmt")), "\t")
cov <- read.csv(file.path(DIR_REF, "Charoentong28_mouse_coverage.csv"))
S7d <- data.frame(Cell_type = gsub("_", " ", sapply(gmt, `[`, 1)), Mouse_genes = sapply(gmt, function(x) paste(x[-(1:2)], collapse = ", ")))
S7d <- merge(S7d, cov, by.x = "Cell_type", by.y = "cell_type")

ev <- read.csv(rp("07_intersection", "intersection_evidence.csv")); ic <- read.csv(rp("07_intersection", "intersection_counts.csv"))

en <- read.csv(rp("08_elastic_net", "elastic_net_summary.csv"))
coefs <- do.call(rbind, lapply(seq_len(nrow(en)), function(i) { x <- read.csv(rp("08_elastic_net", sprintf("coefficients_%s_%s.csv", en$feature_space[i], en$task[i])))
  cbind(Feature_space = en$feature_space[i], Task = en$task[i], x) }))
rat <- read.csv(rp("09_rat_roc", "rat_ROC.csv"))

wb <- createWorkbook(); hs <- createStyle(textDecoration = "bold", fgFill = "#DDEBF7", wrapText = TRUE)
put <- function(sheet, title, df_list, notes = NULL) {
  addWorksheet(wb, sheet); writeData(wb, sheet, title, startRow = 1); addStyle(wb, sheet, createStyle(textDecoration = "bold"), rows = 1, cols = 1)
  r <- 2; if (!is.null(notes)) { writeData(wb, sheet, notes, startRow = 2); r <- 2 + length(notes) }
  for (nm in names(df_list)) { r <- r + 1; if (nzchar(nm)) { writeData(wb, sheet, nm, startRow = r); addStyle(wb, sheet, createStyle(textDecoration = "italic"), rows = r, cols = 1); r <- r + 1 }
    writeData(wb, sheet, df_list[[nm]], startRow = r, headerStyle = hs); r <- r + nrow(df_list[[nm]]) + 1 }
  setColWidths(wb, sheet, cols = 1:30, widths = "auto")
}
put("S1_Datasets", "Table S1. GEO datasets analysed (title, platform, PubMed ID and sample numbers read from the GEO records) and the primary publication of each dataset.", list(` ` = S1))
put("S2_Samples", "Table S2. Samples and their labels. Columns of each GEO file were assigned to GSM accessions by sample title.",
    list(`A. Merged cohort (GSE131712, GSE202659)` = S2a, `B. GSE131193 (brain microvessel fragments); the last two columns document the labelling error of the original analysis` = S2b,
         `C. GSE23163 samples used for Elastic Net` = S2c, `D. GSE227317` = S2d))
put("S3_DEGs", "Table S3. Differentially expressed genes of the merged cohort (tMCAO vs Sham, 6 vs 6; ComBat with group covariate, limma; FDR < 0.05 and |log2FC| > 0.585), with the batch-correction sensitivity analyses.",
    list(` ` = S3))
put("S4_Batch_sensitivity", "Table S4. Batch-correction sensitivity analyses and composite IFN/STING score in the merged cohort.",
    list(`A. DEG counts by strategy` = S4, `B. IFN/STING score (z-mean of Mb21d1, Tmem173, Irf3, Isg15, Cxcl10, Trex1)` = st))
put("S5_Enrichment_PPI", "Table S5. GO/KEGG enrichment of the DEGs (clusterProfiler; BH FDR < 0.05) and STRING v11.5 network (combined score >= 0.7).",
    list(`A. Significant terms` = S5, `B. PPI nodes` = S5b, `C. DEG symbol to Entrez mapping` = S5c),
    notes = c("All tested terms (including non-significant) are in results/04_enrichment_ppi/*_all_tested.csv of the code repository."))
put("S6_WGCNA", "Table S6. WGCNA of GSE131193 (24 samples; unsigned network; soft power 12).",
    list(`A. Module-trait correlations (Pearson r and p; FDR = Benjamini-Hochberg adjustment of the primary-trait p values across the 23 module eigengenes)` = S6a, `B. Module membership of selected genes` = S6b, `C. Soft-threshold scan` = S6c))
put("S7_Immune_STING", "Table S7. GSE131193 immune signatures (ssGSEA, Charoentong et al. 2017, 28 cell types) and composite IFN/STING score (type I IFN/STING-associated score).",
    list(`A. Gene and IFN/STING-score group tests (exact Wilcoxon; Welch t)` = gt, `B. ssGSEA group tests` = ct,
         `C. Spearman correlations of IFN/STING-score versions with ssGSEA scores` = cr, `D. Mouse immune gene sets` = S7d))
put("S8_Intersection", "Table S8. Three-set intersection: DEGs, WGCNA core genes and curated cGAS-STING pathway genes.", list(`A. Counts` = ic, `B. Genes` = ev))
seed <- read.csv(rp("08_elastic_net", "seed_sensitivity_summary.csv"))
put("S9_Elastic_net", "Table S9. Elastic Net logistic regression in GSE23163 (alpha 0.4, 10-fold CV, lambda.min, seed 12345). Predictors were standardized internally (glmnet default); coefficients are reported on the original log2 expression scale, as returned by coef().",
    list(`A. Performance on the held-out test sets (pre-specified seed 12345)` = en,
         `B. Dependence on the random CV fold assignment (seeds 1-200; full transcriptome seeds 1-30)` = seed,
         `C. Non-zero coefficients (seed 12345)` = coefs),
    notes = c("With 16 training samples, 10-fold CV leaves 1-2 samples per fold; glmnet therefore uses ungrouped CV error (grouped = FALSE) and warns that one class has fewer than 8 observations."))
put("S10_Rat", "Table S10. GSE227317 rat cortex (5 MCAO 2 h/R 24 h vs 5 Sham): AUC with 2,000 stratified bootstrap 95% CI and exact Wilcoxon test.", list(` ` = rat))
pr <- read.csv(file.path(DIR_REF, "qPCR_primers_TableS11.csv"), check.names = FALSE, fileEncoding = "UTF-8")
put("S11_qPCR_primers", "Table S11. Primer sequences used for RT-qPCR of mouse brain tissue (Supplementary Fig. S4A). Tmem173 encodes STING (official symbol Sting1); beta-actin was the reference gene.", list(` ` = pr))
out <- file.path(DIR_RES, "Supplementary_Tables_S1-S11.xlsx"); saveWorkbook(wb, out, overwrite = TRUE)
cat("written", out, "\n"); for (s in names(wb)) cat(sprintf("  %-22s\n", s))
cat(sprintf("rows: S3 %d | S5 terms %d, nodes %d | S6 modules %d | S7 cell tests %d | S9 coefficients %d | S10 %d\n", nrow(S3), nrow(S5), nrow(S5b), nrow(S6a), nrow(ct), nrow(coefs), nrow(rat)))
