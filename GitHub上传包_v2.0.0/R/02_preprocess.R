# =============================================================================
# 02_preprocess.R  -  build every analysis matrix from the raw GEO files.
# Sample columns are assigned to GEO accessions BY NAME, using the sample titles in the GEO SOFT files
# (never by position). The mapping is written to data/processed/sample_column_map.csv.
#  GSE131712, GSE202659 : gene-level FPKM (symbol column "gene_name"); duplicate symbols -> median;
#                         log2(FPKM + 1); keep genes with mean log2 > 0.1
#  merged cohort        : genes common to both filtered matrices (6 Sham + 6 tMCAO, Table S2)
#  GSE131193            : normalised counts supplied by the submitters (symbol column "Ensembl gene");
#                         Excel-converted symbols restored (e.g. "1-Mar" -> "March1"); duplicates -> median;
#                         log2(x + 1); keep genes with mean log2 > 1.
#                         NOTE: the columns of this file are grouped by condition and are NOT in accession order.
#  GSE23163             : GeneSpring median-normalised ratios (series matrix; columns are GSM accessions);
#                         log2(x + 1) per probe; probes -> gene symbol (GPL6885) by median; keep mean >= 1
# =============================================================================
source("R/00_config.R")
rd <- function(f) read.delim(gzfile(file.path(DIR_RAW, f)), check.names = FALSE)
soft_titles <- function(acc) {
  s <- readLines(gzfile(file.path(DIR_RAW, paste0(acc, "_family.soft.gz"))), warn = FALSE)
  data.frame(ID = sub("^\\^SAMPLE = ", "", grep("^\\^SAMPLE", s, value = TRUE)),
             title = sub("^!Sample_title = ", "", grep("^!Sample_title", s, value = TRUE)))
}
colmap <- list()

## ---- merged validation cohort ---------------------------------------------
s2 <- read.csv(file.path(DIR_REF, "TableS2_merged_cohort_samples.csv"))
## Gene-symbol harmonisation across datasets. GSE202659 uses the current symbol "Cgas" for ENSMUSG00000032344,
## whereas GSE131712, GSE131193 and GSE23163 use the former symbol "Mb21d1"; the former symbol is used throughout
## (matching the manuscript). Excel-converted symbols ("1-Mar" ...) are restored to "March1" ... .
harmonise <- function(sym, ensembl = NULL) {
  dm <- grepl("^\\d+-Mar$", sym); sym[dm] <- paste0("March", sub("-Mar$", "", sym[dm]))
  if (!is.null(ensembl)) sym[ensembl == "ENSMUSG00000032344"] <- "Mb21d1"
  sym
}
fpkm <- function(acc, f) {
  r <- rd(f); st <- soft_titles(acc); ids <- s2$Sample_ID[s2$Dataset == acc]
  ttl <- st$title[match(ids, st$ID)]; cols <- paste0("FPKM.", ttl)
  stopifnot(!anyNA(ttl), all(cols %in% colnames(r)))
  colmap[[acc]] <<- data.frame(dataset = acc, GSM = ids, geo_title = ttl, file_column = cols)
  m <- as.matrix(r[, cols]); colnames(m) <- ids
  g <- log2(collapse_median(m, harmonise(r$gene_name, r$gene_id)) + 1)
  g[rowMeans(g) > 0.1, , drop = FALSE]
}
e1 <- fpkm("GSE131712", "GSE131712_1_genes_fpkm_expression.txt.gz")
e2 <- fpkm("GSE202659", "GSE202659_genes_fpkm_expression.txt.gz")
common <- intersect(rownames(e1), rownames(e2))
merged <- cbind(e1[common, ], e2[common, ])[, s2$Sample_ID]
write.csv(merged, file.path(DIR_PRO, "merged_log2FPKM.csv"))
cat(sprintf("GSE131712 %d genes | GSE202659 %d genes | merged %d genes x %d samples\n", nrow(e1), nrow(e2), nrow(merged), ncol(merged)))

## ---- GSE131193 -------------------------------------------------------------
ann <- soft_titles("GSE131193")
ann$injury <- ifelse(grepl("^tMCAO", ann$title), "tMCAO", "Sham")
ann$time   <- ifelse(grepl("_1d_", ann$title), "1d", "7d")
ann$side   <- ifelse(grepl("ipsilateral$", ann$title), "ipsilateral", "contralateral")
write.csv(ann, file.path(DIR_PRO, "GSE131193_sample_annotation.csv"), row.names = FALSE)
r <- rd("GSE131193_counts.matrix.norm_anno.txt.gz")
pat <- "^(Sham|Stroke)-(1d|7d)-(contralateral|ipsilateral)_([1-3])$"
cols <- grep(pat, colnames(r), value = TRUE); stopifnot(length(cols) == 24)
z <- regmatches(cols, regexec(pat, cols))
ttl <- vapply(z, function(v) sprintf("%s_%s_%s_%s", ifelse(v[2] == "Stroke", "tMCAO", "sham"), v[3], v[5], v[4]), "")
gsm <- ann$ID[match(ttl, ann$title)]; stopifnot(!anyNA(gsm), !anyDuplicated(gsm))
colmap$GSE131193 <- data.frame(dataset = "GSE131193", GSM = gsm, geo_title = ttl, file_column = cols)
m <- as.matrix(r[, cols]); colnames(m) <- gsm; m <- m[, ann$ID]
sym <- harmonise(r[["Ensembl gene"]], r[["Ensembl gene id"]])
g <- log2(collapse_median(m, sym) + 1); g <- g[rowMeans(g) > 1, ]
write.csv(g, file.path(DIR_PRO, "GSE131193_log2counts.csv"))
cat(sprintf("GSE131193 %d genes x %d samples | file columns not in accession order: %d of 24 (assigned by name)\n",
            nrow(g), ncol(g), sum(gsm != ann$ID)))
grp <- paste(ann$injury, ann$time, ann$side, sep = "_")
cat("sanity check, mean log2 of injury markers by GEO group:\n")
print(round(sapply(split(ann$ID, grp), function(s) rowMeans(g[c("Hmox1", "Timp1", "Cd44"), s, drop = FALSE])), 2))

## ---- GSE23163 --------------------------------------------------------------
sm <- readLines(gzfile(file.path(DIR_RAW, "GSE23163_series_matrix.txt.gz")))
b0 <- grep("^!series_matrix_table_begin", sm); e0 <- grep("^!series_matrix_table_end", sm)
x <- read.delim(text = sm[(b0 + 1):(e0 - 1)], check.names = FALSE); rownames(x) <- x$ID_REF; x <- as.matrix(x[, -1])
an <- readLines(gzfile(file.path(DIR_RAW, "GPL6885.annot.gz"))); hb <- grep("^ID\t", an)[1]
pa <- read.delim(text = an[hb:length(an)], check.names = FALSE, quote = "")
sym <- pa[["Gene symbol"]][match(rownames(x), pa$ID)]
## Symbol clean-up (reproduces the original matrix exactly):
##  - multi-mapping probe "Lgals6///Lgals4" keeps its first symbol (Lgals6);
##  - single-character symbols ("a", "T") and symbols containing "." (e.g. H2-M10.3, Tex19.2) are dropped.
sym <- sub("///.*$", "", sym)
sym[!is.na(sym) & (nchar(sym) == 1 | grepl(".", sym, fixed = TRUE))] <- NA
g <- collapse_median(log2(x + 1), sym)
g <- g[rowMeans(g) >= 1, ]
write.csv(g, file.path(DIR_PRO, "GSE23163_log2ratio.csv"))
ttl <- gsub('"', "", strsplit(sub("^!Sample_title\t", "", grep("^!Sample_title", sm, value = TRUE)), "\t")[[1]])
acc <- gsub('"', "", strsplit(sub("^!Sample_geo_accession\t", "", grep("^!Sample_geo_accession", sm, value = TRUE)), "\t")[[1]])
stopifnot(identical(acc, colnames(x)))
meta <- data.frame(ID = acc, title = ttl)
meta$region   <- ifelse(grepl("^Cortex", meta$title), "Cortex", ifelse(grepl("^Striatum", meta$title), "Striatum", NA))
meta$genotype <- ifelse(grepl("Gpx1-/-", meta$title, fixed = TRUE), "Gpx1KO", "WT")
meta$group    <- ifelse(grepl("Sham", meta$title), "Sham", ifelse(grepl("I/R-24h", meta$title, fixed = TRUE), "IR24h", "other"))
write.csv(meta, file.path(DIR_PRO, "GSE23163_sample_annotation.csv"), row.names = FALSE)
cat(sprintf("GSE23163 %d genes x %d samples\n", nrow(g), ncol(g)))
write.csv(do.call(rbind, colmap), file.path(DIR_PRO, "sample_column_map.csv"), row.names = FALSE)
save_session("02_preprocess")
