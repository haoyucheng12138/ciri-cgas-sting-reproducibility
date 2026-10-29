# =============================================================================
# 01_download_geo.R  -  download the public GEO files used in this study.
# Every file is fetched from the official NCBI GEO FTP; MD5 sums are written to
# data/raw/MD5SUMS.txt so readers can confirm they used identical inputs.
# =============================================================================
source("R/00_config.R")
ftp <- "https://ftp.ncbi.nlm.nih.gov/geo"
files <- c(
  GSE131193_counts.matrix.norm_anno.txt.gz = "series/GSE131nnn/GSE131193/suppl/GSE131193_counts.matrix.norm_anno.txt.gz",
  GSE131193_family.soft.gz                 = "series/GSE131nnn/GSE131193/soft/GSE131193_family.soft.gz",
  GSE131712_1_genes_fpkm_expression.txt.gz = "series/GSE131nnn/GSE131712/suppl/GSE131712_1_genes_fpkm_expression.txt.gz",
  GSE131712_family.soft.gz                 = "series/GSE131nnn/GSE131712/soft/GSE131712_family.soft.gz",
  GSE202659_genes_fpkm_expression.txt.gz   = "series/GSE202nnn/GSE202659/suppl/GSE202659_genes_fpkm_expression.txt.gz",
  GSE202659_family.soft.gz                 = "series/GSE202nnn/GSE202659/soft/GSE202659_family.soft.gz",
  GSE23163_series_matrix.txt.gz            = "series/GSE23nnn/GSE23163/matrix/GSE23163_series_matrix.txt.gz",
  GPL6885.annot.gz                         = "platforms/GPL6nnn/GPL6885/annot/GPL6885.annot.gz",
  GSE227317_mRNA.FPKM.txt.gz               = "series/GSE227nnn/GSE227317/suppl/GSE227317_mRNA.FPKM.txt.gz",
  GSE227317_family.soft.gz                 = "series/GSE227nnn/GSE227317/soft/GSE227317_family.soft.gz")
for (f in names(files)) {
  dest <- file.path(DIR_RAW, f)
  if (!file.exists(dest)) { cat("downloading", f, "\n"); download.file(file.path(ftp, files[[f]]), dest, mode = "wb", quiet = TRUE) }
}
md5 <- tools::md5sum(file.path(DIR_RAW, names(files)))
writeLines(sprintf("%s  %s", md5, basename(names(md5))), file.path(DIR_RAW, "MD5SUMS.txt"))
cat("downloaded", length(files), "files; MD5 sums in data/raw/MD5SUMS.txt\n")
