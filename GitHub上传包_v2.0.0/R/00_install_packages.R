# =============================================================================
# 00_install_packages.R  -  install the packages used by the pipeline (run once).
# Versions used for the reported results (R 4.6.1, Windows x64) are listed in README.md and in
# results/<step>/sessionInfo.txt.
# =============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))
cran <- c("BiocManager", "glmnet", "pROC", "WGCNA", "httr", "igraph")
bioc <- c("limma", "sva", "GSVA", "clusterProfiler", "org.Mm.eg.db", "GO.db", "impute", "preprocessCore", "AnnotationDbi")
for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE, ask = FALSE)
for (p in c(cran, bioc)) cat(sprintf("%-16s %s\n", p, tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")))
