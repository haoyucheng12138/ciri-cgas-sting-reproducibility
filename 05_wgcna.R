# =============================================================================
# 05_wgcna.R
# WGCNA on the GSE131193 discovery cohort.
#   - signed network
#   - soft-thresholding power = 9 (chosen via scale-free topology, R^2 > 0.85)
#   - minimum module size = 30
#   - deepSplit = 2
#   - merge cut height = 0.3
# Module eigengenes are correlated with the binary ischemia trait.
# =============================================================================

source("00_setup.R")

suppressPackageStartupMessages({
  library(WGCNA)
  library(dplyr)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

expr  <- readRDS(file.path(DIR_DATA, "GSE131193", "GSE131193_expr.rds"))
pheno <- readRDS(file.path(DIR_DATA, "GSE131193", "GSE131193_pheno.rds"))

# ----------------------------------------------------------------------------
# Trait vector: 1 = ischemia, 0 = sham
# Edit the regular expression below if your phenoData uses different labels.
# ----------------------------------------------------------------------------
group_field <- grep("characteristics_ch1|title|source_name",
                    colnames(pheno), value = TRUE, ignore.case = TRUE)
if (length(group_field) == 0) {
  stop("Cannot find a phenoData field encoding the experimental group. ",
       "Inspect pheno columns and edit 05_wgcna.R accordingly.")
}
group_label <- as.character(pheno[, group_field[1]])
trait <- as.integer(grepl("(?i)ischemia|stroke|MCAO", group_label, perl = TRUE))
names(trait) <- rownames(pheno)
trait <- trait[colnames(expr)]

# ----------------------------------------------------------------------------
# Filter to most variable genes for tractability
# ----------------------------------------------------------------------------
gene_var <- apply(expr, 1, var)
top_n    <- min(8000, nrow(expr))
expr_top <- expr[order(-gene_var)[seq_len(top_n)], , drop = FALSE]
datExpr  <- t(expr_top)

# ----------------------------------------------------------------------------
# Soft-thresholding power selection (kept here for transparency; we use
# power = 9 as reported in the manuscript Methods).
# ----------------------------------------------------------------------------
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = "signed", verbose = 0)
write.csv(sft$fitIndices,
          file.path(DIR_RESULTS, "WGCNA_soft_threshold_fit.csv"),
          row.names = FALSE)

soft_power <- 9

# ----------------------------------------------------------------------------
# Block-wise network construction
# ----------------------------------------------------------------------------
set.seed(GLOBAL_SEED)
net <- blockwiseModules(
  datExpr,
  power           = soft_power,
  networkType     = "signed",
  TOMType         = "signed",
  minModuleSize   = 30,
  deepSplit       = 2,
  mergeCutHeight  = 0.3,
  numericLabels   = FALSE,
  saveTOMs        = FALSE,
  pamRespectsDendro = FALSE,
  verbose         = 0
)

module_colors <- net$colors
MEs <- net$MEs
saveRDS(net, file.path(DIR_RESULTS, "WGCNA_network.rds"))

# ----------------------------------------------------------------------------
# Module-trait correlation
# ----------------------------------------------------------------------------
moduleTraitCor   <- cor(MEs, trait, use = "p")
moduleTraitPval  <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

mt_table <- data.frame(
  module      = rownames(moduleTraitCor),
  correlation = moduleTraitCor[, 1],
  pvalue      = moduleTraitPval[, 1]
)
write.csv(mt_table,
          file.path(DIR_RESULTS, "WGCNA_module_trait.csv"),
          row.names = FALSE)

# Module membership table for downstream inspection
gene_module <- data.frame(
  gene   = colnames(datExpr),
  module = module_colors
)
write.csv(gene_module,
          file.path(DIR_RESULTS, "WGCNA_gene_module_assignment.csv"),
          row.names = FALSE)

cat("\nWGCNA complete.\n")
cat(sprintf("Number of modules detected: %d\n", length(unique(module_colors))))
