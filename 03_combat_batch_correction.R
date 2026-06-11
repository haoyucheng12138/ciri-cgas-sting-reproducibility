# =============================================================================
# 03_combat_batch_correction.R
# - Build merged validation cohort from GSE131712 + GSE202659
#   (n = 6 Sham, n = 6 Stroke; sample IDs as listed in Supplementary Table S2)
# - Apply parametric ComBat with biological group preserved as covariate
# - Save before/after PCA, UMAP, and sample expression boxplot diagnostics
# =============================================================================

source("00_setup.R")

suppressPackageStartupMessages({
  library(sva)
  library(umap)
  library(ggplot2)
  library(dplyr)
  library(tibble)
})

# ----------------------------------------------------------------------------
# Sample composition of the merged validation cohort
# (mirrors Supplementary Table S2)
# ----------------------------------------------------------------------------
cohort_samples <- tribble(
  ~sample_id,    ~dataset,    ~group,
  "GSM3814504",  "GSE131712", "Sham",
  "GSM3814505",  "GSE131712", "Sham",
  "GSM3814506",  "GSE131712", "Sham",
  "GSM3814510",  "GSE131712", "Stroke",
  "GSM3814511",  "GSE131712", "Stroke",
  "GSM3814512",  "GSE131712", "Stroke",
  "GSM6127878",  "GSE202659", "Sham",
  "GSM6127879",  "GSE202659", "Sham",
  "GSM6127880",  "GSE202659", "Sham",
  "GSM6127881",  "GSE202659", "Stroke",
  "GSM6127882",  "GSE202659", "Stroke",
  "GSM6127883",  "GSE202659", "Stroke"
)

# ----------------------------------------------------------------------------
# Load preprocessed expression matrices
# ----------------------------------------------------------------------------
load_expr <- function(acc) {
  readRDS(file.path(DIR_DATA, acc, paste0(acc, "_expr.rds")))
}
expr1 <- load_expr("GSE131712")
expr2 <- load_expr("GSE202659")

# Restrict to genes present in both matrices
common_genes <- intersect(rownames(expr1), rownames(expr2))
if (length(common_genes) == 0) {
  stop("No overlapping gene symbols between GSE131712 and GSE202659.")
}

expr1 <- expr1[common_genes, intersect(colnames(expr1), cohort_samples$sample_id), drop = FALSE]
expr2 <- expr2[common_genes, intersect(colnames(expr2), cohort_samples$sample_id), drop = FALSE]

merged_expr <- cbind(expr1, expr2)
merged_expr <- merged_expr[, cohort_samples$sample_id, drop = FALSE]

batch <- factor(cohort_samples$dataset)
group <- factor(cohort_samples$group, levels = c("Sham", "Stroke"))

# ----------------------------------------------------------------------------
# ComBat: preserve biological group as covariate
# ----------------------------------------------------------------------------
mod <- model.matrix(~ group)
set.seed(GLOBAL_SEED)
combat_expr <- ComBat(
  dat       = as.matrix(merged_expr),
  batch     = batch,
  mod       = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

saveRDS(merged_expr, file.path(DIR_DATA, "merged_validation_expr_raw.rds"))
saveRDS(combat_expr, file.path(DIR_DATA, "merged_validation_expr_combat.rds"))
saveRDS(cohort_samples, file.path(DIR_DATA, "merged_validation_pheno.rds"))

# ----------------------------------------------------------------------------
# Diagnostic plots
# ----------------------------------------------------------------------------
make_pca_plot <- function(mat, title, color_var) {
  pca <- prcomp(t(mat), scale. = TRUE)
  df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    sample_id = colnames(mat),
    stringsAsFactors = FALSE
  )
  df <- merge(df, cohort_samples, by = "sample_id")
  df$color <- df[[color_var]]
  ggplot(df, aes(PC1, PC2, color = color, shape = group)) +
    geom_point(size = 4) +
    labs(title = title, color = color_var) +
    theme_minimal(base_size = 13)
}

p_before <- make_pca_plot(merged_expr, "PCA: before ComBat", "dataset")
p_after  <- make_pca_plot(combat_expr, "PCA: after ComBat",  "dataset")

ggsave(file.path(DIR_FIG, "Fig2A_PCA_before_combat.pdf"), p_before, width = 6, height = 5)
ggsave(file.path(DIR_FIG, "Fig2B_PCA_after_combat.pdf"),  p_after,  width = 6, height = 5)

# UMAP after correction
set.seed(GLOBAL_SEED)
umap_res <- umap(t(combat_expr), n_neighbors = 4, min_dist = 0.3)
umap_df <- data.frame(
  UMAP1 = umap_res$layout[, 1],
  UMAP2 = umap_res$layout[, 2],
  sample_id = colnames(combat_expr)
)
umap_df <- merge(umap_df, cohort_samples, by = "sample_id")

p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = group, shape = dataset)) +
  geom_point(size = 4) +
  labs(title = "UMAP: after ComBat") +
  theme_minimal(base_size = 13)
ggsave(file.path(DIR_FIG, "Fig2D_UMAP_after_combat.pdf"), p_umap, width = 6, height = 5)

# Sample-level expression distribution boxplot
make_box_df <- function(mat, label) {
  df <- as.data.frame(mat) |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "sample_id", values_to = "expr")
  df$stage <- label
  df
}
df_box <- rbind(make_box_df(merged_expr, "Before ComBat"),
                make_box_df(combat_expr, "After ComBat"))
df_box <- merge(df_box, cohort_samples, by = "sample_id")

p_box <- ggplot(df_box, aes(sample_id, expr, fill = dataset)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ stage, ncol = 1, scales = "free") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sample expression distribution before/after ComBat")
ggsave(file.path(DIR_FIG, "FigS1_expr_boxplot_combat.pdf"), p_box, width = 8, height = 6)

cat("\nComBat batch correction complete.\n")
cat("Outputs saved to data/ and figures/.\n")
