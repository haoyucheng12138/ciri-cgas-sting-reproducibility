# =============================================================================
# 04_limma_deg.R
# Differential expression analysis on the ComBat-corrected merged validation
# cohort (Stroke vs Sham). Thresholds: FDR < 0.05, |log2FC| > 0.585.
# =============================================================================

source("00_setup.R")

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

expr  <- readRDS(file.path(DIR_DATA, "merged_validation_expr_combat.rds"))
pheno <- readRDS(file.path(DIR_DATA, "merged_validation_pheno.rds"))
stopifnot(all(colnames(expr) == pheno$sample_id))

group <- factor(pheno$group, levels = c("Sham", "Stroke"))

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

contrast_matrix <- makeContrasts(StrokeVsSham = Stroke - Sham, levels = design)

fit  <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

deg_table <- topTable(fit2, coef = "StrokeVsSham", number = Inf, adjust.method = "BH") |>
  tibble::rownames_to_column("gene") |>
  mutate(
    Regulation = case_when(
      adj.P.Val < 0.05 & logFC >  0.585 ~ "Up",
      adj.P.Val < 0.05 & logFC < -0.585 ~ "Down",
      TRUE                              ~ "NotSig"
    )
  )

write.csv(deg_table,
          file.path(DIR_RESULTS, "DEG_full_table.csv"),
          row.names = FALSE)

deg_sig <- deg_table |> filter(Regulation %in% c("Up", "Down"))
write.csv(deg_sig,
          file.path(DIR_RESULTS, "DEG_significant_FDR05_logFC0.585.csv"),
          row.names = FALSE)

cat(sprintf("Total significant DEGs: %d (Up: %d, Down: %d)\n",
            nrow(deg_sig),
            sum(deg_sig$Regulation == "Up"),
            sum(deg_sig$Regulation == "Down")))

# Volcano plot
p_volcano <- ggplot(deg_table,
                    aes(logFC, -log10(P.Value), color = Regulation)) +
  geom_point(size = 0.7, alpha = 0.7) +
  scale_color_manual(values = c(Up = "#E64B35", Down = "#3182BD",
                                NotSig = "#BFBFBF")) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano plot: Stroke vs Sham (merged validation cohort)",
       x = "log2 fold change", y = "-log10 P value") +
  theme_minimal(base_size = 12)

ggsave(file.path(DIR_FIG, "Fig3A_volcano.pdf"),
       p_volcano, width = 7, height = 5.5)

cat("DEG analysis complete. Outputs in results/ and figures/.\n")
