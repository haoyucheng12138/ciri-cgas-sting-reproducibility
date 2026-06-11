# =============================================================================
# run_all.R
# Convenience wrapper: execute the full reproducibility pipeline in order.
# Run individual scripts manually if you want to inspect intermediate output.
# =============================================================================

scripts <- c(
  "00_setup.R",
  "01_download_geo.R",
  "02_preprocessing.R",
  "03_combat_batch_correction.R",
  "04_limma_deg.R",
  "05_wgcna.R",
  "06_ssgsea_sting_score.R",
  "07_elasticnet_modeling.R",
  "08_roc_bootstrap.R"
)

for (s in scripts) {
  cat(sprintf("\n========== Running %s ==========\n", s))
  source(s, echo = FALSE)
}

cat("\nFull pipeline complete. See results/ and figures/ for outputs.\n")
