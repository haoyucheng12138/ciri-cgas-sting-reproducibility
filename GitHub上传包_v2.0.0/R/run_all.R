# =============================================================================
# run_all.R  -  run the complete pipeline from the raw GEO files.
# From the repository root:   Rscript R/run_all.R
# Approximate run time on a desktop PC: 1-1.5 h (WGCNA network construction takes ~35 min single-threaded).
# Steps 04 (KEGG, STRING) and 09 (Ensembl, only if the bundled map is missing) query online services.
# =============================================================================
steps <- c("R/01_download_geo.R", "R/02_preprocess.R", "R/03_combat_limma.R", "R/04_enrichment_ppi.R", "R/05_wgcna.R",
           "R/06_ssgsea_sting.R", "R/07_intersection.R", "R/08_elastic_net.R", "R/09_rat_roc.R", "R/10_number_ledger.R")
for (s in steps) {
  cat(sprintf("\n==================== %s  (%s) ====================\n", s, format(Sys.time(), "%H:%M:%S")))
  rc <- system2(file.path(R.home("bin"), "Rscript"), s)
  if (rc != 0) stop("step failed: ", s)
}
cat("\nPipeline complete. All manuscript numbers: results/number_ledger.csv\n")
