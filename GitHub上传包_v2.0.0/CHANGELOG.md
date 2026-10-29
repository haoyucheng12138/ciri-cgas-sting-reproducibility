# Changelog

## v2.0.0 (2026-10-29) - revision of manuscript HELIYON-D-26-05075

This release replaces the previous code (v1). The previous repository described an independent glmnet
re-implementation of analyses originally run on a commercial analysis platform; its results did not correspond
exactly to the manuscript figures. v2 rebuilds every analysis from the raw GEO files, and every number in the revised
manuscript is written by the code to `results/number_ledger.csv`.

### Problems found in the original analysis and corrected in v2

1. **GSE131193 sample labels.** The columns of the GEO supplementary count file are grouped by condition, whereas the
   original analysis attached GSM accessions to them in accession order. 22 of 24 samples carried the wrong label; the
   original "Sham vs Stroke" comparison was in fact day 1 vs day 7. Columns are now assigned by name from the GEO
   sample titles (`data/processed/sample_column_map.csv`; audit in `reference/GSE131193_LKFlow_label_vs_true_content.csv`).
   All WGCNA, ssGSEA and STING-score results for this cohort were recomputed.
2. **GSE131193 tissue.** The samples are brain microvessel fragments (GEO), not whole hemisphere.
3. **Primary comparison in GSE131193** is now ipsilateral tMCAO vs ipsilateral sham (6 vs 6); the original grouping
   mixed ipsilateral and contralateral hemispheres.
4. **PPI input.** Gene symbols passed to STRING had been split at the letters "r", "n" and "t" (e.g. "Irgm1" -> "I", "gm1").
   The PPI is rebuilt from the 1,421 DEGs.
5. **GO/KEGG input** had been 276 genes (DEG intersected with a WGCNA list built on mislabelled traits), not the 1,421 DEGs
   stated in the Methods. Enrichment is rebuilt from the 1,421 DEGs.
6. **Three-set intersection inputs** contained truncated symbols ("rb1", "f10", numeric fragments). Rebuilt from clean sets.
7. **WGCNA** had used 17 of 24 samples (cut height 100) and a network type/power that differ from the Methods. v2 keeps all
   24 samples (no sample with Z.k < -2.5) and uses an unsigned network at the lowest power with scale-free R^2 >= 0.85 (12).
8. **Immune gene sets.** The original 27 signatures came from an undocumented platform list. v2 uses the published
   Charoentong et al. 2017 Table S6 (28 cell types), mapped to mouse orthologs.
9. **Elastic Net coefficient tables (former Table S3/S4)** did not match the coefficients of the fitted models. v2 reports
   the coefficients of the fitted glmnet models.
10. **Preprocessing description.** The merged cohort is log2(FPKM + 1) (not TPM). Gene filters, duplicate handling,
    Excel-damaged symbols ("1-Mar" -> "March1") and the Cgas/Mb21d1 synonym are now explicit in `R/02_preprocess.R`.
11. **GSE227317** labels are taken from GEO (T01-T05 = MCAO); Psmb10 is recomputed from the transcript-level FPKM.

### Analyses added in response to reviewers

- Batch-correction sensitivity: ComBat without the group covariate, batch as a limma covariate, each dataset alone (R2.5).
- STING score with interferon-stimulated genes removed; immune gene sets without STING-score genes (R2.5).
- Elastic Net on the full transcriptome (R1.1).
- WGCNA with several trait definitions (R1.6).
