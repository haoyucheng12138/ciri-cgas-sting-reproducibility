# Changelog

## v2.2.1 (2026-11) - supplementary table numbering

No analysis was changed: all files under `results/` except the supplementary table file, and all figures, are identical to v2.2.0.

- Supplementary tables renumbered in order of first citation in the revised manuscript: former Table S5 (GO/KEGG enrichment
  and PPI) is now Table S7, former Table S6 (WGCNA) is now Table S5, and former Table S7 (immune signatures and IFN/STING
  score) is now Table S6. Tables S1-S4 and S8-S11 are unchanged. Only sheet names, table titles and sheet order of
  `Supplementary_Tables_S1-S11.xlsx` changed; cell contents are identical.

## v2.2.0 (2026-11) - revision with experiments; documentation corrections

No analysis was changed: all transcriptomic results, `results/number_ledger.csv` and Figures 2-6, S1-S3 are identical in content to v2.1.0.

- Figure 1: the bottom box now refers to the experiments of the revision (Supplementary Figs. S4-S6) and lists the measured
  readouts; "signalling" changed to "signaling". Figure 1 and the other main figures are written with a white background
  (Figure 1 previously had a transparent background).
- Supplementary tables: one file `Supplementary_Tables_S1-S11.xlsx`; new Table S11 (qPCR primers, from
  `reference/qPCR_primers_TableS11.csv`); Table S6A (Table S5A since v2.2.1) gains the Benjamini-Hochberg FDR of the primary-trait correlation
  (the values reported in the manuscript, previously only in `results/05_wgcna/module_trait_primary_sorted.csv`);
  "merged validation cohort" and "discovery cohort" replaced by "merged cohort" and "GSE131193 (brain microvessel fragments)".
- Wording: `results/number_ledger.csv` lists the main reported numbers, not every number; all other values are in the
  result files under `results/`.
- Correction to the v2.0.0 notes below: 22 of the 24 GSE131193 columns were attached to the wrong GSM accession; 20 of the
  24 samples thereby received an incorrect condition (injury, time point or hemisphere); 4 kept a correct condition.

## v2.1.0 (2026-11) - figures, supplementary tables, documentation

- Added scripts 04_enrichment_check, 08b (seed sensitivity, convergence diagnostics), 11-14 (all figures and the
  supplementary table file).
- The composite score is labelled "IFN/STING score" in figures and tables (reviewer 1, comment 5).
- `results/number_ledger.csv` reports PPI degrees instead of degree ranks (ranks were ambiguous because of ties) and
  the Charoentong gene-set coverage (709 mouse genes, 618 in the GSE131193 matrix).
- Corrected comments: B2m is on the GSE23163 array and is removed by the expression filter; the Excel-symbol
  safeguard in 02_preprocess.R changes nothing in the current GEO files.
- Added `reference/dataset_primary_publications.csv`.

## v2.0.0 (2026-10-29) - revision of manuscript HELIYON-D-26-05075

This release replaces the previous code (v1). The previous repository described an independent glmnet
re-implementation of analyses originally run on a commercial analysis platform; its results did not correspond
exactly to the manuscript figures. v2 rebuilds every analysis from the raw GEO files, and every number in the revised
manuscript is written by the code to `results/number_ledger.csv`.

### Problems found in the original analysis and corrected in v2

1. **GSE131193 sample labels.** The columns of the GEO supplementary count file are grouped by condition, whereas the
   original analysis attached GSM accessions to them in accession order. 22 of 24 columns carried the wrong accession (see v2.2.0); the
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
10. **Preprocessing description.** The merged cohort is log2(FPKM + 1) (not TPM). Gene filters, duplicate handling
    and the Cgas/Mb21d1 synonym are now explicit in `R/02_preprocess.R`.
11. **GSE227317** labels are taken from GEO (T01-T05 = MCAO); Psmb10 is recomputed from the transcript-level FPKM.

### Analyses added in response to reviewers

- Batch-correction sensitivity: ComBat without the group covariate, batch as a limma covariate, each dataset alone (R2.5).
- STING score with interferon-stimulated genes removed; immune gene sets without STING-score genes (R2.5).
- Elastic Net on the full transcriptome (R1.1).
- WGCNA with several trait definitions (R1.6).
