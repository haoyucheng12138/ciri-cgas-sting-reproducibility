# Reproducibility Code for: Systems Biology and Machine Learning Suggest Coordinated Activation of the cGAS-STING Axis and Antigen Presentation in Cerebral Ischemia-Reperfusion Injury

This repository provides the analytical code required to independently reproduce the bioinformatics and machine learning analyses reported in the manuscript. All analyses were performed on publicly available transcriptomic datasets retrieved from the NCBI Gene Expression Omnibus (GEO).

## Data sources (all public, no access restriction)

| Accession | Species | Role |
|-----------|---------|------|
| GSE131193 | Mus musculus | Discovery cohort (WGCNA, ssGSEA) |
| GSE131712 | Mus musculus | Merged validation cohort (with GSE202659) |
| GSE202659 | Mus musculus | Merged validation cohort (with GSE131712) |
| GSE23163  | Mus musculus | Machine learning hold-out tasks |
| GSE227317 | Rattus norvegicus | Cross-species external evaluation |

## Software environment

- R version: 4.3.0 or higher
- Required CRAN packages: `glmnet`, `pROC`, `WGCNA`, `umap`, `ggplot2`, `dplyr`, `tibble`, `readr`, `Matrix`
- Required Bioconductor packages: `GEOquery`, `limma`, `sva`, `GSVA`, `clusterProfiler`, `org.Mm.eg.db`
- The exact package versions used in the original analysis are recorded in `sessionInfo.txt`.

## Repository layout

```
ciri-repro/
├── README.md                       # this file
├── 00_setup.R                      # package installation and global config
├── 01_download_geo.R               # download all 5 GEO datasets
├── 02_preprocessing.R              # log2(TPM+1), missing imputation, ortholog mapping
├── 03_combat_batch_correction.R    # ComBat + diagnostic plots
├── 04_limma_deg.R                  # limma differential expression
├── 05_wgcna.R                      # WGCNA module identification
├── 06_ssgsea_sting_score.R         # ssGSEA + STING composite score
├── 07_elasticnet_modeling.R        # Elastic Net (alpha=0.4, 10-fold CV, seed=12345)
├── 08_roc_bootstrap.R              # ROC curves and bootstrap 95% CIs
├── run_all.R                       # one-click full pipeline
├── data/                           # downloaded GEO data (created by 01)
├── results/                        # tabular outputs
├── figures/                        # diagnostic plots
└── sessionInfo.txt                 # frozen package versions
```

## How to run

```bash
# Step 1: install dependencies (one-time)
Rscript 00_setup.R

# Step 2: full pipeline
Rscript run_all.R

# OR run scripts individually in order
Rscript 01_download_geo.R
Rscript 02_preprocessing.R
Rscript 03_combat_batch_correction.R
Rscript 04_limma_deg.R
Rscript 05_wgcna.R
Rscript 06_ssgsea_sting_score.R
Rscript 07_elasticnet_modeling.R
Rscript 08_roc_bootstrap.R
```

## Reproducibility notes

1. All random operations are seeded with `set.seed(12345)`.
2. Elastic Net hyperparameters (alpha = 0.4, 10-fold cross-validation, lambda.min) are explicitly set in `07_elasticnet_modeling.R`.
3. The candidate feature space (39 genes) used in the machine learning analyses is hard-coded in `07_elasticnet_modeling.R` and corresponds to the curated set of cGAS-STING/IFN axis genes, MHC-I antigen processing components, PPI network hub genes, and Trpc6 as described in the Methods section.
4. Downstream numerical results may differ minimally across R/package versions; the `sessionInfo.txt` file records the exact environment used in the original analysis.
5. The original analyses reported in the manuscript were executed using the LKFlow analytical platform; this repository provides an independent open-source reimplementation using the `glmnet` R package with the parameters reported in the Methods section.

## Citation

If you use this code, please cite the corresponding manuscript (in preparation).

## Contact

For technical questions, please open a GitHub issue or contact the corresponding author.
