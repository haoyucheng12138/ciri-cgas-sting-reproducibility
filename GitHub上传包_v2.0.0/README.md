# Reproducibility code: multi-cohort transcriptomics of cGAS-STING-associated signalling and MHC class I antigen-processing genes in cerebral ischemia-reperfusion injury

Code for manuscript HELIYON-D-26-05075 (revised version). The pipeline starts from the public GEO files and
reproduces every number reported in the revised manuscript; the numbers are collected automatically in
`results/number_ledger.csv`. See `CHANGELOG.md` for the corrections made relative to the original submission.

## Data (all public, NCBI GEO)

| Accession | Species | Tissue | Samples used | Role |
|---|---|---|---|---|
| GSE131712 | Mus musculus | ipsilateral cortex | 3 Sham + 3 tMCAO 24 h (control shRNA) | merged validation cohort |
| GSE202659 | Mus musculus | hemisphere | 3 Sham + 3 tMCAO | merged validation cohort |
| GSE131193 | Mus musculus | brain microvessel fragments | 24 (Sham/tMCAO x 1 d/7 d x ipsi/contralateral) | discovery cohort (WGCNA, ssGSEA) |
| GSE23163 | Mus musculus | cortex, striatum (WT, Gpx1-/-) | 32 (Sham and I/R 24 h) | Elastic Net hold-out tasks |
| GSE227317 | Rattus norvegicus | cortex | 5 MCAO 2 h/R 24 h + 5 Sham | exploratory cross-species check |

`R/01_download_geo.R` downloads the files from the NCBI GEO FTP and writes their MD5 sums to `data/raw/MD5SUMS.txt`.

## Software

R 4.6.1 (Windows x64). Package versions used for the reported results: limma 3.68.5, sva 3.60.0, WGCNA 1.74,
GSVA 2.6.15, clusterProfiler 4.20.3, org.Mm.eg.db 3.23.0, glmnet 5.1, pROC 1.19.1, igraph 2.3.4, httr 1.4.9.
Install with `Rscript R/00_install_packages.R`. Each step writes its `sessionInfo.txt` to `results/<step>/`.

Online resources queried: KEGG REST (release of 2026-10-28, via clusterProfiler), STRING v11.5 API (the raw response is
saved to `results/04_enrichment_ppi/`), Ensembl BioMart release 116 (rat transcript map, bundled in `reference/`).

## Run

From the repository root:

```
Rscript R/00_install_packages.R   # once
Rscript R/run_all.R               # about 1-1.5 h; WGCNA takes ~35 min
```

## Pipeline

| Script | Content |
|---|---|
| `01_download_geo.R` | download GEO files, write MD5 sums |
| `02_preprocess.R` | build matrices; assign sample columns to GSM accessions by name (`data/processed/sample_column_map.csv`) |
| `03_combat_limma.R` | merged cohort: ComBat + limma (FDR < 0.05, \|log2FC\| > 0.585); batch-correction sensitivity analyses; STING score |
| `04_enrichment_ppi.R` | GO/KEGG (clusterProfiler) and STRING v11.5 PPI (score >= 0.7) of the DEGs |
| `05_wgcna.R` | WGCNA on GSE131193 (24 samples, unsigned, power 12, minModuleSize 30, deepSplit 2, mergeCutHeight 0.3) |
| `06_ssgsea_sting.R` | ssGSEA (Charoentong 2017, 28 immune cell types), STING score, group tests, correlations, overlap sensitivity |
| `07_intersection.R` | DEG x WGCNA core x curated cGAS-STING genes |
| `08_elastic_net.R` | Elastic Net (alpha 0.4, 10-fold CV, lambda.min, seed 12345): original 39-gene panel, full transcriptome, rebuilt panel |
| `09_rat_roc.R` | GSE227317 ROC with 2,000 stratified bootstrap CIs |
| `10_number_ledger.R` | collects all reported numbers into `results/number_ledger.csv` |
| `02_preprocess_check.R` | optional, authors only: compares the matrices with those of the original platform analysis |

## Preprocessing rules (`R/02_preprocess.R`)

- **GSE131712, GSE202659**: gene-level FPKM; duplicate symbols collapsed by median; log2(FPKM + 1); genes with mean log2 > 0.1 kept; merged cohort = genes present in both (20,078 genes).
- **GSE131193**: normalised counts from the submitters; columns assigned to GSM accessions by name (the file's column order differs from accession order); duplicates by median; log2(x + 1); genes with mean log2 > 1 kept (17,087 genes).
- **GSE23163**: GeneSpring median-normalised ratios from the series matrix; log2(x + 1) per probe; probes collapsed to GPL6885 gene symbols by median; genes with mean >= 1 kept (7,693 genes).
- Gene symbols: Excel-converted symbols restored ("1-Mar" -> "March1"); Ensembl gene ENSMUSG00000032344 is named "Mb21d1" in all datasets (GSE202659 uses its current symbol "Cgas").

## Reference files (`reference/`)

| File | Content |
|---|---|
| `TableS2_merged_cohort_samples.csv` | the 12 samples of the merged cohort |
| `Charoentong2017_TableS6_original.xlsx` | original immune metagene table (Charoentong et al., Cell Rep 2017) |
| `Charoentong28_mouse.gmt`, `..._coverage.csv`, `build_Charoentong28_mouse_gmt.R` | mouse ortholog gene sets and how they were built |
| `cGAS_STING_curated_genes.csv` | curated cGAS-STING pathway genes (29 symbols) |
| `original_39gene_candidate_panel.csv` | Elastic Net candidate panel of the original submission |
| `GSE227317_Ensembl_transcript_map.tsv`, `..._release.txt`, `build_GSE227317_transcript_map.R` | rat transcript-to-gene map (Ensembl 116) |
| `GSE131193_column_to_GSM_map.csv`, `GSE131193_LKFlow_label_vs_true_content.csv` | audit of the GSE131193 labelling error in the original analysis |

## Notes on interpretation

- The Elastic Net hold-out tasks use a single study (GSE23163); the AUCs describe within-study generalisation and are not external validation.
- GSE227317 has 5 vs 5 samples; an AUC of 1 has a saturated confidence interval and is reported as a directional check only.
- Results are transcript-level associations; they do not measure protein levels, MHC-I surface expression or antigen presentation.

## Citation

Please cite the article (details will be added on publication) and this repository's Zenodo record.
