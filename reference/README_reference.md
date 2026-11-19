# reference/

Bundled inputs that are not GEO downloads. See the main README for what each file is.

`GSE131193_LKFlow_label_vs_true_content.csv` and `GSE131193_column_to_GSM_map.csv` document the sample-labelling
error of the original (platform-based) analysis: for each GSM accession, the GEO column whose data had been placed
under it, compared with the GEO column that actually belongs to it. They are not used by the pipeline. `GSE131193_LKFlow_label_audit.csv` (written by the optional `R/02_preprocess_check.R`)
lists, for each GSM, the label the original analysis attached to the same data; it is read by `R/13_supplementary_tables.R`
for Supplementary Table S2.
