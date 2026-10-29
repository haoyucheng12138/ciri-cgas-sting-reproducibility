# Builds reference/GSE227317_Ensembl_transcript_map.tsv (rat transcript -> gene symbol) for the 7 genes analysed in
# 09_rat_roc.R by querying Ensembl BioMart (dataset rnorvegicus_gene_ensembl). The bundled file was retrieved on
# 2026-10-29 from Ensembl release 116. Re-running this script overwrites it with the current Ensembl release.
# Run from the repository root.
suppressMessages(library(httr))
genes <- c("Tap1", "Tap2", "Cxcl10", "Psmb8", "Stat1", "Psmb10", "Trpc6")
q <- paste0('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1">',
            '<Dataset name="rnorvegicus_gene_ensembl" interface="default"><Filter name="external_gene_name" value="', paste(genes, collapse = ","), '"/>',
            '<Attribute name="ensembl_transcript_id"/><Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/></Dataset></Query>')
txt <- NULL
for (h in c("https://www.ensembl.org", "https://useast.ensembl.org", "https://asia.ensembl.org")) {
  r <- tryCatch(GET(paste0(h, "/biomart/martservice"), query = list(query = q), timeout(180)), error = function(e) NULL)
  if (!is.null(r) && status_code(r) == 200 && grepl("^Transcript stable ID", content(r, "text", encoding = "UTF-8"))) { txt <- content(r, "text", encoding = "UTF-8"); break }
}
if (is.null(txt)) stop("Ensembl BioMart is unavailable; keep the bundled reference/GSE227317_Ensembl_transcript_map.tsv")
writeLines(txt, "reference/GSE227317_Ensembl_transcript_map.tsv")
rel <- tryCatch(content(GET("https://rest.ensembl.org/info/data/?content-type=application/json", timeout(60)))$releases[[1]], error = function(e) NA)
writeLines(paste("Ensembl release:", rel, "| BioMart rnorvegicus_gene_ensembl | retrieved", format(Sys.Date())), "reference/GSE227317_Ensembl_release.txt")
