# =============================================================================
# 04_enrichment_ppi.R  -  GO / KEGG enrichment and STRING PPI of the 1,421 merged-cohort DEGs.
#  Symbols -> Entrez: official symbols, then unambiguous org.Mm.eg.db aliases (former names such as Tmem173 -> Sting1).
#  clusterProfiler enrichGO / enrichKEGG (organism "mmu", BH, p < 0.05, q < 0.2).
#  STRING v11.5 REST API, Mus musculus (10090), combined score >= 0.7; node degree / betweenness (igraph).
#  KEGG and STRING are queried online; the downloaded STRING responses are saved for reproducibility.
# =============================================================================
source("R/00_config.R")
suppressMessages({ library(clusterProfiler); library(org.Mm.eg.db); library(httr); library(igraph) })
set.seed(SEED)
out <- res_dir("04_enrichment_ppi")
deg <- read.csv(file.path(res_dir("03_combat_limma"), "TableS6_significant_DEGs.csv"), row.names = 1)
genes <- rownames(deg)

m1 <- suppressMessages(suppressWarnings(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)))
left <- setdiff(genes, m1$SYMBOL)
al <- suppressMessages(suppressWarnings(tryCatch(bitr(left, fromType = "ALIAS", toType = "ENTREZID", OrgDb = org.Mm.eg.db),
                                                 error = function(e) data.frame(ALIAS = character(), ENTREZID = character()))))
al <- al[!al$ALIAS %in% al$ALIAS[duplicated(al$ALIAS)] & !al$ENTREZID %in% m1$ENTREZID, ]
map <- rbind(data.frame(input = m1$SYMBOL, ENTREZID = m1$ENTREZID, via = "symbol"), data.frame(input = al$ALIAS, ENTREZID = al$ENTREZID, via = "alias"))
map$official <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, map$ENTREZID, "SYMBOL", "ENTREZID"))
write.csv(map, file.path(out, "DEG_to_Entrez_map.csv"), row.names = FALSE)
write.csv(data.frame(unmapped = setdiff(genes, map$input)), file.path(out, "DEG_unmapped.csv"), row.names = FALSE)
cat(sprintf("DEGs %d | mapped by official symbol %d | by unique alias %d | unmapped %d\n", length(genes), nrow(m1), nrow(al), length(setdiff(genes, map$input))))
cat("alias-mapped:", paste(sprintf("%s->%s", map$input[map$via == "alias"], map$official[map$via == "alias"]), collapse = ", "), "\n")

## Enrichment is run once with no p/q cut-off so that every tested term is saved (*_all_tested.csv) and the BH FDR can
## be recomputed by readers; the significant terms (BH FDR < 0.05 and q < 0.2) are then taken from that table.
sig <- function(df) df[df$p.adjust < 0.05 & df$qvalue < 0.2, ]
for (ont in c("BP", "CC", "MF")) {
  eg <- enrichGO(map$ENTREZID, OrgDb = org.Mm.eg.db, ont = ont, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  all <- as.data.frame(eg); write.csv(all, file.path(out, paste0("GO_", ont, "_all_tested.csv")), row.names = FALSE)
  df <- sig(all); write.csv(df, file.path(out, paste0("GO_", ont, ".csv")), row.names = FALSE)
  cat(sprintf("GO-%s: %d tested, %d significant | top 5: %s\n", ont, nrow(all), nrow(df), paste(head(df$Description, 5), collapse = " | ")))
}
kk <- setReadable(enrichKEGG(map$ENTREZID, organism = "mmu", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1), OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
ka <- as.data.frame(kk); write.csv(ka, file.path(out, "KEGG_all_tested.csv"), row.names = FALSE)
kd <- sig(ka); write.csv(kd, file.path(out, "KEGG.csv"), row.names = FALSE)
cat(sprintf("KEGG: %d tested\n", nrow(ka)))
cat(sprintf("KEGG: %d pathways\n", nrow(kd)))
print(kd[kd$ID %in% c("mmu04612", "mmu04623", "mmu04621", "mmu04620"), c("ID", "Description", "Count", "p.adjust")], row.names = FALSE)

## ---- STRING PPI ----
string_post <- function(endpoint, extra) {
  r <- POST(paste0("https://version-11-5.string-db.org/api/tsv/", endpoint),
            body = c(list(identifiers = paste(genes, collapse = "%0d"), species = 10090, caller_identity = "ciri_revision"), extra),
            encode = "form", timeout(300))
  stop_for_status(r); content(r, "text", encoding = "UTF-8")
}
txt <- string_post("network", list(required_score = 700)); writeLines(txt, file.path(out, "STRING_v11.5_raw_response.tsv"))
idt <- string_post("get_string_ids", list(echo_query = 1, limit = 1)); writeLines(idt, file.path(out, "STRING_v11.5_id_resolution.tsv"))
raw <- read.delim(text = txt); ids <- read.delim(text = idt)
## STRING labels some proteins with older synonyms (e.g. "Irg1" for the submitted DEG "Acod1"); nodes are relabelled by
## STRING protein ID with the DEG symbol that was submitted.
sym_of <- setNames(ids$queryItem, ids$stringId); sym_of <- sym_of[!duplicated(names(sym_of))]
lab <- function(sid, nm) ifelse(sid %in% names(sym_of), unname(sym_of[sid]), nm)
net <- data.frame(gene_A = lab(raw$stringId_A, raw$preferredName_A), gene_B = lab(raw$stringId_B, raw$preferredName_B), score = raw$score)
net <- net[!duplicated(apply(net[, 1:2], 1, function(x) paste(sort(x), collapse = "|"))), ]
cat(sprintf("STRING: %d of %d DEGs resolved | proteins relabelled to the submitted symbol: %d | node labels that are not DEG symbols: %d\n",
            length(unique(ids$queryItem)), length(genes), length(setdiff(unique(c(raw$preferredName_A, raw$preferredName_B)), genes)),
            length(setdiff(unique(c(net$gene_A, net$gene_B)), genes))))
g <- graph_from_data_frame(net[, 1:2], directed = FALSE)
nodes <- data.frame(Gene = V(g)$name, Degree = degree(g), Betweenness = betweenness(g, normalized = TRUE), Closeness = closeness(g, normalized = TRUE))
nodes <- nodes[order(-nodes$Degree, -nodes$Betweenness), ]
write.csv(net, file.path(out, "PPI_edges_score700.csv"), row.names = FALSE); write.csv(nodes, file.path(out, "PPI_nodes_score700.csv"), row.names = FALSE)
cat(sprintf("PPI: %d nodes, %d edges | top 20: %s\n", vcount(g), ecount(g), paste(sprintf("%s(%d)", head(nodes$Gene, 20), head(nodes$Degree, 20)), collapse = " ")))
save_session("04_enrichment_ppi")
