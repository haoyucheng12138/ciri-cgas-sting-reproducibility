# =============================================================================
# 04_enrichment_ppi.R  -  GO / KEGG enrichment and STRING PPI of the 1,421 merged-cohort DEGs.
#  clusterProfiler enrichGO / enrichKEGG (organism "mmu", BH, p < 0.05, q < 0.2).
#  STRING v11.5 REST API, Mus musculus (10090), combined score >= 0.7; node degree / betweenness (igraph).
#  KEGG and STRING are queried online; the exact downloaded responses are saved for reproducibility.
# =============================================================================
source("R/00_config.R")
suppressMessages({ library(clusterProfiler); library(org.Mm.eg.db); library(httr); library(igraph) })
out <- res_dir("04_enrichment_ppi")
deg <- read.csv(file.path(res_dir("03_combat_limma"), "TableS6_significant_DEGs.csv"), row.names = 1)
genes <- rownames(deg)
m <- suppressMessages(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
cat(sprintf("DEGs %d | mapped to Entrez %d\n", length(genes), length(unique(m$SYMBOL))))
for (ont in c("BP", "CC", "MF")) {
  eg <- enrichGO(m$ENTREZID, OrgDb = org.Mm.eg.db, ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  df <- as.data.frame(eg); write.csv(df, file.path(out, paste0("GO_", ont, ".csv")), row.names = FALSE)
  cat(sprintf("GO-%s: %d terms | top 5: %s\n", ont, nrow(df), paste(head(df$Description, 5), collapse = " | ")))
}
kk <- setReadable(enrichKEGG(m$ENTREZID, organism = "mmu", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2), OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
kd <- as.data.frame(kk); write.csv(kd, file.path(out, "KEGG.csv"), row.names = FALSE)
cat(sprintf("KEGG: %d pathways\n", nrow(kd)))
print(kd[kd$ID %in% c("mmu04612", "mmu04623", "mmu04621", "mmu04620"), c("ID", "Description", "Count", "p.adjust")], row.names = FALSE)

r <- POST("https://version-11-5.string-db.org/api/tsv/network",
          body = list(identifiers = paste(genes, collapse = "%0d"), species = 10090, required_score = 700, caller_identity = "ciri_revision"),
          encode = "form", timeout(300))
stop_for_status(r)
txt <- content(r, "text", encoding = "UTF-8"); writeLines(txt, file.path(out, "STRING_v11.5_raw_response.tsv"))
net <- read.delim(text = txt)[, c("preferredName_A", "preferredName_B", "score")]
net <- net[!duplicated(apply(net[, 1:2], 1, function(x) paste(sort(x), collapse = "|"))), ]
g <- graph_from_data_frame(net[, 1:2], directed = FALSE)
nodes <- data.frame(Gene = V(g)$name, Degree = degree(g), Betweenness = betweenness(g, normalized = TRUE), Closeness = closeness(g, normalized = TRUE))
nodes <- nodes[order(-nodes$Degree, -nodes$Betweenness), ]
write.csv(net, file.path(out, "PPI_edges_score700.csv"), row.names = FALSE); write.csv(nodes, file.path(out, "PPI_nodes_score700.csv"), row.names = FALSE)
cat(sprintf("PPI: %d nodes, %d edges | top 20: %s\n", vcount(g), ecount(g), paste(sprintf("%s(%d)", head(nodes$Gene, 20), head(nodes$Degree, 20)), collapse = " ")))
save_session("04_enrichment_ppi")
