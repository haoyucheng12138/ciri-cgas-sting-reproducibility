# =============================================================================
# build_Charoentong28_mouse_gmt.R  -  how reference/Charoentong28_mouse.gmt was made (run from the repository root).
# Source: Charoentong P, et al. Pan-cancer immunogenomic analyses reveal genotype-immunophenotype relationships and
#   predictors of response to checkpoint blockade. Cell Rep 2017;18:248-262. doi:10.1016/j.celrep.2016.12.019
#   Supplementary Table S6 (file 1-s2.0-S2211124716317090-mmc3.xlsx), bundled as reference/Charoentong2017_TableS6_original.xlsx
# Human symbols -> mouse orthologs with babelgene (HGNC HCOP consensus; min_support = 3; top = TRUE).
# Result: 28 immune cell types, 782 human genes, 709 mapped (coverage per set: reference/Charoentong28_mouse_coverage.csv).
# =============================================================================
for (p in c("readxl", "babelgene")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
s6 <- readxl::read_excel("reference/Charoentong2017_TableS6_original.xlsx", skip = 2)
names(s6)[1:3] <- c("gene", "cell", "immunity")
s6 <- s6[!is.na(s6$gene) & !is.na(s6$cell), ]; s6$cell <- gsub("\\s+", " ", trimws(s6$cell))
orth <- babelgene::orthologs(genes = unique(s6$gene), species = "mouse", human = TRUE, min_support = 3, top = TRUE)
map <- unique(orth[, c("human_symbol", "symbol")])
gs <- lapply(split(s6$gene, s6$cell), function(h) sort(unique(map$symbol[map$human_symbol %in% h])))
con <- file("reference/Charoentong28_mouse.gmt", "w", encoding = "UTF-8")
for (k in names(gs)) writeLines(paste(c(gsub(" ", "_", k), "Charoentong2017_TableS6_mouse_ortholog", gs[[k]]), collapse = "\t"), con)
close(con)
cat("cell types:", length(gs), "| human genes:", length(unique(s6$gene)), "| mapped:", length(unique(map$human_symbol)),
    "| babelgene", as.character(packageVersion("babelgene")), "\n")
