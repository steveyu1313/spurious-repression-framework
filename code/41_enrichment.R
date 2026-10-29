# ============================================================
# 41_enrichment.R
# 目标: GO(BP/CC/MF) + KEGG + Hallmark GSEA + 上游 TF/激酶
#       背景基因集统一使用本芯片全部可注释基因(21355), 避免默认全基因组导致的 P 值虚高
# 运行: cd 项目根 && Rscript 08_QC_Revision/scripts/41_enrichment.R
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
  library(msigdbr); library(fgsea); library(jsonlite)
})
tab_dir <- "./08_QC_Revision/tables"
tab_dir <- "./08_QC_Revision/tables"
sig  <- fread(file.path(tab_dir, "T12_DEG_significant.csv"))
allg <- fread(file.path(tab_dir, "T6_GSE46602_DEG_rechecked.csv"))
allg <- allg[!is.na(Symbol) & Symbol != "" & !duplicated(Symbol)]

# ---------- ID 映射 ----------
map <- suppressWarnings(bitr(allg$Symbol, "SYMBOL", "ENTREZID", org.Hs.eg.db))
setDT(map); setnames(map, c("SYMBOL", "ENTREZID"), c("Symbol", "ENTREZID"), skip_absent = TRUE)
universe_ez <- unique(map$ENTREZID)

sig2 <- merge(sig, map, by = "Symbol")
sig2 <- sig2[!duplicated(ENTREZID)]
cat(sprintf("DEG %d -> 有 ENTREZ 注释 %d ; universe %d\n",
            nrow(sig), nrow(sig2), length(universe_ez)))

mk <- function(x, tag) {
  if (nrow(x) < 5) return(NULL)
  r <- try(enrichGO(gene = x$ENTREZID, universe = universe_ez, OrgDb = org.Hs.eg.db,
                    ont = tag, keyType = "ENTREZID", pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
           silent = TRUE)
  if (inherits(r, "try-error") || is.null(r) || nrow(as.data.frame(r)) == 0) return(NULL)
  as.data.table(as.data.frame(r))
}

sets <- list(All = sig2, Up = sig2[direction == "Up"], Down = sig2[direction == "Down"])
out <- list()
for (nm in names(sets)) for (ont in c("BP", "CC", "MF")) {
  d <- mk(sets[[nm]], ont)
  if (!is.null(d)) { d[, `:=`(gene_set = nm, ontology = ont)]; out[[paste(nm, ont)]] <- d }
}

kk <- try(enrichKEGG(gene = sig2$ENTREZID, universe = universe_ez, organism = "hsa",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
          silent = TRUE)
kegg <- NULL
if (!inherits(kk, "try-error") && !is.null(kk) && nrow(as.data.frame(kk)) > 0) {
  kegg <- as.data.table(as.data.frame(setReadable(kk, org.Hs.eg.db, "ENTREZID")))
  kegg[, gene_set := "All"]
}
if (!is.null(kegg)) out[["KEGG"]] <- kegg

# ---------- 汇总 GO ----------
go <- rbindlist(out[grep("^(All|Up|Down) (BP|CC|MF)$", names(out))],
                use.names = TRUE, fill = TRUE)
setorder(go, gene_set, ontology, p.adjust)
fwrite(go, file.path(tab_dir, "T19_GO_enrichment.csv"))
fwrite(kegg, file.path(tab_dir, "T20_KEGG_enrichment.csv"))
cat(sprintf("GO 条目 %d ; KEGG 通路 %d\n", nrow(go), if (is.null(kegg)) 0 else nrow(kegg)))

# ---------- GSEA: Hallmark + TFT(上游转录因子) ----------
rnk <- allg[!is.na(t) & !is.na(Symbol)]
rnk <- rnk[order(-abs(t))]
rnk <- rnk[!duplicated(Symbol)]
rk <- setNames(rnk$t, rnk$Symbol)
rk <- sort(rk, decreasing = TRUE)

run_gsea <- function(collection, subcol = NULL, tag) {
  hs <- if (is.null(subcol))
    msigdbr(species = "Homo sapiens", collection = collection)
  else msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcol)
  lst <- split(hs$gene_symbol, hs$gs_name)
  set.seed(20261020)
  r <- suppressWarnings(fgsea(lst, rk, minSize = 15, maxSize = 500, nPermSimple = 10000))
  r <- as.data.table(r)
  r <- r[!is.na(padj)]
  r[, leadingEdge := sapply(leadingEdge, function(x) paste(head(x, 30), collapse = ","))]
  setorder(r, padj, -NES)
  r[, collection := tag]
  r
}
hall <- try(run_gsea("H", NULL, "Hallmark"), silent = TRUE)
tft  <- try(run_gsea("C3", "TFT:TFT_LEGACY", "TFT"), silent = TRUE)
if (!inherits(hall, "try-error")) fwrite(hall, file.path(tab_dir, "T21_GSEA_hallmark.csv"))
if (!inherits(tft, "try-error"))  fwrite(tft,  file.path(tab_dir, "T22_GSEA_TFT.csv"))
cat(sprintf("Hallmark 显著(padj<0.05) %d ; TFT 显著 %d\n",
            if (inherits(hall, "try-error")) 0 else nrow(hall[padj < 0.05]),
            if (inherits(tft, "try-error")) 0 else nrow(tft[padj < 0.05])))

# ---------- ChEA3 (上游 TF) / KEA3 (上游激酶) ----------
post_json <- function(url, genes, name) {
  body <- toJSON(list(gene_set = genes, query_name = name), auto_unbox = TRUE)
  tmp <- tempfile(fileext = ".json")
  st <- system2("curl", c("-s", "-m", "180", "-X", "POST", "-H",
                          shQuote("Content-Type: application/json"),
                          "-d", shQuote(body), "-o", shQuote(tmp), shQuote(url)),
                stdout = FALSE, stderr = FALSE)
  if (!file.exists(tmp) || file.info(tmp)$size == 0) return(NULL)
  j <- try(fromJSON(tmp, simplifyDataFrame = TRUE), silent = TRUE)
  if (inherits(j, "try-error")) return(NULL)
  k <- grep("meanRank", names(j), value = TRUE)[1]
  if (is.na(k)) return(NULL)
  as.data.table(j[[k]])
}

hub <- fread(file.path(tab_dir, "T17_hub_genes.csv"))$hub
chea_all <- post_json("https://maayanlab.cloud/chea3/api/enrich/", sig$Symbol, "DEG_all")
chea_hub <- post_json("https://maayanlab.cloud/chea3/api/enrich/", hub, "Hub_top20")
kea_all  <- post_json("https://maayanlab.cloud/kea3/api/enrich/",  sig$Symbol, "DEG_all")
if (!is.null(chea_all)) fwrite(chea_all, file.path(tab_dir, "T23_upstream_TF_ChEA3.csv"))
if (!is.null(chea_hub)) fwrite(chea_hub, file.path(tab_dir, "T24_upstream_TF_hub.csv"))
if (!is.null(kea_all))  fwrite(kea_all,  file.path(tab_dir, "T25_upstream_kinase_KEA3.csv"))
cat("ChEA3/TF:", if (is.null(chea_all)) "fail" else nrow(chea_all),
    "| KEA3/kinase:", if (is.null(kea_all)) "fail" else nrow(kea_all), "\n")

# ---------- 下游: hub 一级邻居的功能富集 ----------
nb <- fread(file.path(tab_dir, "T18_hub_first_neighbours.csv"))$downstream_neighbour
nb2 <- merge(data.table(Symbol = nb), map, by = "Symbol")
go_nb <- mk(nb2, "BP")
if (!is.null(go_nb)) {
  go_nb[, gene_set := "hub_neighbour"]
  fwrite(go_nb, file.path(tab_dir, "T26_GO_hub_neighbours.csv"))
}
kk_nb <- try(enrichKEGG(nb2$ENTREZID, universe = universe_ez, organism = "hsa",
                        pvalueCutoff = 0.05, qvalueCutoff = 0.2), silent = TRUE)
if (!inherits(kk_nb, "try-error") && !is.null(kk_nb) && nrow(as.data.frame(kk_nb)) > 0) {
  kn <- as.data.table(as.data.frame(setReadable(kk_nb, org.Hs.eg.db, "ENTREZID")))
  fwrite(kn, file.path(tab_dir, "T27_KEGG_hub_neighbours.csv"))
  cat("hub 邻居 KEGG 通路:", nrow(kn), "\n")
}

cat("\n=== 显著 Hallmark (top8) ===\n")
if (!inherits(hall, "try-error"))
  print(hall[padj < 0.05][1:8, .(pathway, NES = round(NES, 2), padj = signif(padj, 2))])
cat("\n=== 显著 TFT 上游 TF (top12) ===\n")
if (!inherits(tft, "try-error"))
  print(tft[padj < 0.05][1:12, .(pathway, NES = round(NES, 2), padj = signif(padj, 2))])
cat("\n完成 -> T19..T27\n")
