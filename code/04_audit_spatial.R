suppressPackageStartupMessages(library(data.table))
out <- "./08_QC_Revision"
m <- fread("./04_Spatial/results/tables/spatial_metadata_integrated.csv")
cat("=== spatial metadata ===\n")
cat("dims:", nrow(m), "x", ncol(m), "\n")
cat("cols:", paste(names(m), collapse = ", "), "\n\n")
print(head(m, 2))

for (cc in c("seurat_clusters", "gse_sample", "Sample", "orig.ident", "sample")) {
  if (cc %in% names(m)) {
    cat("\n---", cc, ": n_unique =", uniqueN(m[[cc]]), "---\n")
    print(m[, .N, by = cc][order(-N)])
  }
}

cat("\n=== 聚类数核对 (README 称 27) ===\n")
cat("cluster 编号数:", uniqueN(m$seurat_clusters), "\n")

cat("\n=== 每样本 spot 数 ===\n")
sf <- intersect(c("gse_sample", "Sample", "sample", "orig.ident"), names(m))[1]
if (!is.na(sf)) {
  print(m[, .N, by = sf])
  cat("总 spot:", nrow(m), "\n")
  cat("理论上限 8 x 4992 =", 8 * 4992, "\n")
}

cat("\n=== 重复样本检测: 每个 spot 的表达指纹是否完全一致 ===\n")
# 用各样本 median nCount / nFeature / percent.mt 指纹判定是否为同一张切片
num_cols <- intersect(c("nCount_SCT", "nFeature_SCT", "nCount_RNA", "nFeature_RNA", "percent.mt"), names(m))
if (!is.na(sf) && length(num_cols)) {
  fp <- m[, lapply(.SD, function(x) round(median(x, na.rm = TRUE), 4)), by = sf, .SDcols = num_cols]
  print(fp)
  key <- apply(as.data.frame(fp)[, num_cols, drop = FALSE], 1, paste, collapse = "|")
  cat("\n指纹分组:\n"); print(split(fp[[sf]], key))
  cat("唯一指纹数:", length(unique(key)), "/", nrow(fp), "\n")
}

cat("\n=== cluster x sample 组成 ===\n")
if (!is.na(sf)) {
  print(dcast(m[, .N, by = c("seurat_clusters", sf)], seurat_clusters ~ get(sf),
              value.var = "N", fill = 0))
}
