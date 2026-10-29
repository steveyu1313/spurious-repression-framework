# 从空间 Seurat 对象导出轻量绘图数据
suppressPackageStartupMessages({ library(Seurat); library(data.table) })
out <- "./08_QC_Revision/tables"

obj <- readRDS("./04_Spatial/data/Visium/spatial_integrated.RDS")
cat("loaded:", ncol(obj), "spots\n")
cat("assays:", paste(Assays(obj), collapse = ","), "\n")
cat("images:", paste(Images(obj), collapse = ","), "\n")
cat("meta cols:", paste(names(obj[[]]), collapse = ","), "\n")

# 表达
genes <- intersect(c("KLK3", "AR", "AMACR", "PTPRC", "ACTA2", "KRT5", "ERG"), rownames(obj))
cat("genes:", paste(genes, collapse = ","), "\n")
ex <- FetchData(obj, vars = genes, layer = "data")
ex <- data.table(barcode = rownames(ex), ex)
fwrite(ex, file.path(out, "X5_spatial_gene_expr.csv"))

# 坐标: 逐 image 提取
coords <- list()
for (im in Images(obj)) {
  cc <- tryCatch(GetTissueCoordinates(obj, image = im), error = function(e) NULL)
  if (is.null(cc)) next
  cc$barcode <- rownames(cc)
  cc$section <- im
  coords[[im]] <- as.data.table(cc)
  cat("image", im, "spots:", nrow(cc), "\n")
}
co <- rbindlist(coords, fill = TRUE)
fwrite(co, file.path(out, "X6_spatial_coords.csv"))

md <- data.table(barcode = colnames(obj), sample = obj$sample,
                 seurat_clusters = obj$seurat_clusters,
                 nCount = obj$nCount_Spatial, nFeature = obj$nFeature_Spatial,
                 mt = obj$percent.mt)
fwrite(md, file.path(out, "X7_spatial_meta.csv"))
cat("DONE spatial export\n")
