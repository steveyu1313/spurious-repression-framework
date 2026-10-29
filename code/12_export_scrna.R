# 从 Seurat 对象导出轻量绘图数据 (单细胞)
suppressPackageStartupMessages({ library(Seurat); library(data.table) })
out <- "./08_QC_Revision/tables"

obj <- readRDS("./03_SingleCell/data/scrna_final.RDS")
cat("loaded:", ncol(obj), "cells\n")
cat("assays:", paste(Assays(obj), collapse = ","), "\n")
print(head(obj[[]], 2))

um <- Embeddings(obj, "umap")
cat("umap dim:", dim(um), "\n")
md <- data.table(barcode = rownames(um), UMAP1 = um[, 1], UMAP2 = um[, 2])
cols <- intersect(c("gse_sample", "seurat_clusters", "celltype_marker",
                    "celltype_singler", "nCount_RNA", "nFeature_RNA"), names(obj[[]]))
cat("export cols:", paste(cols, collapse = ","), "\n")
md <- cbind(md, as.data.table(obj[[]])[, ..cols])
fwrite(md, file.path(out, "X1_scrna_umap.csv"))

# 关键基因表达 (SCT data layer)
genes <- c("KLK3", "AR", "AMACR", "KRT5", "PTPRC", "ACTA2", "CD68", "PECAM1", "MKI67", "MS4A1")
genes <- intersect(genes, rownames(obj))
cat("genes found:", paste(genes, collapse = ","), "\n")
ex <- FetchData(obj, vars = genes, layer = "data")
ex <- data.table(barcode = rownames(ex), ex)
fwrite(ex, file.path(out, "X2_scrna_gene_expr.csv"))

# 每细胞类型 marker 平均表达 (用于 dotplot 复算)
Idents(obj) <- "celltype_marker"
mk <- c("KLK3","KLK2","ACPP","AR","MSMB","KRT5","KRT15","TP63","PECAM1","VWF","CD3D","CD3E",
        "CD68","C1QC","MS4A1","CD79A","COL1A1","DCN","ACTA2","MYH11","TPSAB1","CPA3","MKI67","TOP2A")
mk <- intersect(mk, rownames(obj))
cat("markers found:", length(mk), "\n")
avg <- AverageExpression(obj, features = mk, assays = "SCT", layer = "data")$SCT
fwrite(data.table(gene = rownames(avg), as.data.table(as.matrix(avg))), file.path(out, "X3_scrna_marker_avg.csv"))
pct <- sapply(mk, function(g) tapply(FetchData(obj, vars = g, layer = "data")[[1]],
                                     obj$celltype_marker, function(x) mean(x > 0)))
fwrite(data.table(gene = rownames(pct), as.data.table(as.matrix(pct))), file.path(out, "X4_scrna_marker_pct.csv"))
cat("DONE scRNA export\n")
