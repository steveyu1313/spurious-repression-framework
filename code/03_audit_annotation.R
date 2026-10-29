suppressPackageStartupMessages(library(data.table))
m <- fread("./03_SingleCell/results/tables/scRNA_metadata_integrated.csv")

cat("=== 聚类 -> 细胞类型 映射 ===\n")
map <- m[, .N, by = .(seurat_clusters, celltype_marker)][order(seurat_clusters, -N)]
print(map)
cat("\n每聚类唯一注释数(>1 表示混合):\n")
print(map[, .(n_type = uniqueN(celltype_marker), N = sum(N)), by = seurat_clusters])

cat("\n=== 聚类 x 样本 交叉(检查是否存在单样本驱动聚类) ===\n")
ct <- dcast(m[, .N, by = .(seurat_clusters, gse_sample)], seurat_clusters ~ gse_sample,
            value.var = "N", fill = 0)
print(ct)

cat("\n=== 每个聚类的优势样本占比 ===\n")
x <- m[, .N, by = .(seurat_clusters, gse_sample)]
x[, frac := N / sum(N), by = seurat_clusters]
print(x[order(seurat_clusters, -frac)][, head(.SD, 1), by = seurat_clusters])

cat("\n=== SingleR 结果 ===\n")
sr <- fread("./03_SingleCell/results/tables/scRNA_SingleR_by_cluster.csv")
print(sr)
cat("\nSingleR 与 marker 注释交叉:\n")
print(table(m$celltype_marker, m$celltype_singler))
