# =============================================================================
# 15_spatial_reanalysis.R — 空间数据按规范重跑
#   1) 只保留 in_tissue == 1 的 spot (原分析误用全部 4,992 个 spot)
#   2) 只使用 4 个唯一 section (GSE245388 中 JLR/PRS/RJZ/SLX 矩阵与 CDZ/HB/LGL/CFX 完全重复)
#   3) SCTransform -> Harmony(按 section) -> 聚类 -> UMAP
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(harmony); library(Matrix)
})
set.seed(42)
base <- "./04_Spatial/data/Visium"
out  <- "./08_QC_Revision/tables"

sections <- c("CDZ", "CFX", "HB", "LGL")   # 4 个唯一 section

read_one <- function(s) {
  d <- file.path(base, s, "outs", "filtered_feature_bc_matrix")
  cnt <- Read10X(d, gene.column = 2)
  pos <- fread(file.path(base, s, "outs", "spatial", "tissue_positions_list.csv"), header = FALSE)
  setnames(pos, c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  keep <- pos[in_tissue == 1, barcode]
  keep <- intersect(keep, colnames(cnt))
  cat(sprintf("  %s: total %d, in-tissue %d, retained %d\n", s, ncol(cnt), sum(pos$in_tissue == 1), length(keep)))
  o <- CreateSeuratObject(cnt[, keep], project = s, min.cells = 3, min.features = 100)
  o$section <- s
  o$percent.mt <- PercentageFeatureSet(o, pattern = "^MT-")
  p <- pos[barcode %in% keep]
  o$array_row <- p$array_row[match(colnames(o), p$barcode)]
  o$array_col <- p$array_col[match(colnames(o), p$barcode)]
  list(obj = o, n_total = ncol(cnt))
}

cat("=== 读取并筛除组织外 spot ===\n")
lst <- lapply(sections, read_one)
objs <- lapply(lst, `[[`, "obj")

cat("\n=== QC 概览 (in-tissue subspots) ===\n")
qc <- rbindlist(lapply(objs, function(o) data.table(
  section = o$section[1], spots = ncol(o),
  med_UMI = median(o$nCount_RNA), med_gene = median(o$nFeature_RNA),
  med_mt = round(median(o$percent.mt), 2),
  min_gene = min(o$nFeature_RNA))))
print(as.data.frame(qc))
fwrite(qc, file.path(out, "T8_spatial_QC_corrected.csv"))

cat("\n=== SCTransform + Harmony ===\n")
objs <- lapply(objs, SCTransform, verbose = FALSE)
feats <- SelectIntegrationFeatures(objs, nfeatures = 2000)
m <- merge(objs[[1]], y = objs[-1])
m <- RunPCA(m, features = feats, npcs = 20, verbose = FALSE)
m <- RunHarmony(m, group.by.vars = "section", reduction.use = "pca",
                reduction.save = "harmony", verbose = FALSE)
m <- FindNeighbors(m, reduction = "harmony", dims = 1:15)
m <- FindClusters(m, resolution = 0.5, verbose = FALSE)
m <- RunUMAP(m, reduction = "harmony", dims = 1:15, verbose = FALSE)
cat("clusters:", nlevels(m$seurat_clusters), "\n")
print(table(m$seurat_clusters))

# 导出
um <- Embeddings(m, "umap")
md <- data.table(barcode = rownames(um), UMAP1 = um[, 1], UMAP2 = um[, 2],
                 section = m$section, cluster = as.character(m$seurat_clusters),
                 array_row = m$array_row, array_col = m$array_col,
                 nCount = m$nCount_RNA, nFeature = m$nFeature_RNA, mt = m$percent.mt)
genes <- intersect(c("KLK3", "AR", "AMACR", "PTPRC", "ACTA2", "KRT5", "ERG", "MSMB"), rownames(m))
ex <- as.data.frame(FetchData(m, vars = genes, layer = "data"))
md <- cbind(md, as.data.table(ex))

# 空间域注释: 依据 marker 均值
Idents(m) <- "seurat_clusters"
av <- AverageExpression(m, features = genes, assays = "SCT", layer = "data")$SCT
fwrite(data.table(cluster = colnames(av), t(round(as.matrix(av), 2))),
       file.path(out, "T9_spatial_cluster_markers.csv"))
fwrite(md, file.path(out, "X8_spatial_reanalysed.csv"))
cat("\nDONE spatial reanalysis, spots:", nrow(md), "\n")
