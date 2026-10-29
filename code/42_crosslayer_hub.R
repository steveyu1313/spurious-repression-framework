# ============================================================
# 42_crosslayer_hub.R
# 目标: hub 基因在单细胞各细胞型中的表达 -> 判断 bulk 信号来源
# 输出: T28_hub_scrna_by_celltype.csv
# 运行: cd 项目根 && Rscript 08_QC_Revision/scripts/42_crosslayer_hub.R
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(Seurat); library(Matrix)
})
tab_dir <- "./08_QC_Revision/tables"
hub <- fread(file.path(tab_dir, "T17_hub_genes.csv"))$hub
extra <- c("KLK3", "AR", "AMACR", "KRT5", "KRT14", "TP63", "PCA3", "ERG", "MSMB", "PTPRC", "ACTA2")
genes <- unique(c(hub, extra))

sc <- readRDS("./03_SingleCell/data/scrna_final.RDS")
cat("对象:", ncol(sc), "细胞 x", nrow(sc), "基因; assays:", paste(Assays(sc), collapse = ","), "\n")
DefaultAssay(sc) <- if ("RNA" %in% Assays(sc)) "RNA" else Assays(sc)[1]

meta <- sc@meta.data
ct_col <- intersect(c("celltype_marker", "celltype", "cell_type", "SingleR.labels"), names(meta))[1]
if (is.na(ct_col)) stop("未找到细胞类型列")
cat("细胞类型列:", ct_col, "\n"); print(table(meta[[ct_col]]))

genes <- intersect(genes, rownames(sc))
cat("可查到基因", length(genes), "/", length(unique(c(hub, extra))), "\n")

sc <- subset(sc, cells = colnames(sc)[!is.na(meta[[ct_col]])])
# Seurat 会把 identity 中的下划线替换为短横线, 先统一命名避免后续列名错位
sc$ct <- gsub("_", "-", as.character(sc@meta.data[[ct_col]]))
Idents(sc) <- "ct"
avg <- AverageExpression(sc, features = genes, assays = DefaultAssay(sc),
                         slot = "data", verbose = FALSE)[[1]]
avg <- as.data.table(as.data.frame(avg), keep.rownames = "gene")

# 表达比例
cnt <- LayerData(sc, assay = DefaultAssay(sc), layer = "counts")[genes, , drop = FALSE]
grp <- sc$ct
pct <- sapply(levels(factor(grp)), function(g) {
  m <- cnt[, grp == g, drop = FALSE]
  Matrix::rowMeans(m > 0) * 100
})
pct <- as.data.table(as.data.frame(pct), keep.rownames = "gene")

setnames(pct, setdiff(names(pct), "gene"), paste0("pct_", setdiff(names(pct), "gene")))
out <- merge(avg, pct, by = "gene")
out[, is_hub := gene %in% hub]
setorder(out, -is_hub, gene)
fwrite(out, file.path(tab_dir, "T28_hub_scrna_by_celltype.csv"))

cat("\n=== hub 基因最高表达细胞型 ===\n")
ctcols <- setdiff(names(avg), "gene")
top <- out[is_hub == TRUE]
top[, `:=`(max_celltype = ctcols[max.col(as.matrix(.SD), ties.method = "first")]),
    .SDcols = ctcols]
top[, max_avg := apply(.SD, 1, max), .SDcols = ctcols]
print(top[order(-max_avg), .(gene, max_celltype, max_avg = round(max_avg, 2))][1:20])
cat("\n完成 -> T28\n")
