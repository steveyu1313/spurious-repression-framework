# 52_panel_scrna.R — 四个 panel 基因在单细胞各细胞型的表达（验证摘要断言）
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix)
})
root <- "."
tab_dir <- file.path(root, "08_QC_Revision/tables")

sc <- readRDS(file.path(root, "03_SingleCell/data/scrna_final.RDS"))
cat("cells:", ncol(sc), " features:", nrow(sc), "\n")
cat("meta cols:", paste(grep("celltype|cluster|sample", names(sc@meta.data), value = TRUE), collapse = ", "), "\n")

ct_col <- intersect(c("celltype_marker", "celltype_singler", "celltype", "CellType"), names(sc@meta.data))[1]
if (is.na(ct_col)) stop("no celltype column")
ct <- sc@meta.data[[ct_col]]
cat("using column:", ct_col, "\n")
print(table(ct))

genes <- c("AMACR", "PCA3", "HOXC6", "HPN", "KLK3")
genes <- intersect(genes, rownames(sc))
cat("found:", paste(genes, collapse = ", "), "\n")

avg <- AverageExpression(sc, features = genes, assays = DefaultAssay(sc),
                         group.by = ct_col, layer = "data", verbose = FALSE)[[1]]
cnt <- LayerData(sc, assay = DefaultAssay(sc), layer = "counts")[genes, , drop = FALSE]
pct <- sapply(levels(factor(ct)), function(g) {
  m <- cnt[, ct == g, drop = FALSE]
  Matrix::rowMeans(m > 0) * 100
})

out <- data.table(gene = genes,
                  as.data.table(as.data.frame(pct), keep.rownames = FALSE))
avg_dt <- data.table(gene = rownames(avg), as.data.table(as.data.frame(avg)))
setnames(out, c("gene", colnames(pct)))
setnames(avg_dt, c("gene", colnames(avg)))

fwrite(out, file.path(tab_dir, "X9_panel_scrna_pct.csv"))
fwrite(avg_dt, file.path(tab_dir, "X10_panel_scrna_avg.csv"))

cat("\n=== 检出率 % ===\n"); print(out)
cat("\n=== 平均表达 (data layer) ===\n"); print(avg_dt)

# 每个基因表达最高的细胞型
for (g in genes) {
  v <- as.numeric(avg_dt[gene == g, -1])
  names(v) <- colnames(avg)[seq_along(v)]
  cat(sprintf("%-6s 最高: %-16s (%.2f) | 次高: %s (%.2f)\n", g,
              names(which.max(v)), max(v), names(sort(v, decreasing = TRUE))[2],
              sort(v, decreasing = TRUE)[2]))
}
