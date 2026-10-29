# =============================================================================
# 01_audit_scRNA.R — 单细胞数据合理性核查
# =============================================================================
suppressPackageStartupMessages(library(data.table))
out <- "./08_QC_Revision"
dir.create(file.path(out, "tables"), showWarnings = FALSE, recursive = TRUE)

m <- fread("./03_SingleCell/results/tables/scRNA_metadata_integrated.csv")
cat("================ 1. metadata 结构 ================\n")
cat("dims:", nrow(m), "x", ncol(m), "\n")
cat("columns:", paste(names(m), collapse = ", "), "\n\n")
print(head(m, 3))

for (cc in c("seurat_clusters", "celltype_marker", "gse_sample", "Sample",
             "orig.ident", "celltype_SingleR", "singler_label")) {
  if (cc %in% names(m)) {
    cat("\n---", cc, ": n_unique =", uniqueN(m[[cc]]), "---\n")
    print(m[, .N, by = cc][order(-N)])
  }
}

cat("\n================ 2. 聚类数核对 (正文称 21) ================\n")
if ("seurat_clusters" %in% names(m)) {
  cl <- sort(unique(m$seurat_clusters))
  cat("cluster 编号:", paste(cl, collapse = ", "), "\n")
  cat("实际聚类数:", length(cl), "\n")
  cat("编号范围:", min(cl), "-", max(cl), "\n")
}

cat("\n================ 3. 细胞总数核对 (正文 36,424) ================\n")
cat("metadata 行数:", nrow(m), "\n")
if ("celltype_marker" %in% names(m)) {
  tb <- m[, .N, by = celltype_marker][order(-N)]
  tb[, pct := round(100 * N / sum(N), 2)]
  print(tb)
  cat("\n合计:", sum(tb$N), "\n")
  fwrite(tb, file.path(out, "tables/T3_scRNA_celltype_counts.csv"))
}

cat("\n================ 4. 每样本细胞数 vs QC 表 ================\n")
qcfile <- "./03_SingleCell/results/tables/scRNA_QC_stats.csv"
if (file.exists(qcfile)) {
  qc <- fread(qcfile)
  print(qc)
  cat("QC 表合计:", sum(qc$cells), "\n")
  samp_col <- intersect(c("gse_sample", "Sample", "orig.ident"), names(m))[1]
  if (!is.na(samp_col)) {
    per <- m[, .N, by = samp_col]
    setnames(per, "N", "n_meta")
    mm <- merge(qc[, .(sample, cells)], per, by.x = "sample", by.y = samp_col, all = TRUE)
    mm[, diff := n_meta - cells]
    print(mm)
    cat("标签错位(diff != 0)样本数:", mm[!is.na(diff) & diff != 0, .N], "\n")
    fwrite(mm, file.path(out, "tables/T4_scRNA_sample_label_check.csv"))
  }
}

cat("\n================ 5. QC 指标合理性 ================\n")
if (exists("qc")) {
  cat("median_mt 全部为 0 的样本数:", qc[median_mt == 0, .N], "/", nrow(qc), "\n")
  cat("median_UMI 范围:", range(qc$median_UMI), "\n")
  cat("median_genes 范围:", range(qc$median_genes), "\n")
  qc[, genes_per_umi := median_genes / median_UMI]
  print(qc[, .(sample, cells, median_UMI, median_genes, ratio = round(genes_per_umi, 2))])
  cat("典型 10x 数据 genes/UMI 参考区间约 0.3-0.5; 超出者提示指标异常\n")
}

cat("\n================ 6. 模块评分与注释一致性 ================\n")
ms <- "./03_SingleCell/results/tables/scRNA_module_scores_by_cluster.csv"
if (file.exists(ms)) {
  s <- fread(ms)
  cat("dims:", nrow(s), "x", ncol(s), "\n")
  print(names(s))
  print(head(s, 3))
}
sr <- "./03_SingleCell/results/tables/scRNA_SingleR_by_cluster.csv"
if (file.exists(sr)) {
  s2 <- fread(sr)
  cat("\nSingleR dims:", nrow(s2), "x", ncol(s2), "\n")
  print(head(s2, 5))
}
