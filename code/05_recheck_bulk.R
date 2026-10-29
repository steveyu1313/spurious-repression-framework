# =============================================================================
# 05_recheck_bulk.R — GSE46602 独立复算 (limma)
# 目的: 从 series_matrix 重新走一遍流程, 验证 704 DEG / 260 上调 / 444 下调
# =============================================================================
suppressPackageStartupMessages({
  library(limma); library(hgu133plus2.db); library(data.table)
})
out <- "./08_QC_Revision"

# ---- 1. 解析 series_matrix ----
rows <- readLines(gzfile("./02_Bulk_RNAseq/data/GSE46602_series_matrix.txt.gz"))
meta_lines <- rows[grepl("^!Sample_", rows)]
meta <- list()
for (ln in meta_lines) {
  field <- sub("^!([^\t]+).*", "\\1", ln)
  vals <- strsplit(ln, "\t")[[1]][-1]
  vals <- gsub('"', "", trimws(vals))
  meta[[field]] <- c(meta[[field]], list(vals))
}
samples <- meta[["Sample_geo_accession"]][[1]]
char_list <- meta[["Sample_characteristics_ch1"]]
tissue <- char_list[[grep("^tissue:", vapply(char_list, `[`, character(1), 1))[1]]]
group <- ifelse(grepl("tumor", tissue, ignore.case = TRUE), "Tumor",
                ifelse(grepl("benign", tissue, ignore.case = TRUE), "Benign", NA))

bi <- grep("^!series_matrix_table_begin", rows)[1]
be <- grep("^!series_matrix_table_end", rows)[1]
expr <- read.table(text = rows[(bi + 1):(be - 1)], sep = "\t", header = TRUE,
                   row.names = 1, check.names = FALSE)
cat("=== 数据 ===\n")
cat("probes:", nrow(expr), " samples:", ncol(expr), "\n")
cat("分组:\n"); print(table(group))
chk <- data.table(sample = samples, group = group)
fwrite(chk, file.path(out, "tables/T5_GSE46602_group_assignment.csv"))

# ---- 2. 探针注释 + 基因层聚合 ----
sym <- AnnotationDbi::mapIds(hgu133plus2.db, keys = rownames(expr), column = "SYMBOL",
                             keytype = "PROBEID", multiVals = "first")
keep <- !is.na(sym)
e <- as.data.frame(expr[keep, ]); e$Symbol <- sym[keep]
agg <- aggregate(e[, setdiff(names(e), "Symbol")], by = list(Symbol = e$Symbol), mean)
mat <- as.matrix(agg[, -1]); rownames(mat) <- agg$Symbol
cat("\ngene-level matrix:", nrow(mat), "x", ncol(mat), "\n")

# ---- 3. 标准化 + limma ----
mat <- normalizeBetweenArrays(mat, method = "quantile")
g <- factor(ifelse(group == "Tumor", "Tumor", "Benign"), levels = c("Benign", "Tumor"))
design <- model.matrix(~0 + g); colnames(design) <- levels(g)
fit <- eBayes(contrasts.fit(lmFit(mat, design), makeContrasts(Tumor - Benign, levels = design)))
res <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
res$Symbol <- rownames(res)
sig <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ]
cat("\n=== 复算结果 ===\n")
cat("DEG:", nrow(sig), " up:", sum(sig$logFC > 0), " down:", sum(sig$logFC < 0), "\n")
fwrite(as.data.table(res), file.path(out, "tables/T6_GSE46602_DEG_rechecked.csv"))

# ---- 4. 与原结果比对 ----
old <- fread("./02_Bulk_RNAseq/results/tables/GSE46602_DEG_all.csv")
m <- merge(old[, .(Symbol, logFC_old = logFC, p_old = adj.P.Val)],
           as.data.table(res)[, .(Symbol, logFC_new = logFC, p_new = adj.P.Val)],
           by = "Symbol")
cat("\n=== 与原结果比对 ===\n")
cat("共同基因:", nrow(m), "\n")
cat("logFC 相关系数 r =", round(cor(m$logFC_old, m$logFC_new), 4), "\n")
cat("logFC 最大绝对差 =", round(max(abs(m$logFC_old - m$logFC_new)), 4), "\n")
cat("-log10P 相关系数 r =", round(cor(-log10(m$p_old), -log10(m$p_new)), 4), "\n")
oldsig <- old[adj.P.Val < 0.05 & abs(logFC) > 1, Symbol]
newsig <- sig$Symbol
cat("原 DEG 数:", length(oldsig), " 复算 DEG 数:", length(newsig), "\n")
cat("交集:", length(intersect(oldsig, newsig)),
    " 仅原结果:", length(setdiff(oldsig, newsig)),
    " 仅复算:", length(setdiff(newsig, oldsig)), "\n")

# ---- 5. 标志基因方向核验 ----
mk <- c("AMACR", "PCA3", "ERG", "MYC", "GSTP1", "PTEN", "KLK3", "AR", "NKX3-1",
        "KRT5", "KRT15", "KRT17", "TP63", "HOXC6", "FOLH1", "SLC45A3")
mk <- intersect(mk, res$Symbol)
cat("\n=== 标志基因 (复算) ===\n")
print(res[mk, c("Symbol", "logFC", "adj.P.Val")][order(-res[mk, "logFC"]), ], digits = 3)

# ---- 6. 设计效应量检查 ----
cat("\n=== 效应量分布 ===\n")
print(summary(abs(sig$logFC)))
cat("|logFC|>2 的 DEG:", sum(abs(sig$logFC) > 2), "\n")
cat("adj.P < 0.01 的 DEG:", sum(sig$adj.P.Val < 0.01), "\n")

# 保存复算矩阵供绘图
saveRDS(list(mat = mat, group = group, res = as.data.table(res), sig = sig),
        file.path(out, "tables/bulk_recheck.rds"))
cat("\n完成\n")
