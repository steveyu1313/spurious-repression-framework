suppressPackageStartupMessages({library(data.table); library(limma)})
B <- readRDS("./08_QC_Revision/tables/bulk_recheck.rds")
mat <- B$mat; grp <- B$group
run <- function(keep) {
  design <- model.matrix(~0 + factor(grp[keep], levels = c("Benign", "Tumor")))
  colnames(design) <- c("Benign", "Tumor")
  fit <- eBayes(contrasts.fit(lmFit(mat[, keep], design), makeContrasts(Tumor - Benign, levels = design)))
  tt <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
  tt$Symbol <- rownames(tt)
  tt
}
full <- run(rep(TRUE, ncol(mat)))
out  <- !(colnames(mat) == "GSM1133140")
drop <- run(out)

sig <- function(x) x[x$adj.P.Val < 0.05 & abs(x$logFC) > 1, "Symbol"]
cat("全样本 DEG:", length(sig(full)), " (up", sum(full[full$adj.P.Val<0.05 & abs(full$logFC)>1,"logFC"]>0),
    ", down", sum(full[full$adj.P.Val<0.05 & abs(full$logFC)>1,"logFC"]<0), ")\n")
cat("剔除 GSM1133140 后 DEG:", length(sig(drop)), "\n")
cat("重叠:", length(intersect(sig(full), sig(drop))), "\n")
m <- merge(full[, c("Symbol","logFC")], drop[, c("Symbol","logFC")], by = "Symbol", suffixes = c("_full","_drop"))
cat("logFC 相关 r =", round(cor(m$logFC_full, m$logFC_drop), 4), " 最大绝对差 =",
    round(max(abs(m$logFC_full - m$logFC_drop)), 3), "\n")
fwrite(m, "./08_QC_Revision/tables/T11_bulk_sensitivity_outlier.csv")
