suppressPackageStartupMessages(library(data.table))
m <- fread("./03_SingleCell/results/tables/scRNA_metadata_integrated.csv")
cat("percent.mt summary:\n"); print(summary(m$percent.mt))
cat("n distinct mt values:", uniqueN(m$percent.mt), "\n")
cat("frac mt>25:", mean(m$percent.mt > 25), "\n\n")
cat("per gse_sample medians:\n")
print(m[, .(cells = .N, medUMI = median(nCount_RNA), medGene = median(nFeature_RNA),
            medMT = median(percent.mt)),
        by = gse_sample][order(as.integer(sub("S", "", gse_sample)))])
cat("\nnCount_RNA range:", range(m$nCount_RNA), " nFeature range:", range(m$nFeature_RNA), "\n")

cat("\n=== 匹配检验: metadata 每样本中位UMI 与 QC 表 median_UMI 的对应 ===\n")
qc <- fread("./03_SingleCell/results/tables/scRNA_QC_stats.csv")
per <- m[, .(cells = .N, medUMI = median(nCount_RNA)), by = gse_sample]
for (i in seq_len(nrow(per))) {
  hit <- qc[abs(medUMI - per$medUMI[i]) < 1]
  cat(sprintf("metadata %-4s (cells=%5d, medUMI=%8.1f) -> QC 表中匹配样本: %s\n",
              per$gse_sample[i], per$cells[i], per$medUMI[i],
              if (nrow(hit)) paste(hit$sample, collapse = ",") else "无"))
}
