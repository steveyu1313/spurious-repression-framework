suppressPackageStartupMessages(library(data.table))
m <- fread("./04_Spatial/results/tables/spatial_metadata_integrated.csv")
cat("=== 每个 section 的 spot 计数分布 ===\n")
print(as.data.frame(m[, .(spots = .N, min_gene = min(nFeature_Spatial), q05_gene = quantile(nFeature_Spatial, .05),
            med_gene = median(nFeature_Spatial), max_gene = max(nFeature_Spatial),
            min_UMI = min(nCount_Spatial), med_UMI = median(nCount_Spatial)),
        by = sample][order(sample)]))
cat("\n基因数 <200 的 spot 占比:", round(mean(m$nFeature_Spatial < 200) * 100, 2), "%\n")
cat("count 为 0 的 spot 数:", sum(m$nCount_Spatial == 0), "\n")
cat("percent.mt 范围:", range(m$percent.mt), "\n")
cat("\n=== 唯一 section 数(按表达指纹) ===\n")
fp <- m[, .(medUMI = median(nCount_Spatial), medGene = median(nFeature_Spatial),
            medMT = median(percent.mt)), by = sample]
fp[, key := paste(medUMI, medGene, medMT)]
print(as.data.frame(fp))
cat("唯一指纹:", uniqueN(fp$key), "\n")
cat("\n=== 去重后应报告的样本量 ===\n")
cat("名义 section:", uniqueN(m$sample), " 唯一 section:", uniqueN(fp$key), "\n")
cat("名义 spot:", nrow(m), " 唯一 spot:", nrow(m) / 2, "\n")
