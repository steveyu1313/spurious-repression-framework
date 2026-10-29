#!/usr/bin/env Rscript
# =============================================================
# 50_validation_diagnostic.R
# 目的：把"单队列描述"升级为"三队列验证的诊断/复发分层标志物"
#   1. GSE46602（发现，36 瘤 + 14 良性，LMD）+ GSE69223（15 配对）
#      + GSE55945（13 + 8）三队列独立 limma
#   2. 跨队列方向一致的交集基因
#   3. 单基因 ROC / 预设 4 基因 panel / LASSO panel，发现→验证
#   4. GSE46602 自带 BCR + BCR-free time → Cox / KM 复发分层
# 说明：三队列均为 GPL570，探针层面可直接比较
# =============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(hgu133plus2.db)
  library(AnnotationDbi)
  library(pROC)
  library(glmnet)
  library(survival)
})

root    <- "."
data_d  <- file.path(root, "08_QC_Revision/data")
tab_dir <- file.path(root, "08_QC_Revision/tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------- 解析 series matrix ----------------
get_field <- function(lines, key) {
  l <- grep(paste0("^", key, "\t"), lines, value = TRUE)
  if (!length(l)) return(NULL)
  v <- strsplit(l[1], "\t")[[1]][-1]
  gsub('^"|"$', "", v)
}

read_series <- function(path, tag) {
  con <- gzfile(path, "rt"); lines <- readLines(con, warn = FALSE); close(con)

  gsm  <- get_field(lines, "!Sample_geo_accession")
  ttl  <- get_field(lines, "!Sample_title")
  plats <- get_field(lines, "!Series_platform_id")

  ch_lines <- grep("^!Sample_characteristics_ch1\t", lines, value = TRUE)
  # 保留键名：解析成 样本 × 字段 的宽表
  ch <- if (length(ch_lines)) {
    rows <- lapply(ch_lines, function(l) gsub('^"|"$', "", strsplit(l, "\t")[[1]][-1]))
    n <- min(lengths(rows)); rows <- lapply(rows, `[`, seq_len(n))
    keys <- unique(unlist(lapply(rows, function(v)
      sub(":.*$", "", v[nzchar(v) & !is.na(v)]))))
    mat <- matrix(NA_character_, nrow = n, ncol = length(keys),
                  dimnames = list(NULL, keys))
    # rows 是"特征行"，每行含 n 个样本的取值；mat 为 样本 × 字段
    for (ri in seq_along(rows)) {
      for (j in seq_len(n)) {
        kv <- rows[[ri]][j]
        if (is.na(kv) || !nzchar(kv)) next
        k <- sub(":.*$", "", kv)
        if (k %in% keys) mat[j, k] <- sub("^[^:]*:\\s*", "", kv)
      }
    }
    as.data.frame(mat, stringsAsFactors = FALSE)
  } else {
    data.frame(matrix(NA_character_, nrow = length(gsm), ncol = 0))
  }

  i1 <- grep("^!series_matrix_table_begin", lines)
  i2 <- grep("^!series_matrix_table_end", lines)
  tab <- read.delim(text = lines[(i1 + 1):(i2 - 1)], header = TRUE,
                    row.names = 1, check.names = FALSE, quote = "\"")

  list(tag = tag, gsm = gsm, title = ttl, platform = plats,
       chars = ch, expr = as.matrix(tab))
}

# ---------------- 探针 → 符号 ----------------
to_symbol <- function(expr) {
  sym <- suppressMessages(AnnotationDbi::mapIds(
    hgu133plus2.db, keys = rownames(expr), column = "SYMBOL",
    keytype = "PROBEID", multiVals = "first"))
  keep <- !is.na(sym)
  expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
  # 同一符号多探针：保留平均表达最高者
  o <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
  expr <- expr[o, , drop = FALSE]; sym <- sym[o]
  dup <- duplicated(sym)
  expr <- expr[!dup, , drop = FALSE]
  rownames(expr) <- sym[!dup]
  expr
}

# ---------------- 读取三队列 ----------------
cohorts <- list(
  GSE46602 = file.path(root, "02_Bulk_RNAseq/data/GSE46602_series_matrix.txt.gz"),
  GSE69223 = file.path(data_d, "GSE69223_series_matrix.txt.gz"),
  GSE55945 = file.path(data_d, "GSE55945_series_matrix.txt.gz")
)

raw <- lapply(names(cohorts), function(k) read_series(cohorts[[k]], k))
names(raw) <- names(cohorts)

# ---------------- 分组 ----------------
meta <- list()

# GSE46602：LMD 上皮，非配对；自带 BCR / BCR-free time / Gleason
m  <- raw$GSE46602
ch <- m$chars
getcol <- function(df, k, default = NA_character_) {
  if (!is.null(df[[k]])) df[[k]] else rep(default, nrow(df))
}
meta$GSE46602 <- data.table(
  gsm      = m$gsm,
  group    = ifelse(grepl("tumor", getcol(ch, "tissue")), "tumor", "benign"),
  patient  = NA_character_,
  bcr      = getcol(ch, "bcr"),
  bcr_time = suppressWarnings(as.numeric(getcol(ch, "bcr_free_time"))),
  gleason  = suppressWarnings(as.numeric(getcol(ch, "gleason_grade")))
)

# GSE69223：15 例配对（癌 vs 癌旁）
m  <- raw$GSE69223
ch <- m$chars
meta$GSE69223 <- data.table(
  gsm      = m$gsm,
  group    = ifelse(getcol(ch, "tissue") == "PCa", "tumor", "benign"),
  patient  = getcol(ch, "individual"),
  bcr      = NA_character_, bcr_time = NA_real_, gleason = NA_real_)

# GSE55945：13 癌 + 8 良性
m  <- raw$GSE55945
ch <- m$chars
meta$GSE55945 <- data.table(
  gsm      = m$gsm,
  group    = ifelse(grepl("^Malignant", getcol(ch, "tissue type")), "tumor", "benign"),
  patient  = NA_character_,
  bcr      = NA_character_, bcr_time = NA_real_, gleason = NA_real_)

cat("=== 队列概览 ===\n")
man <- rbindlist(lapply(names(meta), function(k) {
  md <- meta[[k]]
  data.table(cohort = k, platform = raw[[k]]$platform[1],
             n_total = nrow(md),
             n_tumor = sum(md$group == "tumor"),
             n_benign = sum(md$group == "benign"),
             n_paired = ifelse(all(!is.na(md$patient)), length(unique(md$patient)), NA_integer_),
             n_bcr_event = sum(md$bcr == "YES", na.rm = TRUE),
             n_bcr_time = sum(!is.na(md$bcr_time)))
}))
print(man)
fwrite(man, file.path(tab_dir, "T29_cohort_manifest.csv"))

# ---------------- 各队列独立 limma ----------------
run_limma <- function(tag) {
  md <- meta[[tag]]
  e  <- to_symbol(raw[[tag]]$expr)
  e  <- e[, md$gsm, drop = FALSE]
  # 若为原始非 log 值则 log2
  if (max(e, na.rm = TRUE) > 50) e <- log2(e + 1)
  e <- normalizeBetweenArrays(e, method = "quantile")

  grp <- factor(md$group, levels = c("benign", "tumor"))
  d <- data.frame(grp = grp, stringsAsFactors = FALSE)
  if (all(!is.na(md$patient))) {
    d$patient <- factor(md$patient)
    design <- model.matrix(~ patient + grp, data = d)
  } else {
    design <- model.matrix(~ grp, data = d)
  }
  fit <- eBayes(lmFit(e, design))
  co <- ncol(design)
  tt <- topTable(fit, coef = co, number = Inf, sort.by = "none")
  tt$symbol <- rownames(tt)
  as.data.table(tt)[, .(symbol, logFC, AveExpr, t, P.Value, adj.P.Val)]
}

degs <- lapply(names(meta), run_limma)
names(degs) <- names(meta)
for (k in names(degs))
  cat(sprintf("%s: %d 基因, 显著 %d (上调 %d / 下调 %d)\n", k, nrow(degs[[k]]),
              sum(degs[[k]]$adj.P.Val < 0.05 & abs(degs[[k]]$logFC) > 1),
              sum(degs[[k]]$adj.P.Val < 0.05 & degs[[k]]$logFC > 1),
              sum(degs[[k]]$adj.P.Val < 0.05 & degs[[k]]$logFC < -1)))

fwrite(rbindlist(lapply(names(degs), function(k) cbind(cohort = k, degs[[k]]))),
       file.path(tab_dir, "T30_cohort_DEG_all.csv"))

# ---------------- 跨队列交集 ----------------
merge_degs <- Reduce(function(a, b) merge(a, b, by = "symbol", suffixes = c("", ".y")),
                     lapply(names(degs), function(k) {
                       d <- degs[[k]][, .(symbol, logFC, adj.P.Val)]
                       setnames(d, c("symbol", paste0("logFC_", k), paste0("padj_", k)))
                     }))
setnames(merge_degs, gsub("\\.y$", "", names(merge_degs)))

fc_cols <- paste0("logFC_", names(degs))
p_cols  <- paste0("padj_",  names(degs))

sig <- function(d, dir) {
  ok <- rep(TRUE, nrow(d))
  for (i in seq_along(fc_cols)) {
    ok <- ok & d[[p_cols[i]]] < 0.05 &
      (if (dir == "up") d[[fc_cols[i]]] > 1 else d[[fc_cols[i]]] < -1)
  }
  ok
}
merge_degs[, up_all   := sig(.SD, "up")]
merge_degs[, down_all := sig(.SD, "down")]
merge_degs[, consistent := up_all | down_all]

# 方向一致但放宽 FC 阈值（P<0.05，同向）
merge_degs[, up_loose := Reduce(`&`, lapply(seq_along(fc_cols), function(i)
  merge_degs[[p_cols[i]]] < 0.05 & merge_degs[[fc_cols[i]]] > 0))]
merge_degs[, down_loose := Reduce(`&`, lapply(seq_along(fc_cols), function(i)
  merge_degs[[p_cols[i]]] < 0.05 & merge_degs[[fc_cols[i]]] < 0))]

cat(sprintf("\n跨队列交集: 严格(|FC|>1 三队列) %d 个（上调 %d / 下调 %d）；放宽 %d 个\n",
            sum(merge_degs$consistent), sum(merge_degs$up_all), sum(merge_degs$down_all),
            sum(merge_degs$up_loose | merge_degs$down_loose)))

fwrite(merge_degs[order(-(up_all | down_all), padj_GSE46602)],
       file.path(tab_dir, "T31_crosscohort_intersection.csv"))

cat("\n严格交集中上调:\n")
print(head(merge_degs[up_all == TRUE][order(padj_GSE46602), .(symbol, logFC_GSE46602, logFC_GSE69223, logFC_GSE55945)], 25))
cat("\n严格交集中下调:\n")
print(head(merge_degs[down_all == TRUE][order(padj_GSE46602), .(symbol, logFC_GSE46602, logFC_GSE69223, logFC_GSE55945)], 30))

# ---------------- 表达矩阵（各队列 z-score）----------------
get_expr <- function(tag, genes) {
  e <- to_symbol(raw[[tag]]$expr)
  md <- meta[[tag]]
  e <- e[, md$gsm, drop = FALSE]
  if (max(e, na.rm = TRUE) > 50) e <- log2(e + 1)
  e <- normalizeBetweenArrays(e, method = "quantile")
  g <- intersect(genes, rownames(e))
  z <- t(scale(t(e[g, , drop = FALSE])))
  list(z = z, y = factor(md$group, levels = c("benign", "tumor")), meta = md, raw = e[g, ])
}

# ---------------- 单基因 AUC ----------------
auc_one <- function(tag, genes) {
  o <- get_expr(tag, genes)
  rbindlist(lapply(rownames(o$z), function(g) {
    x <- as.numeric(o$z[g, ])
    r <- suppressMessages(roc(o$y, x, quiet = TRUE, direction = "<"))
    ci <- as.numeric(ci.auc(r))
    data.table(cohort = tag, symbol = g, n = length(o$y),
               auc = as.numeric(auc(r)), ci_lo = ci[1], ci_hi = ci[3], p = NA_real_)
  }))
}

cand <- merge_degs[consistent == TRUE, symbol]
cat("\n候选基因数:", length(cand), "\n")
auc_all <- rbindlist(lapply(names(meta), function(k) auc_one(k, cand)))
fwrite(auc_all, file.path(tab_dir, "T32_AUC_single_gene.csv"))

wide <- dcast(auc_all, symbol ~ cohort, value.var = "auc")
setnames(wide, sub("^GSE", "AUC_GSE", names(wide)))
wide[, min_auc := pmin(AUC_GSE46602, AUC_GSE69223, AUC_GSE55945)]
wide[, mean_auc := rowMeans(.SD), .SDcols = paste0("AUC_", names(meta))]
cat("\n单基因 AUC 稳健性最好的 15 个:\n")
print(head(wide[order(-min_auc)], 15))

# ---------------- 预设 4 基因 panel（生物学先验，无训练）----------------
panel4 <- c("AMACR", "PCA3", "KRT15", "TP63")
score4 <- function(tag) {
  o <- get_expr(tag, panel4)
  z <- o$z
  miss <- setdiff(panel4, rownames(z))
  sc <- rowMeans(cbind(z["AMACR", ], z["PCA3", ], -z["KRT15", ], -z["TP63", ]))
  data.table(cohort = tag, score = as.numeric(sc), y = o$y, gsm = colnames(z),
             missing = if (length(miss)) paste(miss, collapse = ",") else "")
}

res4 <- rbindlist(lapply(names(meta), score4))
auc4 <- rbindlist(lapply(names(meta), function(k) {
  d <- res4[cohort == k]
  r <- suppressMessages(roc(d$y, d$score, quiet = TRUE, direction = "<"))
  ci <- as.numeric(ci.auc(r))
  data.table(cohort = k, n = nrow(d), auc = as.numeric(auc(r)),
             ci_lo = ci[1], ci_hi = ci[3])
}))
auc4[, panel := "AMACR+PCA3-(KRT15+TP63)"]
cat("\n=== 预设 4 基因 panel ===\n"); print(auc4)

# ---------------- 共识上调基因 panel（不含基底角蛋白，bulk 可复现）----------------
panel5 <- c("AMACR", "PCA3", "HOXC6", "HPN")
res5 <- rbindlist(lapply(names(meta), function(tag) {
  o <- get_expr(tag, panel5)
  sc <- colMeans(o$z[intersect(panel5, rownames(o$z)), , drop = FALSE])
  data.table(cohort = tag, score = as.numeric(sc), y = o$y, gsm = colnames(o$z))
}))
auc5 <- rbindlist(lapply(names(meta), function(k) {
  d <- res5[cohort == k]
  r <- suppressMessages(roc(d$y, d$score, quiet = TRUE, direction = "<"))
  ci <- as.numeric(ci.auc(r))
  data.table(cohort = k, n = nrow(d), auc = as.numeric(auc(r)),
             ci_lo = ci[1], ci_hi = ci[3], panel = "AMACR+PCA3+HOXC6+HPN")
}))
cat("\n=== 共识上调基因 panel ===\n"); print(auc5)

# ---------------- LASSO panel（在 GSE46602 训练，另两队列验证）----------------
o466 <- get_expr("GSE46602", cand)
y466 <- as.numeric(o466$y == "tumor")
set.seed(20261020)
cvfit <- cv.glmnet(t(o466$z), y466, family = "binomial", alpha = 1, nfolds = 10)
coefs <- as.matrix(coef(cvfit, s = "lambda.1se"))
sel <- rownames(coefs)[coefs[, 1] != 0 & rownames(coefs) != "(Intercept)"]
w   <- coefs[sel, 1]
cat("\nLASSO 选中基因 (lambda.1se):", length(sel), "\n")
print(sort(w, decreasing = TRUE))

auc_lasso <- rbindlist(lapply(names(meta), function(k) {
  o <- get_expr(k, sel)
  z <- o$z[sel, , drop = FALSE]
  lp <- as.numeric(t(z) %*% w)
  r <- suppressMessages(roc(o$y, lp, quiet = TRUE, direction = "<"))
  ci <- as.numeric(ci.auc(r))
  data.table(cohort = k, n = length(o$y), auc = as.numeric(auc(r)),
             ci_lo = ci[1], ci_hi = ci[3], n_gene = length(sel))
}))
auc_lasso[, panel := sprintf("LASSO (%d genes, trained on GSE46602)", length(sel))]
cat("\n=== LASSO panel ===\n"); print(auc_lasso)

fwrite(rbindlist(list(auc4, auc5, auc_lasso), fill = TRUE), file.path(tab_dir, "T33_AUC_panel.csv"))
fwrite(res4, file.path(tab_dir, "T33b_panel4_scores.csv"))
fwrite(res5, file.path(tab_dir, "T33c_panel5_scores.csv"))

# ---------------- BCR 预后（GSE46602 肿瘤样本）----------------
tm <- meta$GSE46602[group == "tumor"]
surv_d <- tm[!is.na(bcr_time) & !is.na(bcr) & bcr %in% c("YES", "NO")]
cat(sprintf("\nBCR 可分析肿瘤样本: %d（事件 %d）\n", nrow(surv_d), sum(surv_d$bcr == "YES")))

o <- get_expr("GSE46602", unique(c(panel4, sel)))
gz <- o$z[, surv_d$gsm, drop = FALSE]
sc <- rowMeans(cbind(gz["AMACR", ], gz["PCA3", ], -gz["KRT15", ], -gz["TP63", ]))
surv_d[, score := as.numeric(sc)]
surv_d[, high := score > median(score)]

cox_panel <- coxph(Surv(bcr_time, bcr == "YES") ~ score, data = surv_d)
km <- survdiff(Surv(bcr_time, bcr == "YES") ~ high, data = surv_d)
p_lr <- 1 - pchisq(km$chisq, df = 1)
cat(sprintf("panel Cox: HR=%.2f (95%%CI %.2f-%.2f), P=%.4f | log-rank P=%.4f\n",
            exp(coef(cox_panel)), exp(confint(cox_panel))[1], exp(confint(cox_panel))[2],
            summary(cox_panel)$coefficients[5], p_lr))

cox_tab <- rbindlist(lapply(c(panel4, sel), function(g) {
  if (!g %in% rownames(gz)) return(NULL)
  d <- surv_d[gsm %in% colnames(gz)]
  d[, x := as.numeric(gz[g, d$gsm])]
  f <- coxph(Surv(bcr_time, bcr == "YES") ~ x, data = d)
  s <- summary(f)
  data.table(gene = g, HR = s$coefficients[2], lo = exp(confint(f))[1],
             hi = exp(confint(f))[2], p = s$coefficients[5])
}))
cox_tab <- rbind(data.table(gene = "PANEL_4gene", HR = exp(coef(cox_panel)),
                            lo = exp(confint(cox_panel))[1], hi = exp(confint(cox_panel))[2],
                            p = summary(cox_panel)$coefficients[5]), cox_tab)
cat("\n=== BCR Cox（单因素）===\n"); print(cox_tab[order(p)])
cat(sprintf("log-rank P = %.4f\n", p_lr))
fwrite(cox_tab, file.path(tab_dir, "T34_BCR_cox.csv"))
fwrite(surv_d[, .(gsm, group, bcr, bcr_time, gleason, score, high)],
       file.path(tab_dir, "T35_BCR_patient_level.csv"))

saveRDS(list(meta = meta, panel4 = panel4, panel5 = panel5, lasso = sel, w = w,
             auc5 = auc5, res5 = res5,
             surv = surv_d, auc4 = auc4, auc_lasso = auc_lasso,
             res4 = res4, wide = wide),
        file.path(tab_dir, "validation_objects.rds"))
cat("\n完成，输出 T29–T35\n")
