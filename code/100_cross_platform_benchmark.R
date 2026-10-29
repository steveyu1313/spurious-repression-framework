#!/usr/bin/env Rscript
# =============================================================
# 100_cross_platform_benchmark.R  (JBI 投稿新增)
# 跨平台基准对比: 零参数固定权重模型 vs 监督学习
#   训练: GSE46602 (Affymetrix GPL570, 发现队列)
#   外测: GSE70770 Cambridge (Illumina HumanHT-12, 125T/74B)
#   模型: ZeroParam (within-cohort z 均值) / Logistic / LASSO(1-SE) /
#         Random Forest(500 树) / SVM-RBF
# 指标: 内部 AUC(监督=10折CV袋外; ZeroParam=无参数直接算),
#       外部 AUC(95%CI), 泛化衰减 ΔAUC, 配对 DeLong vs ZeroParam,
#       外部 Brier(ZeroParam 用训练集冻结的 2 参数 logistic 重标定),
#       外部校准斜率/截距
# 输出:
#   tables/T46_benchmark_summary.csv
#   tables/T46b_benchmark_roc_curves.csv
# =============================================================
suppressPackageStartupMessages({
  library(data.table); library(limma); library(hgu133plus2.db)
  library(illuminaHumanv4.db)
  library(AnnotationDbi); library(pROC); library(glmnet)
  library(randomForest); library(e1071)
})
set.seed(2026)
root    <- "."
data_d  <- file.path(root, "08_QC_Revision/data")
ext_d   <- file.path(root, "08_QC_Revision/data/external")
tab_dir <- file.path(root, "08_QC_Revision/tables")
panel4  <- c("AMACR", "PCA3", "HOXC6", "HPN")

# ---------------- series matrix 解析（同 80/82 脚本） ----------------
get_field <- function(lines, key) {
  l <- grep(paste0("^", key, "\t"), lines, value = TRUE)
  if (!length(l)) return(NULL)
  gsub('^"|"$', "", strsplit(l[1], "\t")[[1]][-1])
}
read_series <- function(path) {
  con <- gzfile(path, "rt"); lines <- readLines(con, warn = FALSE); close(con)
  gsm <- get_field(lines, "!Sample_geo_accession")
  ch_lines <- grep("^!Sample_characteristics_ch1\t", lines, value = TRUE)
  ch <- if (length(ch_lines)) {
    rows <- lapply(ch_lines, function(l) gsub('^"|"$', "", strsplit(l, "\t")[[1]][-1]))
    n <- min(lengths(rows)); rows <- lapply(rows, `[`, seq_len(n))
    keys <- unique(unlist(lapply(rows, function(v) sub(":.*$", "", v[nzchar(v) & !is.na(v)]))))
    mat <- matrix(NA_character_, nrow = n, ncol = length(keys), dimnames = list(NULL, keys))
    for (ri in seq_along(rows)) for (j in seq_len(n)) {
      kv <- rows[[ri]][j]
      if (is.na(kv) || !nzchar(kv)) next
      k <- sub(":.*$", "", kv)
      if (k %in% keys) mat[j, k] <- sub("^[^:]*:\\s*", "", kv)
    }
    as.data.frame(mat, stringsAsFactors = FALSE)
  } else data.frame(matrix(NA_character_, nrow = length(gsm), ncol = 0))
  i1 <- grep("^!series_matrix_table_begin", lines)
  i2 <- grep("^!series_matrix_table_end", lines)
  tab <- read.delim(text = lines[(i1 + 1):(i2 - 1)], header = TRUE,
                    row.names = 1, check.names = FALSE, quote = "\"")
  list(gsm = gsm, chars = ch, expr = as.matrix(tab))
}
to_symbol <- function(expr, db_name) {
  db <- get(db_name)
  sym <- suppressMessages(AnnotationDbi::mapIds(
    db, keys = rownames(expr), column = "SYMBOL",
    keytype = "PROBEID", multiVals = "first"))
  keep <- !is.na(sym); expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
  o <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
  expr <- expr[o, , drop = FALSE]; sym <- sym[o]
  expr <- expr[!duplicated(sym), , drop = FALSE]
  sym <- sym[!duplicated(sym)]
  rownames(expr) <- sym
  expr
}
norm_log <- function(e) {
  if (max(e, na.rm = TRUE) > 50) e <- log2(e + 1)
  normalizeBetweenArrays(e, method = "quantile")
}
zscore_rows <- function(e) {
  mu <- rowMeans(e); s <- apply(e, 1, sd); s[s == 0] <- NA
  sweep(sweep(e, 1, mu), 1, s, "/")
}

# ---------------- 训练队列 GSE46602 ----------------
raw_tr <- read_series(file.path(root, "02_Bulk_RNAseq/data/GSE46602_series_matrix.txt.gz"))
tr_grp <- ifelse(grepl("tumor", raw_tr$chars$tissue), "tumor", "benign")
E_tr <- norm_log(to_symbol(raw_tr$expr, "hgu133plus2.db")[, raw_tr$gsm, drop = FALSE])
stopifnot(all(panel4 %in% rownames(E_tr)))
cat("train GSE46602:", ncol(E_tr), "samples | tumour", sum(tr_grp == "tumor"),
    "/ benign", sum(tr_grp == "benign"), "\n")

# ---------------- 外测队列 GSE70770 Cambridge ----------------
raw_te <- read_series(file.path(ext_d, "GSE70770-GPL10558_series_matrix.txt.gz"))
te_meta <- raw_te$chars
te_grp <- ifelse(!is.na(te_meta$`sample type`) & te_meta$`sample type` == "Tumour",
                 "tumor", ifelse(!is.na(te_meta$`sample type`), "benign", NA))
E_te <- norm_log(to_symbol(raw_te$expr, "illuminaHumanv4.db")[, raw_te$gsm, drop = FALSE])
stopifnot(all(panel4 %in% rownames(E_te)))
ok <- !is.na(te_grp)
E_te <- E_te[, raw_te$gsm[ok], drop = FALSE]; te_grp <- te_grp[ok]
cat("test GSE70770 (Cambridge):", ncol(E_te), "samples | tumour", sum(te_grp == "tumor"),
    "/ benign", sum(te_grp == "benign"), "\n")

# ---------------- 设计矩阵 ----------------
Xtr <- t(E_tr[panel4, , drop = FALSE]); ytr <- factor(tr_grp, levels = c("benign", "tumor"))
Xte <- t(E_te[panel4, , drop = FALSE]); yte <- factor(te_grp, levels = c("benign", "tumor"))

# ---------------- 模型训练（监督） ----------------
fit_glm  <- suppressWarnings(glm(ytr ~ ., data = as.data.frame(Xtr), family = binomial))
cv_lasso <- cv.glmnet(as.matrix(Xtr), as.integer(ytr) - 1, family = "binomial",
                      alpha = 1, nfolds = 10, type.measure = "auc")
fit_rf   <- randomForest(x = Xtr, y = ytr, ntree = 500, importance = FALSE)
fit_svm  <- svm(x = Xtr, y = ytr, kernel = "radial", probability = TRUE, cost = 1)

# 内部 10 折 CV（监督模型袋外预测; ZeroParam 无参数, 全样本直接算）
# 判别用保序分数: Logistic/LASSO 用线性预测值(link), RF 用概率, SVM 用决策值 ——
# 避免 MLE 在分离数据上系数爆炸 -> plogis 饱和 -> 平坦并列毁掉排序
folds <- sample(rep(1:10, length.out = nrow(Xtr)))
cv_pred_mat <- matrix(NA_real_, nrow = nrow(Xtr), ncol = 4,
                      dimnames = list(NULL, c("Logistic", "LASSO", "RF", "SVM")))
for (f in 1:10) {
  tr_i <- folds != f; te_i <- folds == f
  dtr <- as.data.frame(Xtr[tr_i, , drop = FALSE]); dtr$y <- ytr[tr_i]
  g  <- suppressWarnings(glm(y ~ ., data = dtr, family = binomial))
  cv_pred_mat[te_i, "Logistic"] <- as.numeric(predict(g, as.data.frame(Xtr[te_i, , drop = FALSE]), type = "link"))
  l  <- cv.glmnet(as.matrix(Xtr[tr_i, , drop = FALSE]), as.integer(ytr[tr_i]) - 1,
                  family = "binomial", alpha = 1, nfolds = 5, type.measure = "auc")
  cv_pred_mat[te_i, "LASSO"] <- as.numeric(predict(l, as.matrix(Xtr[te_i, , drop = FALSE]),
                                                   s = "lambda.1se", type = "link"))
  rf <- randomForest(x = Xtr[tr_i, , drop = FALSE], y = ytr[tr_i], ntree = 500)
  cv_pred_mat[te_i, "RF"] <- predict(rf, Xtr[te_i, , drop = FALSE], type = "prob")[, "tumor"]
  sv <- svm(x = Xtr[tr_i, , drop = FALSE], y = ytr[tr_i], kernel = "radial",
            probability = TRUE, cost = 1)
  sp <- predict(sv, Xtr[te_i, , drop = FALSE], probability = TRUE, decision.values = TRUE)
  cv_pred_mat[te_i, "SVM"] <- as.numeric(attr(sp, "decision.values")[, 1])
}

# ---------------- 外部预测 ----------------
# score_ext: 保序分数（算 AUC / DeLong）; prob_ext: 概率（算 Brier / 校准）
score_ext <- list(
  Logistic = as.numeric(predict(fit_glm, as.data.frame(Xte), type = "link")),
  LASSO    = as.numeric(predict(cv_lasso, as.matrix(Xte), s = "lambda.1se", type = "link")),
  RF       = predict(fit_rf, Xte, type = "prob")[, "tumor"],
  SVM      = as.numeric(attr(predict(fit_svm, Xte, probability = TRUE,
                                     decision.values = TRUE), "decision.values")[, 1])
)
prob_ext <- list(
  Logistic = as.numeric(predict(fit_glm, as.data.frame(Xte), type = "response")),
  LASSO    = as.numeric(predict(cv_lasso, as.matrix(Xte), s = "lambda.1se", type = "response")),
  RF       = score_ext$RF,
  SVM      = attr(predict(fit_svm, Xte, probability = TRUE), "probabilities")[, "tumor"]
)
# ZeroParam: 各队列内部 z 均值（与全文一致的设计）
Ztr <- zscore_rows(E_tr); Zte <- zscore_rows(E_te)
zp_tr <- colMeans(Ztr[panel4, , drop = FALSE])
zp_te <- colMeans(Zte[panel4, , drop = FALSE])
score_ext$ZeroParam <- zp_te
# ZeroParam 概率化: 训练集冻结的 2 参数岭回归重标定
# （无惩罚 MLE 在分离数据上发散, 系数爆炸 -> 饱和; 岭惩罚保持有限、保序）
ridge_logit <- function(x, y, lambda = 1) {
  xz <- scale(x); xm <- attr(xz, "scaled:center"); xs <- attr(xz, "scaled:scale")
  f <- function(p) {
    eta <- p[1] + p[2] * as.vector(xz)
    -sum(y * eta - log(1 + exp(eta))) + lambda / 2 * p[2]^2
  }
  o <- optim(c(0, 0), f, method = "BFGS")
  c(intercept = unname(o$par[1] - o$par[2] * xm / xs), slope = unname(o$par[2] / xs))
}
cal_coefs <- ridge_logit(zp_tr, as.integer(ytr) - 1, lambda = 1)
prob_ext$ZeroParam <- as.numeric(plogis(cal_coefs[1] + cal_coefs[2] * zp_te))
zp_int_auc <- as.numeric(auc(roc(ytr, zp_tr, levels = c("benign", "tumor"), direction = "<", quiet = TRUE)))

# ---------------- 汇总 ----------------
models <- c("ZeroParam", "Logistic", "LASSO", "RF", "SVM")
lab <- c(ZeroParam = "Zero-parameter fixed-weight panel (proposed)",
         Logistic = "Logistic regression", LASSO = "LASSO (10-fold CV, 1-SE)",
         RF = "Random forest (500 trees)", SVM = "SVM (RBF kernel)")
roc_int <- roc(ytr, cv_pred_mat[, "Logistic"], quiet = TRUE)  # placeholder
int_auc <- c(ZeroParam = zp_int_auc,
             Logistic = as.numeric(auc(roc(ytr, cv_pred_mat[, "Logistic"], quiet = TRUE))),
             LASSO    = as.numeric(auc(roc(ytr, cv_pred_mat[, "LASSO"], quiet = TRUE))),
             RF       = as.numeric(auc(roc(ytr, cv_pred_mat[, "RF"], quiet = TRUE))),
             SVM      = as.numeric(auc(roc(ytr, cv_pred_mat[, "SVM"], quiet = TRUE))))
r_zp <- roc(yte, score_ext$ZeroParam, levels = c("benign", "tumor"), direction = "<", quiet = TRUE)
ci_zp <- as.numeric(ci.auc(r_zp))
rows <- list()
r_list <- list(ZeroParam = r_zp)
for (m in models[-1]) {
  r <- roc(yte, score_ext[[m]], levels = c("benign", "tumor"), direction = "<", quiet = TRUE)
  r_list[[m]] <- r
}
for (m in models) {
  r <- r_list[[m]]
  ci <- as.numeric(ci.auc(r))
  dl <- if (m == "ZeroParam") NA_real_ else
    roc.test(r_zp, r, method = "delong", paired = TRUE)$p.value
  brier <- mean((prob_ext[[m]] - (yte == "tumor"))^2)
  # 外部校准（概率 vs 结局的 logistic 回归斜率/截距）
  cf <- suppressWarnings(glm(yte ~ prob_ext[[m]], family = binomial))
  rows[[m]] <- data.table(
    model = lab[[m]], short = m, n_ext = length(yte),
    internal_auc = round(int_auc[[m]], 3),
    external_auc = round(ci[2], 3), ci_lo = round(ci[1], 3), ci_hi = round(ci[3], 3),
    auc_drop = round(int_auc[[m]] - ci[2], 3),
    p_delong_vs_zeroparam = if (is.na(dl)) NA else signif(dl, 3),
    brier_external = round(brier, 4),
    cal_slope_ext = round(unname(coef(cf)[2]), 2),
    cal_intercept_ext = round(unname(coef(cf)[1]), 2))
}
summ <- rbindlist(rows)
fwrite(summ, file.path(tab_dir, "T46_benchmark_summary.csv"))
roc_rows <- rbindlist(lapply(names(r_list), function(m) {
  r <- r_list[[m]]
  data.table(model = lab[[m]], short = m, fpr = 1 - r$specificities,
             tpr = r$sensitivities, auc = as.numeric(auc(r)), n = length(yte))
}))
fwrite(roc_rows, file.path(tab_dir, "T46b_benchmark_roc_curves.csv"))
cat("\n=== 跨平台基准（GSE46602 训练 -> GSE70770 Cambridge 外测）===\n")
print(summ)
cat("\nlambda.1se =", signif(cv_lasso$lambda.1se, 3), "\n")
