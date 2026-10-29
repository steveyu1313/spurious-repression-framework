#!/usr/bin/env Rscript
# =============================================================
# 80_incremental_value.R
# 审稿人式补强：
#   1. 增量价值：panel vs 单基因（AMACR/PCA3/HOXC6/HPN）配对 DeLong 检验
#   2. AUC 95% CI 已有；校准曲线（3 队列）+ 决策曲线分析（DCA, 手工净获益）
# 输出:
#   tables/T40_incremental_value.csv   — 单基因与 panel AUC(95%CI) + 配对 DeLong P
#   tables/T41_dca_netbenefit.csv      — DCA 净获益曲线数据
#   tables/T42_calibration.csv         — 校准统计（斜率/截距/Hosmer-Lemeshow）
#   figures/FigS3_calibration_dca.png  — 上排校准, 下排 DCA
# =============================================================
suppressPackageStartupMessages({
  library(data.table); library(limma); library(hgu133plus2.db)
  library(AnnotationDbi); library(pROC); library(ggplot2)
})

root    <- "."
data_d  <- file.path(root, "08_QC_Revision/data")
tab_dir <- file.path(root, "08_QC_Revision/tables")
fig_dir <- file.path(root, "08_QC_Revision/figures")

# ---------------- 解析 series matrix（同 50 脚本） ----------------
get_field <- function(lines, key) {
  l <- grep(paste0("^", key, "\t"), lines, value = TRUE)
  if (!length(l)) return(NULL)
  v <- strsplit(l[1], "\t")[[1]][-1]
  gsub('^"|"$', "", v)
}
read_series <- function(path, tag) {
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
  list(tag = tag, gsm = gsm, chars = ch, expr = as.matrix(tab))
}
to_symbol <- function(expr) {
  sym <- suppressMessages(AnnotationDbi::mapIds(
    hgu133plus2.db, keys = rownames(expr), column = "SYMBOL",
    keytype = "PROBEID", multiVals = "first"))
  keep <- !is.na(sym); expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
  o <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
  expr <- expr[o, , drop = FALSE]; sym <- sym[o]
  expr <- expr[!duplicated(sym), , drop = FALSE]
  rownames(expr) <- sym[!duplicated(sym)]
  expr
}

cohorts <- list(
  GSE46602 = file.path(root, "02_Bulk_RNAseq/data/GSE46602_series_matrix.txt.gz"),
  GSE69223 = file.path(data_d, "GSE69223_series_matrix.txt.gz"),
  GSE55945 = file.path(data_d, "GSE55945_series_matrix.txt.gz")
)
raw <- lapply(names(cohorts), function(k) read_series(cohorts[[k]], k))
names(raw) <- names(cohorts)

meta <- list()
m <- raw$GSE46602; ch <- m$chars
gc1 <- function(df, k) if (!is.null(df[[k]])) df[[k]] else rep(NA_character_, nrow(df))
meta$GSE46602 <- data.table(gsm = m$gsm,
  group = ifelse(grepl("tumor", gc1(ch, "tissue")), "tumor", "benign"))
m <- raw$GSE69223; ch <- m$chars
meta$GSE69223 <- data.table(gsm = m$gsm,
  group = ifelse(gc1(ch, "tissue") == "PCa", "tumor", "benign"))
m <- raw$GSE55945; ch <- m$chars
meta$GSE55945 <- data.table(gsm = m$gsm,
  group = ifelse(grepl("^Malignant", gc1(ch, "tissue type")), "tumor", "benign"))

# ---------------- 各队列 z-score 表达 ----------------
prep <- lapply(names(raw), function(k) {
  md <- meta[[k]]
  e <- to_symbol(raw[[k]]$expr)[, md$gsm, drop = FALSE]
  if (max(e, na.rm = TRUE) > 50) e <- log2(e + 1)
  e <- normalizeBetweenArrays(e, method = "quantile")
  mu <- rowMeans(e); s <- apply(e, 1, sd); s[s == 0] <- NA
  z <- sweep(sweep(e, 1, mu), 1, s, "/")
  list(e = e, z = z, md = md)
})
names(prep) <- names(raw)

panel4 <- c("AMACR", "PCA3", "HOXC6", "HPN")

# ---------------- AUC + 配对 DeLong ----------------
res_list <- list(); dl_list <- list()
for (k in names(prep)) {
  p  <- prep[[k]]
  y  <- factor(p$md$group, levels = c("benign", "tumor"))
  z  <- p$z
  scores <- list(panel = colMeans(z[panel4, , drop = FALSE]))
  for (g in panel4) scores[[g]] <- z[g, ]
  auc_row <- list()
  rocs <- list()
  for (nm in names(scores)) {
    v <- scores[[nm]][!is.na(scores[[nm]])]
    yy <- y[!is.na(scores[[nm]])]
    r <- roc(response = yy, predictor = v, levels = c("benign", "tumor"),
             direction = "<", quiet = TRUE)
    rocs[[nm]] <- r
    ci <- as.numeric(ci.auc(r))
    auc_row[[nm]] <- data.table(cohort = k, model = nm, n = length(yy),
                                auc = ci[2], ci_lo = ci[1], ci_hi = ci[3])
  }
  res_list[[k]] <- rbindlist(auc_row)
  # 配对 DeLong：panel vs 各单基因
  for (g in panel4) {
    tt <- roc.test(rocs$panel, rocs[[g]], method = "delong", paired = TRUE)
    dl_list[[length(dl_list) + 1]] <- data.table(
      cohort = k, comparison = paste0("panel vs ", g),
      auc_panel = as.numeric(auc(rocs$panel)), auc_ref = as.numeric(auc(rocs[[g]])),
      delta_auc = as.numeric(auc(rocs$panel)) - as.numeric(auc(rocs[[g]])),
      p_delong = tt$p.value)
  }
}
res <- rbindlist(res_list); res$model <- factor(res$model,
  levels = c("panel", panel4), labels = c("PANEL (fixed-weight)", panel4))
res <- res[order(cohort, model)]
fwrite(res, file.path(tab_dir, "T40_incremental_value.csv"))
dl <- rbindlist(dl_list); dl$p_delong <- signif(dl$p_delong, 3)
fwrite(dl, file.path(tab_dir, "T40b_delong_tests.csv"))
cat("=== AUC（95%CI）===\n"); print(res)
cat("\n=== 配对 DeLong（panel vs 单基因）===\n"); print(dl)

# ---------------- 校准 + DCA ----------------
cal_rows <- list(); hl_rows <- list(); dca_rows <- list()
for (k in names(prep)) {
  p <- prep[[k]]; y <- as.integer(p$md$group == "tumor")
  sc <- colMeans(p$z[panel4, , drop = FALSE]); ok <- !is.na(sc)
  sc <- sc[ok]; yy <- y[ok]
  fit <- suppressWarnings(glm(yy ~ sc, family = binomial))
  pr <- predict(fit, type = "response")
  sep <- !fit$converged || any(pr < 1e-6) || any(pr > 1 - 1e-6)
  # 校准：分组(分位)观测 vs 预期；完全分离时不可估
  if (!sep) {
    qb <- unique(quantile(pr, probs = 0:4/4))
    q  <- cut(pr, breaks = qb, include.lowest = TRUE)
    agg <- aggregate(cbind(obs = yy, pred = pr) ~ q, FUN = mean)
    cal_rows[[k]] <- data.table(cohort = k, predicted = agg$pred, observed = agg$obs, n = as.integer(table(q)))
    hl_stat <- sum((table(q) * (agg$obs - agg$pred)^2) / (agg$pred * (1 - agg$pred)))
    hl_p <- pchisq(hl_stat, df = max(nrow(agg) - 2, 1), lower.tail = FALSE)
    hl_rows[[k]] <- data.table(cohort = k, slope = coef(fit)[2], intercept = coef(fit)[1],
                               HL_stat = hl_stat, HL_df = max(nrow(agg) - 2, 1), HL_p = hl_p)
  } else {
    hl_rows[[k]] <- data.table(cohort = k, slope = NA_real_, intercept = NA_real_,
                               HL_stat = NA_real_, HL_df = NA_integer_, HL_p = NA_real_,
                               note = "degenerate: complete separation (AUC = 1); calibration not estimable")
    cal_rows[[k]] <- data.table(cohort = k, predicted = NA_real_, observed = NA_real_, n = NA_integer_)
  }
  # DCA: 阈值 1%-99%
  ths <- seq(0.01, 0.99, by = 0.01); N <- length(yy)
  prev <- mean(yy)
  for (th in ths) {
    nb_model <- mean((pr >= th) & (yy == 1)) - mean((pr >= th) & (yy == 0)) * (th / (1 - th))
    nb_all   <- prev - (1 - prev) * (th / (1 - th))
    dca_rows[[length(dca_rows) + 1]] <- data.table(cohort = k, threshold = th,
      nb_model = nb_model, nb_treat_all = nb_all, nb_treat_none = 0)
  }
}
cal <- rbindlist(cal_rows, fill = TRUE); hl <- rbindlist(hl_rows, fill = TRUE); dca <- rbindlist(dca_rows)
fwrite(cal, file.path(tab_dir, "T42_calibration.csv"))
fwrite(hl,  file.path(tab_dir, "T42b_calibration_stats.csv"))
fwrite(dca, file.path(tab_dir, "T41_dca_netbenefit.csv"))
cat("\n=== 校准统计 ===\n"); print(hl)

# ---------------- 图 FigS3 ----------------
cal_p <- cal[!is.na(predicted)]
pal <- c(GSE46602 = "#1B3A57", GSE69223 = "#3D8B8B", GSE55945 = "#C94F4F")
p1 <- ggplot(cal_p, aes(predicted, observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(colour = cohort), size = 2.4) +
  geom_line(aes(colour = cohort), linewidth = 0.6) +
  scale_colour_manual(values = pal) +
  xlim(0, 1) + ylim(0, 1) +
  labs(x = "Predicted probability of tumour", y = "Observed proportion of tumour",
       title = "Calibration of the fixed-weight panel (GSE69223)",
       subtitle = "In GSE46602 and GSE55945 the classes separate completely (AUC = 1) and calibration is not estimable") +
  theme_bw(base_size = 9) + theme(legend.position = "none",
  plot.title = element_text(face = "bold", size = 10))

dca_l <- melt(dca, id.vars = c("cohort", "threshold"),
              measure.vars = c("nb_model", "nb_treat_all", "nb_treat_none"))
dca_l[, model := fcase(variable == "nb_model", "Four-gene panel",
                       variable == "nb_treat_all", "Treat all",
                       default = "Treat none")]
p2 <- ggplot(dca_l, aes(threshold, value, colour = model, linewidth = variable == "nb_model")) +
  geom_line() + facet_wrap(~cohort, nrow = 1) +
  scale_colour_manual(values = c("Four-gene panel" = "#1B3A57",
                                 "Treat all" = "grey55", "Treat none" = "grey80")) +
  scale_linewidth_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.45), guide = "none") +
  coord_cartesian(ylim = c(-0.05, 0.80)) +
  labs(x = "Threshold probability", y = "Net benefit",
       title = "Decision-curve analysis",
       subtitle = "Net benefit of biopsying according to the panel score") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom",
  plot.title = element_text(face = "bold", size = 10))

tmp_png <- file.path(tempdir(), "FigS3_calibration_dca.png")
ggsave(tmp_png,
  gridExtra::grid.arrange(p1, p2, nrow = 2), width = 7.2, height = 5.6, dpi = 600)
file.copy(tmp_png, file.path(fig_dir, "FigS3_calibration_dca.png"), overwrite = TRUE)
cat("FigS3 written\n")
