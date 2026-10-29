#!/usr/bin/env Rscript
# =============================================================
# 82_external_validation.R  (v2)
# 外部验证（GSE70770, Ross-Adams; 双中心, Illumina HumanHT-12 V4）:
#   剑桥子集      : 诊断 AUC（125 瘤 vs 74 良性）+ panel vs 诊断时 PSA（配对 DeLong）
#   斯德哥尔摩子集: 生存 Cox（94 瘤, 45 BCR 事件; 单变量 + PSA/Gleason 多变量）
# 输出: tables/T43_external_auc.csv, T43b_external_delong.csv,
#       T44_external_cox.csv
# =============================================================
suppressPackageStartupMessages({
  library(data.table); library(limma); library(pROC); library(survival)
  library(illuminaHumanv4.db); library(AnnotationDbi)
})

root    <- "."
ext_d   <- file.path(root, "08_QC_Revision/data/external")
tab_dir <- file.path(root, "08_QC_Revision/tables")
panel4  <- c("AMACR", "PCA3", "HOXC6", "HPN")

get_field <- function(lines, key) {
  l <- grep(paste0("^", key, "\t"), lines, value = TRUE)
  if (!length(l)) return(NULL)
  gsub('^"|"$', "", strsplit(l[1], "\t")[[1]][-1])
}
con <- gzfile(file.path(ext_d, "GSE70770-GPL10558_series_matrix.txt.gz"), "rt")
lines <- readLines(con, warn = FALSE); close(con)
gsm <- get_field(lines, "!Sample_geo_accession")
ch_lines <- grep("^!Sample_characteristics_ch1\t", lines, value = TRUE)
rows <- lapply(ch_lines, function(l) gsub('^"|"$', "", strsplit(l, "\t")[[1]][-1]))
keys <- unique(unlist(lapply(rows, function(v) sub(":.*$", "", v[nzchar(v) & !is.na(v)]))))
mat <- matrix(NA_character_, nrow = length(gsm), ncol = length(keys), dimnames = list(NULL, keys))
for (ri in seq_along(rows)) for (j in seq_along(rows[[ri]])) {
  kv <- rows[[ri]][j]
  if (is.na(kv) || !nzchar(kv)) next
  k <- sub(":.*$", "", kv)
  if (k %in% keys) mat[j, k] <- sub("^[^:]*:\\s*", "", kv)
}
meta <- as.data.table(mat); meta[, gsm := gsm]
setnames(meta, "biochemical relapse (bcr)", "bcr")
meta[, bcr := toupper(substr(bcr, 1, 1))]
i1 <- grep("^!series_matrix_table_begin", lines)
i2 <- grep("^!series_matrix_table_end", lines)
expr_tab <- read.delim(text = lines[(i1 + 1):(i2 - 1)], header = TRUE,
                       row.names = 1, check.names = FALSE, quote = "\"")
expr <- as.matrix(expr_tab)

sym <- suppressMessages(AnnotationDbi::mapIds(
  illuminaHumanv4.db, keys = rownames(expr), column = "SYMBOL",
  keytype = "PROBEID", multiVals = "first"))
keep <- !is.na(sym); expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
o <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
expr <- expr[o, , drop = FALSE]; sym <- sym[o]
expr <- expr[!duplicated(sym), , drop = FALSE]
rownames(expr) <- sym[!duplicated(sym)]
cat("panel genes present:", paste(intersect(panel4, rownames(expr)), collapse = ", "),
    "| missing:", paste(setdiff(panel4, rownames(expr)), collapse = ","), "\n")
if (max(expr, na.rm = TRUE) > 50) expr <- log2(expr + 1)
expr <- normalizeBetweenArrays(expr, method = "quantile")

z_of <- function(gsm_vec) {
  e <- expr[, gsm_vec, drop = FALSE]
  mu <- rowMeans(e); s <- apply(e, 1, sd); s[s == 0] <- NA
  sweep(sweep(e, 1, mu), 1, s, "/")
}

# ---------------- 剑桥子集：诊断 + panel vs PSA ----------------
cam <- meta[!is.na(`sample type`)]
cam[, grp := factor(`sample type` == "Tumour", levels = c(FALSE, TRUE), labels = c("benign", "tumor"))]
cat("Cambridge diagnostic: n =", nrow(cam), "-> tumour", sum(cam$grp == "tumor"),
    "/ benign", sum(cam$grp == "benign"), "| PSA available:", sum(!is.na(cam[["psa at diag"]])), "\n")
z_cam <- z_of(cam$gsm)
cam[, panel_score := colMeans(z_cam[panel4, , drop = FALSE])]
# PSA 只在癌患者记录（良性无 PSA），诊断 AUC 对比不成立；
# panel vs PSA 的增量检验放在生存模型（见下）。
auc_tab <- list(); dl_tab <- list()
rocs <- list()
for (nm in c("PANEL (fixed-weight)", panel4)) {
  v <- if (nm == "PANEL (fixed-weight)") cam$panel_score else z_cam[nm, ]
  ok <- !is.na(v)
  r <- roc(cam$grp[ok], v[ok], levels = c("benign", "tumor"), direction = "<", quiet = TRUE)
  rocs[[nm]] <- list(roc = r, ok = ok)
  ci <- as.numeric(ci.auc(r))
  auc_tab[[length(auc_tab) + 1]] <- data.table(cohort = "GSE70770 (Cambridge)", model = nm,
    n = sum(ok), auc = ci[2], ci_lo = ci[1], ci_hi = ci[3])
}
# 配对 DeLong: panel vs 各单基因（剑桥）
for (g in panel4) {
  okp <- rocs[["PANEL (fixed-weight)"]]$ok & rocs[[g]]$ok
  r1 <- roc(cam$grp[okp], cam$panel_score[okp], levels = c("benign", "tumor"), direction = "<", quiet = TRUE)
  r2 <- roc(cam$grp[okp], z_cam[g, cam$gsm[okp]], levels = c("benign", "tumor"), direction = "<", quiet = TRUE)
  tt <- roc.test(r1, r2, method = "delong", paired = TRUE)
  dl_tab[[length(dl_tab) + 1]] <- data.table(cohort = "GSE70770 (Cambridge)",
    comparison = paste0("panel vs ", g),
    auc_panel = as.numeric(auc(r1)), auc_ref = as.numeric(auc(r2)),
    delta_auc = as.numeric(auc(r1)) - as.numeric(auc(r2)),
    p_delong = signif(tt$p.value, 3), n = sum(okp))
}

# ---------------- 斯德哥尔摩子集：生存 ----------------
sth <- meta[is.na(`sample type`) & bcr %in% c("Y", "N")]
cat("Stockholm survival: n =", nrow(sth), ", events =", sum(sth$bcr == "Y"), "\n")
sth[, is_event := bcr == "Y"]
sth[, t_event := suppressWarnings(as.numeric(`time to bcr (months)`))]
sth[, t_cens  := suppressWarnings(as.numeric(`total follow up (months)`))]
cat("event time (raw, labelled months): range",
    sprintf("%.2f", min(sth[is_event == TRUE]$t_event, na.rm = TRUE)), "-",
    sprintf("%.2f", max(sth[is_event == TRUE]$t_event, na.rm = TRUE)),
    "| censor range",
    sprintf("%.2f", min(sth[is_event == FALSE]$t_cens, na.rm = TRUE)), "-",
    sprintf("%.2f", max(sth[is_event == FALSE]$t_cens, na.rm = TRUE)), "\n")
# 事件时间呈小数（0.36-<30），实为"年"被误标为 months；统一换算为月
sth[, time_months := fcase(is_event, t_event * 12, default = t_cens)]
sth[, event_i := as.integer(is_event)]
sth <- sth[!is.na(time_months) & time_months > 0]
cat("after cleaning: n =", nrow(sth), ", events =", sum(sth$event_i),
    "| time range", sprintf("%.1f-%.1f months", min(sth$time_months), max(sth$time_months)), "\n")
z_sth <- z_of(sth$gsm)
sth[, panel_score := colMeans(z_sth[panel4, , drop = FALSE])]
sth[, psa_c := pmin(as.numeric(`psa at diag`), 200)]
sth[, psa_c := ifelse(is.na(psa_c), median(psa_c, na.rm = TRUE), psa_c)]
sth[, gleason_n := suppressWarnings(as.numeric(sub("=.*$", "", `tumour gleason`)))]
cat("Gleason available:", sum(!is.na(sth$gleason_n)), "\n")

cox_rows <- list()
fit1 <- coxph(Surv(time_months, event_i) ~ panel_score, data = sth)
s1 <- summary(fit1)
cox_rows[[1]] <- data.table(cohort = "GSE70770 (Stockholm)", model = "panel (univariable)",
  n = nrow(sth), events = sum(sth$event_i), HR = s1$conf.int[1],
  lo = s1$conf.int[3], hi = s1$conf.int[4], p = signif(s1$coefficients[5], 3))
fit2 <- coxph(Surv(time_months, event_i) ~ panel_score + psa_c + gleason_n, data = sth)
s2 <- summary(fit2)
cox_rows[[2]] <- data.table(cohort = "GSE70770 (Stockholm)", model = "panel + PSA + Gleason",
  n = nrow(sth), events = sum(sth$event_i), HR = s2$conf.int["panel_score", 1],
  lo = s2$conf.int["panel_score", 3], hi = s2$conf.int["panel_score", 4],
  p = signif(s2$coefficients["panel_score", 5], 3))
fit3 <- coxph(Surv(time_months, event_i) ~ psa_c + gleason_n, data = sth)
s3 <- summary(fit3)
cox_rows[[3]] <- data.table(cohort = "GSE70770 (Stockholm)", model = "clinical only (PSA + Gleason)",
  n = nrow(sth), events = sum(sth$event_i), HR = NA_real_, lo = NA_real_, hi = NA_real_,
  p = signif(summary(fit3)$logtest["pvalue"], 3))
cox_tab <- rbindlist(cox_rows)
# 似然比检验: panel 增量（ clinical vs clinical+panel ）
lr <- anova(fit3, fit2)
cox_rows[[4]] <- data.table(cohort = "GSE70770 (Stockholm)",
  model = "panel added to PSA + Gleason (LRT)", n = nrow(sth), events = sum(sth$event_i),
  HR = NA_real_, lo = NA_real_, hi = NA_real_, p = signif(lr$`Pr(>|Chi|)`[2], 3))
cox_tab <- rbindlist(cox_rows)

# ---------------- 输出 ----------------
auc_all <- rbindlist(auc_tab)[order(model)]
fwrite(auc_all, file.path(tab_dir, "T43_external_auc.csv"))
fwrite(rbindlist(dl_tab), file.path(tab_dir, "T43b_external_delong.csv"))
fwrite(cox_tab, file.path(tab_dir, "T44_external_cox.csv"))
# ROC 曲线坐标（供 Fig6E）
roc_rows <- lapply(names(rocs), function(nm) {
  r <- rocs[[nm]]$roc
  data.table(model = nm, fpr = 1 - r$specificities, tpr = r$sensitivities,
             auc = as.numeric(auc(r)), n = sum(rocs[[nm]]$ok))
})
fwrite(rbindlist(roc_rows), file.path(tab_dir, "T43c_external_roc_curves.csv"))
cat("\n=== 外部诊断 AUC（GSE70770 剑桥）===\n"); print(auc_all)
cat("\n=== panel vs PSA ===\n"); print(rbindlist(dl_tab))
cat("\n=== 外部生存 Cox（GSE70770 斯德哥尔摩）===\n"); print(cox_tab)
