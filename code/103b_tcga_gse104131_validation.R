#!/usr/bin/env Rscript
# =============================================================
# 103b_tcga_gse104131_validation.R  (KJR 投稿新增补充分析)
# 两层公共 RNA-seq 外部验证:
#   A. GSE104131 (Illumina HiSeq 2000, cufflinks FPKM, 16 对
#      肿瘤-癌旁配对, 13 AAM + 19 EAM?): 诊断外部验证
#      - 零参数 panel (队列内 z 均值) AUC + 95%CI
#      - 配对 DeLong vs 单基因
#      - 患者内配对 Wilcoxon (panel score tumor vs normal)
#      - 祖先敏感性 (EA-only / AA-only AUC)
#   B. TCGA-PRAD (pan-cancer atlas 2018, RSEM, 493 肿瘤):
#      - panel score 与 PFS 的 Cox (单变量 + 调整 T/N 分期 + LRT)
#      - panel score 与病理 T / N 分期的关联
#      - KM (panel score 中位分层)
# 输出:
#   tables/T45c_gse104131_auc.csv        T45d_gse104131_delong.csv
#   tables/T45e_gse104131_roc_curves.csv T45f_gse104131_paired.csv
#   tables/T45g_tcga_pfs_cox.csv         T45h_tcga_clinical_assoc.csv
#   tables/T45i_tcga_km_curves.csv       tables/T45j_tcga_panel_scores.csv
#   figures/FigS4_rnaseq_external.pdf / .jpg
# =============================================================
suppressPackageStartupMessages({
  library(data.table); library(pROC); library(survival); library(ggplot2)
  library(patchwork); library(ragg)
})
set.seed(2026)
root   <- "."
ext_d  <- file.path(root, "08_QC_Revision/data/external")
tcga_d <- file.path(ext_d, "tcga")
tab_dir <- file.path(root, "08_QC_Revision/tables")
fig_dir <- file.path(root, "08_QC_Revision/figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
panel4 <- c("AMACR", "PCA3", "HOXC6", "HPN")

get_field <- function(lines, key) {
  l <- grep(paste0("^", key, "\t"), lines, value = TRUE)
  if (!length(l)) return(NULL)
  gsub('^"|"$', "", strsplit(l[1], "\t")[[1]][-1])
}

# =============================================================
# A. GSE104131 — RNA-seq 诊断外部验证
# =============================================================
cat("========== GSE104131 ==========\n")
con <- gzfile(file.path(ext_d, "GSE104131_series_matrix.txt.gz"), "rt")
lines <- readLines(con, warn = FALSE); close(con)
title <- get_field(lines, "!Sample_title")
title <- sub("^Sample_", "", title)                      # -> 334_A7 ...
ch <- grep("^!Sample_characteristics_ch1\t", lines, value = TRUE)
chv <- lapply(ch, function(l) gsub('^"|"$', "", strsplit(l, "\t")[[1]][-1]))
tissue <- character(length(title)); ethnicity <- character(length(title))
for (v in chv) {
  if (all(grepl("^tissue:", v)))      tissue    <- sub("^tissue:\\s*", "", v)
  if (all(grepl("^ethnicity:", v)))   ethnicity <- sub("^ethnicity:\\s*", "", v)
}
stopifnot(length(tissue) == length(title))
meta104 <- data.table(sample = title, tissue = tissue, ethnicity = ethnicity)
meta104[, grp := fcase(tissue == "prostate tumor", "tumor",
                       tissue == "normal prostate", "normal",
                       default = NA_character_)]
cat("tissue label counts:", paste(names(table(meta104$grp, useNA = "ifany")),
    table(meta104$grp, useNA = "ifany"), sep = "=", collapse = ", "), "\n")
meta104[, patient := sub("_.*$", "", sample)]

# FPKM 矩阵 (batch1)
fp <- fread(file.path(ext_d, "GSE104131_batch1.fpkms.txt.gz"))
sym_col <- "GENESYMBOL"
scols <- setdiff(names(fp), c("GID", sym_col, "Locus"))
E104 <- as.matrix(fp[, ..scols]); rownames(E104) <- fp[[sym_col]]
cat("FPKM matrix:", nrow(E104), "genes x", ncol(E104), "samples\n")
stopifnot(all(panel4 %in% rownames(E104)))

# 共同样本 + 标签对齐 (FPKM 列名与 title 一致)
common <- intersect(colnames(E104), meta104$sample)
meta104 <- meta104[match(common, sample)]
E104 <- E104[, common, drop = FALSE]
cat("analytic set:", ncol(E104), "samples |", sum(meta104$grp == "tumor"),
    "tumour /", sum(meta104$grp == "normal"), "normal |\n")

# 零参数 panel: log2(FPKM+1) -> 队列内 z -> 均值
El <- log2(E104 + 1)
mu <- rowMeans(El); s <- apply(El, 1, sd); s[s == 0] <- NA
Z104 <- sweep(sweep(El, 1, mu), 1, s, "/")
meta104[, panel_score := colMeans(Z104[panel4, , drop = FALSE])]
for (g in panel4) meta104[[paste0("z_", g)]] <- Z104[g, ]

# AUC
auc_tab <- list(); rocs <- list()
for (nm in c("panel", panel4)) {
  v <- if (nm == "panel") meta104$panel_score else meta104[[paste0("z_", nm)]]
  r <- roc(meta104$grp, v, levels = c("normal", "tumor"), direction = "<", quiet = TRUE)
  rocs[[nm]] <- r
  ci <- as.numeric(ci.auc(r))
  auc_tab[[length(auc_tab) + 1]] <-
    data.table(model = if (nm == "panel") "Zero-param panel (fixed-weight)" else nm,
               short = nm, n = nrow(meta104),
               auc = round(ci[2], 3), ci_lo = round(ci[1], 3), ci_hi = round(ci[3], 3))
}
auc_tab <- rbindlist(auc_tab)

# 配对 DeLong: panel vs 单基因
dl_tab <- list()
for (g in panel4) {
  tt <- roc.test(rocs$panel, rocs[[g]], method = "delong", paired = TRUE)
  dl_tab[[length(dl_tab) + 1]] <-
    data.table(comparison = paste0("panel vs ", g),
               auc_panel = as.numeric(auc(rocs$panel)),
               auc_ref = as.numeric(auc(rocs[[g]])),
               delta_auc = round(as.numeric(auc(rocs$panel)) - as.numeric(auc(rocs[[g]])), 3),
               p_delong = signif(tt$p.value, 3))
}
dl_tab <- rbindlist(dl_tab)

# 患者内配对 Wilcoxon
p104 <- meta104[!is.na(grp), .(patient, grp, panel_score)]
w <- tidyr_pivot <- dcast(p104, patient ~ grp, value.var = "panel_score")
w <- w[!is.na(tumor) & !is.na(normal)]
pw <- suppressWarnings(wilcox.test(w$tumor, w$normal, paired = TRUE))
cat("paired Wilcoxon: n pairs =", nrow(w), ", P =", signif(pw$p.value, 3), "\n")

# 祖先敏感性
eth_tab <- list()
for (e in unique(meta104$ethnicity)) {
  sub104 <- meta104[ethnicity == e]
  if (length(unique(sub104$grp)) < 2) next
  r <- roc(sub104$grp, sub104$panel_score,
           levels = c("normal", "tumor"), direction = "<", quiet = TRUE)
  ci <- as.numeric(ci.auc(r))
  eth_tab[[length(eth_tab) + 1]] <-
    data.table(ethnicity = e, n = nrow(sub104),
               n_tumor = sum(sub104$grp == "tumor"), n_normal = sum(sub104$grp == "normal"),
               auc = round(ci[2], 3), ci_lo = round(ci[1], 3), ci_hi = round(ci[3], 3))
}
eth_tab <- rbindlist(eth_tab)
cat("\n=== GSE104131 诊断 AUC ===\n"); print(auc_tab)
print(eth_tab); print(dl_tab)

fwrite(auc_tab, file.path(tab_dir, "T45c_gse104131_auc.csv"))
fwrite(dl_tab, file.path(tab_dir, "T45d_gse104131_delong.csv"))
fwrite(eth_tab, file.path(tab_dir, "T45e_gse104131_ethnicity.csv"))
fwrite(rbindlist(lapply(names(rocs), function(nm) {
  r <- rocs[[nm]]
  data.table(model = nm, fpr = 1 - r$specificities, tpr = r$sensitivities,
             auc = as.numeric(auc(r)))
})), file.path(tab_dir, "T45f_gse104131_roc_curves.csv"))
fwrite(w[, .(patient, normal, tumor, delta = tumor - normal)],
       file.path(tab_dir, "T45g_gse104131_paired.csv"))
fwrite(meta104[, .(sample, patient, tissue, ethnicity, grp, panel_score)],
       file.path(tab_dir, "T45h_gse104131_scores.csv"))

# =============================================================
# B. TCGA-PRAD — 493 肿瘤, PFS 结局
# =============================================================
cat("\n========== TCGA-PRAD ==========\n")
ex <- fread(file.path(tcga_d, "T45_tcga_prad_expression.csv"))
cl <- fread(file.path(tcga_d, "T45b_tcga_prad_clinical.csv"))
ex[, patient := patientId]
cl[, patient := patientId]
m <- merge(ex, cl, by = "patient", all.x = TRUE)
cat("merged:", nrow(m), "samples | PFS_STATUS n =",
    sum(m$PFS_STATUS != ""), "\n")
m[, pfs_event := as.integer(startsWith(PFS_STATUS, "1"))]
m[, pfs_time := suppressWarnings(as.numeric(PFS_MONTHS))]
ok <- !is.na(m$pfs_time) & m$pfs_time > 0 & !is.na(m$pfs_event) & m$PFS_STATUS != ""
cat("PFS analytic set:", sum(ok), "| events =", sum(m$pfs_event[ok]), "\n")

# 零参数 panel: log2(RSEM+1) -> 每基因队列内 z -> 样本均值
Etcga <- as.matrix(m[, ..panel4])
Etl <- log2(Etcga + 1)
mu_g <- apply(Etl, 2, mean); sd_g <- apply(Etl, 2, sd); sd_g[sd_g == 0] <- NA
Zt2 <- sweep(sweep(Etl, 2, mu_g), 2, sd_g, "/")   # 按基因(列)标准化, 行(样本)为观测
m[, panel_score := rowMeans(Zt2[, panel4, drop = FALSE])]
cat("panel score range:", sprintf("%.2f", min(m$panel_score)), "-",
    sprintf("%.2f", max(m$panel_score)), "\n")

# Cox: PFS
mm <- m[ok]
mm[, t_stage_hi := as.integer(grepl("^T3|^T4", PATH_T_STAGE))]
mm[, n_stage_hi := as.integer(PATH_N_STAGE == "N1")]
cat("t_stage n =", sum(!is.na(mm$t_stage_hi)), "(T3/T4:", sum(mm$t_stage_hi == 1, na.rm = TRUE),
    ") | n_stage n =", sum(!is.na(mm$n_stage_hi)), "(N1:", sum(mm$n_stage_hi == 1, na.rm = TRUE), ")\n")
fit1 <- coxph(Surv(pfs_time, pfs_event) ~ panel_score, data = mm)
s1 <- summary(fit1)
# 调整模型: T/N 分期在 93 事件下稀疏 -> 系数分离风险, 用分层基线风险估计 panel 效应
fit0 <- coxph(Surv(pfs_time, pfs_event) ~ strata(t_stage_hi, n_stage_hi), data = mm)
fit2 <- coxph(Surv(pfs_time, pfs_event) ~ panel_score + strata(t_stage_hi, n_stage_hi), data = mm)
s2 <- summary(fit2)
lr <- anova(fit0, fit2)
cox_tab <- rbindlist(list(
  data.table(model = "panel (univariable)", n = nrow(mm), events = sum(mm$pfs_event),
             HR = round(s1$conf.int[1], 2), lo = round(s1$conf.int[3], 2),
             hi = round(s1$conf.int[4], 2), p = signif(s1$coefficients[5], 3)),
  data.table(model = "panel, stratified by T and N stage", n = nrow(mm), events = sum(mm$pfs_event),
             HR = round(s2$conf.int["panel_score", 1], 2),
             lo = round(s2$conf.int["panel_score", 3], 2),
             hi = round(s2$conf.int["panel_score", 4], 2),
             p = signif(s2$coefficients["panel_score", 5], 3)),
  data.table(model = "panel added to T/N strata (LRT)", n = nrow(mm), events = sum(mm$pfs_event),
             HR = NA_real_, lo = NA_real_, hi = NA_real_,
             p = signif(lr$`Pr(>|Chi|)`[2], 3))))
cat("\n=== TCGA-PRAD PFS Cox ===\n"); print(cox_tab)

# 临床关联: T 分期 / N 分期 (向量法, 自动丢弃 NA)
x1 <- mm[!is.na(t_stage_hi)]$panel_score; g1 <- mm[!is.na(t_stage_hi)]$t_stage_hi
x2 <- mm[!is.na(n_stage_hi)]$panel_score; g2 <- mm[!is.na(n_stage_hi)]$n_stage_hi
tt <- suppressWarnings(wilcox.test(x1[g1 == 1], x1[g1 == 0]))
nn <- suppressWarnings(wilcox.test(x2[g2 == 1], x2[g2 == 0]))
assoc <- rbindlist(list(
  data.table(comparison = "T3/T4 vs T2 (panel score)", n = length(x1),
             median_hi = round(median(x1[g1 == 1]), 2),
             median_lo = round(median(x1[g1 == 0]), 2),
             p = signif(tt$p.value, 3)),
  data.table(comparison = "N1 vs N0 (panel score)", n = length(x2),
             median_hi = round(median(x2[g2 == 1]), 2),
             median_lo = round(median(x2[g2 == 0]), 2),
             p = signif(nn$p.value, 3))))
cat("\n=== 临床关联 ===\n"); print(assoc)

# KM (panel 中位分层)
mm[, grp_km := fcase(panel_score >= median(panel_score, na.rm = TRUE), "panel high",
                     default = "panel low")]
kmf <- survfit(Surv(pfs_time, pfs_event) ~ grp_km, data = mm)
km_tab <- rbindlist(lapply(seq_along(kmf$strata), function(i) {
  st <- names(kmf$strata)[i]
  data.table(group = sub("^grp_km=", "", st), time = kmf$time, surv = kmf$surv, stratum = i)
}))
fwrite(km_tab, file.path(tab_dir, "T45i_tcga_km_curves.csv"))
fwrite(m[, .(sampleId, patient, AMACR, PCA3, HOXC6, HPN, panel_score,
             pfs_time, pfs_event, PFS_STATUS, PATH_T_STAGE, PATH_N_STAGE)],
       file.path(tab_dir, "T45j_tcga_panel_scores.csv"))
fwrite(cox_tab, file.path(tab_dir, "T45k_tcga_pfs_cox.csv"))
fwrite(assoc, file.path(tab_dir, "T45l_tcga_clinical_assoc.csv"))

# =============================================================
# FigS4 — 补充图: RNA-seq 外部验证
# =============================================================
theme_pub <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 8, hjust = 0),
        plot.subtitle = element_text(size = 7, hjust = 0),
        legend.position = "bottom",
        legend.title = element_text(size = 7), legend.text = element_text(size = 7),
        strip.background = element_rect(fill = "grey92", colour = NA))

roc104 <- fread(file.path(tab_dir, "T45f_gse104131_roc_curves.csv"))
roc104[, lab := sprintf("%s %.2f", model, auc)]
roc104[model == "panel", lab := sprintf("Panel %.2f", auc)]
roc104[model == "panel", model := "Panel"]
lv <- roc104[!duplicated(lab)]$lab
p1 <- ggplot(roc104, aes(fpr, tpr, colour = factor(lab, levels = lv))) +
  geom_path(linewidth = 0.4) + geom_abline(linetype = 2, colour = "grey55", linewidth = 0.3) +
  scale_colour_manual(values = c("#B2182B", "#2166AC", "#35978F", "#762A83", "#D95F02")) +
  guides(colour = guide_legend(ncol = 3, byrow = TRUE, override.aes = list(linewidth = 0.6))) +
  labs(title = "GSE104131 RNA-seq", subtitle = "15 paired tumour-normal specimens",
       x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
  theme_pub + coord_equal()

pw2 <- fread(file.path(tab_dir, "T45g_gse104131_paired.csv"))
wlong <- melt(pw2, id.vars = "patient", measure.vars = c("normal", "tumor"),
              variable.name = "grp", value.name = "score")
p2 <- ggplot(wlong, aes(grp, score, group = patient)) +
  geom_line(colour = "grey60", linewidth = 0.25) +
  geom_point(aes(colour = grp), size = 1) +
  scale_colour_manual(values = c(normal = "#8073AC", tumor = "#B2182B")) +
  labs(title = "Paired panel score", subtitle = sprintf("Wilcoxon P = %.4f", pw$p.value),
       x = NULL, y = "Zero-param panel score (z)", colour = NULL) +
  theme_pub

km <- fread(file.path(tab_dir, "T45i_tcga_km_curves.csv"))
km[, lab := group]
p3 <- ggplot(km, aes(time, surv, colour = lab)) +
  geom_step(linewidth = 0.45) +
  scale_colour_manual(values = c("panel high" = "#B2182B", "panel low" = "#2166AC")) +
  labs(title = "TCGA-PRAD progression-free survival", subtitle = sprintf("n = %d, 93 events; Cox P = %.2f", cox_tab[2]$n, cox_tab[2]$p),
       x = "Months", y = "Progression-free survival", colour = NULL) +
  theme_pub + ylim(0, 1)

ass <- fread(file.path(tab_dir, "T45l_tcga_clinical_assoc.csv"))
sc <- fread(file.path(tab_dir, "T45j_tcga_panel_scores.csv"))
sc[, Tgrp := fcase(PATH_T_STAGE %in% c("T2A", "T2B", "T2C"), "T2",
                   PATH_T_STAGE %in% c("T3A", "T3B", "T4"), "T3/T4",
                   default = NA_character_)]
p4 <- ggplot(sc[!is.na(Tgrp)], aes(Tgrp, panel_score, fill = Tgrp)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.3, width = 0.45) +
  scale_fill_manual(values = c("T2" = "#BDD7E7", "T3/T4" = "#3182BD")) +
  labs(title = "Panel score by pathologic T stage", subtitle = sprintf("Mann-Whitney P = %.2f", ass[comparison %like% "T3"]$p),
       x = NULL, y = "Zero-param panel score (z)", fill = NULL) +
  theme_pub + theme(legend.position = "none")

pw_fig <- (p1 | p2) / (p3 | p4) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

pdf(file.path(tempdir(), "FigS4_rnaseq_external.pdf"), width = 7.2, height = 5.4,
    onefile = TRUE)
print(pw_fig)
dev.off()
okc <- file.copy(file.path(tempdir(), "FigS4_rnaseq_external.pdf"),
                 file.path(fig_dir, "FigS4_rnaseq_external.pdf"), overwrite = TRUE)
cat("FigS4 saved:", okc, "\n")

# JPEG 600dpi (投稿嵌入用)
ragg::agg_jpeg(file.path(tempdir(), "FigS4_rnaseq_external.jpg"),
               width = 7.2 * 300, height = 5.4 * 300, res = 300)
print(pw_fig)
dev.off()
okj <- file.copy(file.path(tempdir(), "FigS4_rnaseq_external.jpg"),
                 file.path(fig_dir, "FigS4_rnaseq_external.jpg"), overwrite = TRUE)
cat("FigS4 jpg saved:", okj, "\n")
cat("DONE\n")
