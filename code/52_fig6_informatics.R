# =============================================================================
# 52_fig6_informatics.R — JBI 版 Figure 6（诊断性能与临床效用）
#   A  三队列 ROC（固定权重 panel）
#   B  跨平台基准 ROC：零参数 vs Logistic/LASSO/RF/SVM（GSE46602 训练 -> Cambridge 外测）
#   C  外部队列 ROC：panel vs 单基因（GSE70770 Cambridge）
#   D  校准曲线（GSE69223；其余队列完全分离不可估）
#   E  决策曲线分析 DCA（三队列）
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
  library(ragg); library(pROC)
})
FONT <- "Arial"
fig_dir <- "./08_QC_Revision/figures"
tab_dir <- "./08_QC_Revision/tables"
theme_pub <- function(base = 7) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks.length = unit(1.2, "mm"),
          axis.text = element_text(colour = "black", size = base),
          axis.title = element_text(colour = "black", size = base + 0.5),
          plot.title = element_text(size = base + 1, face = "bold", hjust = 0),
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = base - 0.8),
          legend.title = element_text(size = base - 0.5),
          plot.margin = margin(2, 3, 2, 2))
}
save_fig <- function(p, file, w, h, dpi = 600) {
  tmp <- file.path(tempdir(), basename(file))
  ggsave(tmp, p, width = w, height = h, dpi = dpi, device = ragg::agg_png, bg = "white")
  file.copy(tmp, file, overwrite = TRUE)
}
COH <- c(GSE46602 = "#C0392B", GSE69223 = "#2E86AB", GSE55945 = "#E08214")

# ---------------- A: 三队列 ROC ----------------
sc <- fread(file.path(tab_dir, "T33c_panel5_scores.csv"))
rocs <- rbindlist(lapply(unique(sc$cohort), function(k) {
  s <- sc[cohort == k]
  r <- suppressMessages(roc(s$y, s$score, quiet = TRUE, direction = "<"))
  data.table(cohort = k, fpr = 1 - r$specificities, tpr = r$sensitivities,
             auc = as.numeric(auc(r)), n = nrow(s))
}))
lab <- unique(rocs[, .(cohort, auc, n)])[order(match(cohort, names(COH)))]
lab[, txt := sprintf("%s (n=%d)  AUC %.2f", cohort, n, auc)]
lab_txt <- setNames(lab$txt, lab$cohort)
pA <- ggplot(rocs, aes(fpr, tpr, colour = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70", linewidth = 0.3) +
  geom_step(linewidth = 0.75) +
  scale_colour_manual(values = COH, labels = lab_txt, name = NULL) +
  scale_x_continuous(expand = c(0.01, 0), labels = percent_format(accuracy = 1)) +
  scale_y_continuous(expand = c(0.01, 0), labels = percent_format(accuracy = 1)) +
  labs(title = "Within-cohort performance",
       subtitle = "Fixed-weight 4-gene panel, z-score mean",
       x = "1 - specificity", y = "Sensitivity") +
  theme_pub() +
  theme(legend.position = c(0.99, 0.03), legend.justification = c(1, 0),
        legend.text = element_text(size = 5.2),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

# ---------------- B: 跨平台基准 ROC（零参数 vs 监督） ----------------
bm <- fread(file.path(tab_dir, "T46b_benchmark_roc_curves.csv"))
bm[, short := factor(short, levels = c("ZeroParam", "Logistic", "LASSO", "RF", "SVM"))]
bm[, short := factor(short, levels = rev(levels(short)))]  # 图例自上而下
PAL_BM <- c(ZeroParam = "#C0392B", Logistic = "#2E86AB", LASSO = "#7B5EA7",
            RF = "#E08214", SVM = "#3D8B8B")
labB <- bm[, .(auc = auc[1], n = n[1]), by = short]
labB[, txt := sprintf("%s - AUC %.3f", short, auc)]
labB_txt <- setNames(labB$txt, as.character(labB$short))
labB_txt <- gsub("ZeroParam", "Zero-parameter (proposed)", labB_txt)
labB_txt <- gsub("LASSO", "LASSO (1-SE)", labB_txt)
labB_txt <- gsub("RF", "Random forest", labB_txt)
labB_txt <- gsub("SVM", "SVM (RBF)", labB_txt)
pB <- ggplot(bm, aes(fpr, tpr, colour = short)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70", linewidth = 0.3) +
  geom_step(linewidth = 0.65) +
  scale_colour_manual(values = PAL_BM, labels = labB_txt, name = NULL) +
  scale_x_continuous(expand = c(0.01, 0), labels = percent_format(accuracy = 1)) +
  scale_y_continuous(expand = c(0.01, 0), labels = percent_format(accuracy = 1)) +
  labs(title = "Cross-platform benchmark",
       subtitle = "GSE46602 -> GSE70770 across platforms, no retraining",
       x = "1 - specificity", y = "Sensitivity") +
  theme_pub() +
  theme(legend.position = c(0.99, 0.03), legend.justification = c(1, 0),
        legend.text = element_text(size = 5.2),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

# ---------------- C: 外部队列（GSE70770 Cambridge）panel vs 单基因 ----------------
exc <- fread(file.path(tab_dir, "T43c_external_roc_curves.csv"))
exc[, model := factor(model, levels = c("PANEL (fixed-weight)", "AMACR", "PCA3", "HOXC6", "HPN"))]
exc[, model := relevel(model, "PANEL (fixed-weight)")]
exc[, model := factor(model, levels = rev(levels(model)))]
PAL5 <- c("PANEL (fixed-weight)" = "#C0392B", AMACR = "#2E86AB", PCA3 = "#E08214",
          HOXC6 = "#7B5EA7", HPN = "#3D8B8B")
lab5 <- exc[, .(auc = auc[1], n = n[1]), by = model]
lab5[, txt := sprintf("%s - AUC %.3f (n=%d)", sub(" \\(fixed-weight\\)", ", fixed-weight", model), auc, n)]
lab5_txt <- setNames(lab5$txt, as.character(lab5$model))
pC <- ggplot(exc, aes(fpr, tpr, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70", linewidth = 0.3) +
  geom_step(linewidth = 0.65) +
  scale_colour_manual(values = PAL5, labels = lab5_txt, name = NULL) +
  scale_x_continuous(expand = c(0.01, 0), labels = percent_format(accuracy = 1)) +
  scale_y_continuous(expand = c(0.01, 0), labels = percent_format(accuracy = 1)) +
  labs(title = "Against markers in current use",
       subtitle = "GSE70770 (Cambridge); DeLong P: AMACR 0.0011, PCA3 0.00064",
       x = "1 - specificity", y = "Sensitivity") +
  theme_pub() +
  theme(legend.position = c(0.99, 0.03), legend.justification = c(1, 0),
        legend.text = element_text(size = 5.2),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

# ---------------- D: 校准（GSE69223） ----------------
cal <- fread(file.path(tab_dir, "T42_calibration.csv"))
cal_p <- cal[!is.na(predicted)]
pD <- ggplot(cal_p, aes(predicted, observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_point(colour = COH[["GSE69223"]], size = 2.0) +
  geom_line(colour = COH[["GSE69223"]], linewidth = 0.6) +
  xlim(0, 1) + ylim(0, 1) +
  labs(title = "Calibration",
       subtitle = "GSE69223; other cohorts separate completely (AUC = 1)",
       x = "Predicted probability of tumour", y = "Observed proportion") +
  theme_pub()

# ---------------- E: DCA（三队列） ----------------
dca <- fread(file.path(tab_dir, "T41_dca_netbenefit.csv"))
dca_l <- melt(dca, id.vars = c("cohort", "threshold"),
              measure.vars = c("nb_model", "nb_treat_all", "nb_treat_none"))
dca_l[, model := fcase(variable == "nb_model", "Four-gene panel",
                       variable == "nb_treat_all", "Treat all",
                       default = "Treat none")]
dca_l[, cohort := factor(cohort, levels = names(COH))]
pE <- ggplot(dca_l, aes(threshold, value, colour = model, linewidth = variable == "nb_model")) +
  geom_line() + facet_wrap(~cohort, nrow = 1) +
  scale_colour_manual(values = c("Four-gene panel" = "#1B3A57",
                                 "Treat all" = "grey55", "Treat none" = "grey80"), name = NULL) +
  scale_linewidth_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.45), guide = "none") +
  coord_cartesian(ylim = c(-0.05, 0.80)) +
  labs(title = "Decision-curve analysis",
       subtitle = "Net benefit of biopsying according to the panel score, by threshold probability",
       x = "Threshold probability", y = "Net benefit") +
  theme_pub(base = 7) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(size = 6.5, face = "bold", family = FONT))

fig <- (pA | pB) / (pC | pD) / pE +
  plot_layout(heights = c(1, 1, 1.08)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))

save_fig(fig, file.path(fig_dir, "Fig6_informatics.png"), 7.2, 7.9)
cat("saved Fig6_informatics (A-E)\n")
