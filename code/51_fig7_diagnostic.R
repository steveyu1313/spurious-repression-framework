# =============================================================================
# 51_fig7_diagnostic.R — Figure 7 三队列验证的诊断/复发分层标志物
#   A  各队列 DEG 数与跨队列共识集（76 个）
#   B  共识 4 基因 panel 的 ROC（发现 → 两独立队列验证）
#   C  同一批基因在三队列中的 logFC 一致性（含基底角蛋白的 bulk 稀释）
#   D  BCR 单因素 Cox（探索性，n = 36）
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
COL <- c(Up = "#C0392B", Down = "#2E6DA4")

# ---------------- A: DEG 数与共识集 ----------------
d <- fread(file.path(tab_dir, "T30_cohort_DEG_all.csv"))
d[, sig := adj.P.Val < 0.05 & abs(logFC) > 1]
cnt <- d[sig == TRUE, .(n = .N,
                        up = sum(logFC > 0), down = sum(logFC < 0)), by = cohort]
cons <- fread(file.path(tab_dir, "T31_crosscohort_intersection.csv"))
cnt <- rbind(cnt, data.table(cohort = "Consensus\n(3 cohorts)", n = sum(cons$consistent),
                             up = sum(cons$up_all), down = sum(cons$down_all)))
cnt[, cohort := factor(cohort, levels = c("GSE46602", "GSE69223", "GSE55945", "Consensus\n(3 cohorts)"))]
cl <- melt(cnt, id.vars = c("cohort", "n"), measure.vars = c("up", "down"),
           variable.name = "dir", value.name = "k")
cl[, dir := factor(ifelse(dir == "up", "Up", "Down"), levels = c("Down", "Up"))]
pA <- ggplot(cl, aes(cohort, k, fill = dir)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.2) +
  geom_text(data = cnt, aes(cohort, n + 28, label = n), inherit.aes = FALSE,
            size = 2.1, family = FONT, fontface = "bold") +
  scale_fill_manual(values = COL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = "Differentially expressed genes",
       subtitle = "|log2FC| > 1, BH FDR < 0.05; 21,355 genes tested in all three",
       x = NULL, y = "Genes") +
  theme_pub() +
  theme(legend.position = c(0.14, 0.86), legend.direction = "horizontal",
        axis.text.x = element_text(size = 5.8, lineheight = 0.9))

# ---------------- B: ROC ----------------
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
pB <- ggplot(rocs, aes(fpr, tpr, colour = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70", linewidth = 0.3) +
  geom_step(linewidth = 0.75) +
  scale_colour_manual(values = COH, labels = lab_txt, name = NULL) +
  scale_x_continuous(expand = c(0.01, 0), labels = percent_format(accuracy = 1)) +
  scale_y_continuous(expand = c(0.01, 0), labels = percent_format(accuracy = 1)) +
  labs(title = "Diagnostic performance",
       subtitle = "Fixed-weight 4-gene panel",
       x = "1 - specificity", y = "Sensitivity") +
  theme_pub() +
  theme(legend.position = c(0.99, 0.03), legend.justification = c(1, 0),
        legend.text = element_text(size = 5.2),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

# ---------------- C: logFC 一致性 ----------------
gl <- c("AMACR", "PCA3", "HOXC6", "HPN", "KRT15", "TP63")
fc <- cons[symbol %in% gl, .(symbol,
                             GSE46602 = logFC_GSE46602,
                             GSE69223 = logFC_GSE69223,
                             GSE55945 = logFC_GSE55945)]
fc <- melt(fc, id.vars = "symbol", variable.name = "cohort", value.name = "logFC")
fc[, cohort := factor(cohort, levels = names(COH))]
fc[, symbol := factor(symbol, levels = rev(gl))]
pC <- ggplot(fc, aes(logFC, symbol)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_vline(xintercept = c(-1, 1), linewidth = 0.25, linetype = 3, colour = "grey65") +
  geom_line(aes(group = symbol), colour = "grey80", linewidth = 0.35) +
  geom_point(aes(colour = cohort), size = 1.7) +
  scale_colour_manual(values = COH, name = NULL) +
  scale_x_continuous(limits = c(-6, 5.2), breaks = seq(-6, 4, 2)) +
  labs(title = "Reproducibility across cohorts",
       subtitle = "Basal keratins collapse in bulk tissue; luminal markers hold",
       x = expression(log[2]~"fold change (tumour vs benign)"), y = NULL) +
  theme_pub() +
  theme(legend.position = c(0.80, 0.22), legend.direction = "vertical",
        legend.text = element_text(size = 5.2),
        axis.text.y = element_text(size = 5.8, face = "italic"))

# ---------------- D: BCR forest ----------------
cx <- fread(file.path(tab_dir, "T34_BCR_cox.csv"))
cx <- cx[gene %in% c("PCA3", "NEK5", "AMACR", "KRT15", "TP63", "HPN", "PANEL_4gene")]
cx[, gene := factor(gene, levels = rev(c("PCA3", "NEK5", "AMACR", "KRT15", "TP63", "HPN", "PANEL_4gene")))]
cx[, gene := factor(gsub("PANEL_4gene", "Panel score", gene),
                    levels = rev(c("PCA3", "NEK5", "AMACR", "KRT15", "TP63", "HPN", "Panel score")))]
cx[, sig := ifelse(p < 0.05, "P < 0.05", "n.s.")]
pD <- ggplot(cx, aes(HR, gene, colour = sig)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.22, linewidth = 0.45) +
  geom_point(size = 1.8) +
  geom_text(aes(x = 4.6, label = sprintf("P=%.3f", p)), hjust = 0, size = 1.85,
            family = FONT, colour = "grey30") +
  scale_x_log10(limits = c(0.12, 9), breaks = c(0.25, 0.5, 1, 2, 4),
                labels = c("0.25", "0.5", "1", "2", "4")) +
  scale_colour_manual(values = c("P < 0.05" = "#C0392B", "n.s." = "grey50"), guide = "none") +
  labs(title = "Biochemical recurrence (exploratory)",
       subtitle = "GSE46602, 36 tumours / 22 events; grey = not significant",
       x = "Hazard ratio (log scale)", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 5.6))

# ---------------- E: 外部队列（GSE70770 Cambridge）ROC ----------------
exc <- fread(file.path(tab_dir, "T43c_external_roc_curves.csv"))
exc[, model := factor(model, levels = c("PANEL (fixed-weight)", "AMACR", "PCA3", "HOXC6", "HPN"))]
exc[, model := relevel(model, "PANEL (fixed-weight)")]
exc[, model := factor(model, levels = rev(levels(model)))]  # 图例自上而下
PAL5 <- c("PANEL (fixed-weight)" = "#C0392B", AMACR = "#2E86AB", PCA3 = "#E08214",
          HOXC6 = "#7B5EA7", HPN = "#3D8B8B")
lab5 <- exc[, .(auc = auc[1], n = n[1]), by = model]
lab5[, txt := sprintf("%s - AUC %.3f (n=%d)", sub(" \\(fixed-weight\\)", ", fixed-weight", model), auc, n)]
lab5_txt <- setNames(lab5$txt, as.character(lab5$model))
pE <- ggplot(exc, aes(fpr, tpr, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70", linewidth = 0.3) +
  geom_step(linewidth = 0.65) +
  scale_colour_manual(values = PAL5, labels = lab5_txt, name = NULL) +
  scale_x_continuous(expand = c(0.01, 0), labels = percent_format(accuracy = 1)) +
  scale_y_continuous(expand = c(0.01, 0), labels = percent_format(accuracy = 1)) +
  labs(title = "External validation on an independent platform",
       subtitle = "GSE70770 (Cambridge, Illumina HumanHT-12); 125 tumours vs 74 benign; paired DeLong: panel vs AMACR P = 0.0011, vs PCA3 P = 0.00064",
       x = "1 - specificity", y = "Sensitivity") +
  theme_pub() +
  theme(legend.position = c(0.99, 0.03), legend.justification = c(1, 0),
        legend.text = element_text(size = 5.2),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

fig <- (pA | pB) / (pC | pD) / pE +
  plot_layout(heights = c(1, 1.02, 1.25)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))

save_fig(fig, file.path(fig_dir, "Fig6_validation.png"), 7.2, 8.1)
cat("saved Fig6 (with external panel E)\n")
