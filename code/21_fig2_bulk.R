# =============================================================================
# 21_fig2_bulk.R — Figure 2 重绘 (GSE46602, 50 例)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
  library(ragg); library(ggrepel); library(grid)
})
FONT <- "Arial"
fig_dir <- "./08_QC_Revision/figures"
theme_pub <- function(base = 7) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks.length = unit(1.2, "mm"),
          axis.text = element_text(colour = "black", size = base),
          axis.title = element_text(colour = "black", size = base + 0.5),
          plot.title = element_text(size = base + 1, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = base - 0.7, colour = "grey25", hjust = 0),
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

B <- readRDS("./08_QC_Revision/tables/bulk_recheck.rds")
mat <- B$mat; grp <- factor(ifelse(B$group == "Tumor", "Tumour", "Benign"), levels = c("Benign", "Tumour"))
res <- B$res
COL <- c(Tumour = "#C0392B", Benign = "#2E6DA4")

# ---------------- A: PCA ----------------
pca <- prcomp(t(mat), scale. = TRUE)
ve <- 100 * summary(pca)$importance[2, 1:2]
pd <- data.table(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Group = grp)
xr <- quantile(pd$PC1, c(0.01, 0.99)); yr <- quantile(pd$PC2, c(0.01, 0.99))
pad <- function(v, k = 0.35) c(v[1] - diff(v) * k, v[2] + diff(v) * k)
out_id <- pd$sample[which.max(abs(pd$PC1))]
pd_plot <- pd[pd$sample != out_id]
n_out <- 1
pA <- ggplot(pd_plot, aes(PC1, PC2, fill = Group)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
  stat_ellipse(aes(colour = Group), level = 0.95, linetype = 2, linewidth = 0.35, show.legend = FALSE) +
  geom_point(shape = 21, size = 1.9, stroke = 0.25, colour = "white") +
  scale_fill_manual(values = COL, name = NULL) +
  scale_colour_manual(values = COL) +
  coord_cartesian() +
  annotate("text", x = pad(xr)[1], y = pad(yr)[2], hjust = 0, vjust = 1,
           label = sprintf("1 sample (%s) off-scale, omitted", out_id), family = FONT, size = 1.8, colour = "grey40") +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
       title = "Transcriptome-wide separation",
       subtitle = sprintf("n = %d (%d tumour, %d benign)", ncol(mat), sum(grp == "Tumour"), sum(grp == "Benign"))) +
  theme_pub() + theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0))

# ---------------- B: Volcano ----------------
v <- copy(res)
v[, sig := fcase(adj.P.Val < 0.05 & logFC > 1, "Up",
                 adj.P.Val < 0.05 & logFC < -1, "Down", default = "NS")]
v[, sig := factor(sig, levels = c("Up", "Down", "NS"))]
lab <- rbind(v[sig == "Up"][order(-logFC)][1:6], v[sig == "Down"][order(logFC)][1:6])
pB <- ggplot(v, aes(logFC, -log10(adj.P.Val))) +
  geom_point(data = v[sig == "NS"], colour = "grey80", size = 0.45, stroke = 0) +
  geom_point(data = v[sig != "NS"], aes(colour = sig), size = 0.65, stroke = 0) +
  geom_vline(xintercept = c(-1, 1), linetype = 3, linewidth = 0.25, colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = 3, linewidth = 0.25, colour = "grey50") +
  geom_text_repel(data = lab, aes(label = Symbol), size = 1.9, family = FONT,
                  segment.size = 0.15, min.segment.length = 0, max.overlaps = 20,
                  box.padding = 0.18, show.legend = FALSE) +
  scale_colour_manual(values = c(Up = "#C0392B", Down = "#2E6DA4"), name = NULL,
                      labels = c(sprintf("Up (%d)", sum(v$sig == "Up")),
                                 sprintf("Down (%d)", sum(v$sig == "Down")))) +

  labs(x = expression(log[2]~"fold change (tumour / benign)"), y = expression(-log[10]~"adjusted"~italic(P)),
       title = "Differential expression",
       subtitle = sprintf("%d DEGs (260 up, 444 down); |log2FC| > 1, adjusted P < 0.05", sum(v$sig != "NS"))) +
  theme_pub() + theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1))

# ---------------- C: Heatmap ----------------
top <- res[adj.P.Val < 0.05 & abs(logFC) > 1][order(-abs(logFC))][1:30, Symbol]
hm <- mat[top, order(grp)]
z <- t(scale(t(hm))); z[z > 2] <- 2; z[z < -2] <- -2
hd <- as.data.table(as.table(z)); setnames(hd, c("Gene", "Sample", "Z"))
hd[, Gene := factor(Gene, levels = rev(top))]
ord <- data.table(Sample = colnames(hm), Group = grp[order(grp)])
hd[, Sample := factor(Sample, levels = colnames(hm))]
pC <- ggplot(hd, aes(Sample, Gene, fill = Z)) +
  geom_tile() +
  geom_vline(xintercept = sum(grp[order(grp)] == "Benign") + 0.5, colour = "white", linewidth = 1.2) +
  scale_fill_gradient2(low = "#2E6DA4", mid = "white", high = "#C0392B", midpoint = 0,
                       limits = c(-2, 2), name = "Z-score",
                       breaks = c(-2, 0, 2), labels = c("\u2264-2", "0", "\u22652")) +
  scale_x_discrete(expand = c(0, 0), labels = NULL, breaks = NULL) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Samples (left: benign, right: tumour)", y = NULL,
       title = "Top 30 differentially expressed genes",
       subtitle = "Row-wise Z-score of quantile-normalised expression") +
  theme_pub(6.2) +
  theme(axis.text.y = element_text(size = 5.4, face = "italic"),
        axis.ticks.y = element_blank(), axis.line = element_blank(),
        legend.position = "right", legend.key.height = unit(6, "mm"),
        legend.key.width = unit(2, "mm"))

# ---------------- D: GO dotplot ----------------
go <- fread("./02_Bulk_RNAseq/results/tables/GSE46602_GO_BP_enrichment.csv")
go[, gr := as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))]
g12 <- go[order(p.adjust)][1:12]
g12[, Description := factor(Description, levels = rev(Description))]
g12[, lab := sprintf("%d/%d", Count, as.integer(sub(".*/", "", GeneRatio)))]
pD <- ggplot(g12, aes(gr, Description)) +
  geom_segment(aes(x = 0, xend = gr, yend = Description), colour = "grey85", linewidth = 0.3) +
  geom_point(aes(size = Count, colour = -log10(p.adjust))) +
  geom_text(aes(label = lab), size = 1.75, family = FONT, hjust = -0.45, vjust = 0.35, colour = "grey30") +
  scale_colour_gradient(low = "#F5C6B0", high = "#9B2226", name = expression(-log[10]~"FDR")) +
  scale_size_continuous(range = c(0.7, 3.0), name = "Genes") +
  scale_x_continuous(limits = c(0, 0.115), breaks = seq(0, 0.10, 0.05)) +
  labs(x = "Gene ratio", y = NULL, title = "GO biological processes",
       subtitle = sprintf("Top 12 of %d enriched terms", nrow(go))) +
  theme_pub(6.2) +
  theme(axis.text.y = element_text(size = 5.6))

# ---------------- E: marker boxplots ----------------
mk <- c("KLK3", "AR", "AMACR", "ERG", "GSTP1", "NKX3-1")
bd <- rbindlist(lapply(mk, function(g) {
  data.table(Gene = g, Expr = mat[g, ], Group = grp)
}))
bd[, Gene := factor(Gene, levels = mk)]
pv <- bd[, .(p = wilcox.test(Expr ~ Group)$p.value), by = Gene]
pv[, lab := sprintf("P = %s", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3)))]
pE <- ggplot(bd, aes(Group, Expr, fill = Group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3,
               colour = "grey25", alpha = 0.85) +
  geom_jitter(width = 0.12, size = 0.35, stroke = 0, colour = "grey30", alpha = 0.55) +
  geom_text(data = pv, aes(x = 1.5, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.4, size = 1.85, family = FONT) +
  facet_wrap(~Gene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = COL, guide = "none") +
  scale_x_discrete(labels = c("Benign" = "B", "Tumour" = "T")) +
  labs(x = NULL, y = "Normalised expression (log2)", title = "Marker genes",
       subtitle = "Wilcoxon rank-sum test, benign (B) vs tumour (T)") +
  theme_pub(6.2) +
  theme(strip.text = element_text(size = 6, face = "italic"),
        axis.text.x = element_text(size = 5.8))

fig2 <- (pA | pB) / (pC | pD) / pE +
  plot_layout(heights = c(1, 1.55, 0.72)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(family = FONT, face = "bold", size = 9)))

save_fig(fig2, file.path(fig_dir, "Fig1_transcriptome.png"), 7.09, 8.6)
cat("Fig2 saved\n")
