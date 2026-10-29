# =============================================================================
# 22_fig3_scrna.R — Figure 3 重绘 (GSE141445, 36,424 cells / 13 tumours)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales); library(ragg)
})
FONT <- "Arial"; fig_dir <- "./08_QC_Revision/figures"
theme_pub <- function(base = 7) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks.length = unit(1.1, "mm"),
          axis.text = element_text(colour = "black", size = base),
          axis.title = element_text(colour = "black", size = base + 0.5),
          plot.title = element_text(size = base + 1, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = base - 0.8, colour = "grey25", hjust = 0),
          legend.key.size = unit(2.6, "mm"),
          legend.text = element_text(size = base - 1.2),
          legend.title = element_text(size = base - 0.8),
          plot.margin = margin(2, 3, 2, 2))
}
save_fig <- function(p, file, w, h, dpi = 600) {
  tmp <- file.path(tempdir(), basename(file))
  ggsave(tmp, p, width = w, height = h, dpi = dpi, device = ragg::agg_png, bg = "white")
  file.copy(tmp, file, overwrite = TRUE)
}

um  <- fread("./08_QC_Revision/tables/X1_scrna_umap.csv")
ex  <- fread("./08_QC_Revision/tables/X2_scrna_gene_expr.csv")
avg <- fread("./08_QC_Revision/tables/X3_scrna_marker_avg.csv")
pct <- fread("./08_QC_Revision/tables/X4_scrna_marker_pct.csv")

ct_levels <- c("Luminal", "Basal", "T_cells", "B_cells", "Macrophage", "Endothelial",
               "Fibroblast", "Smooth_muscle", "Mast_cells", "Cycling")
ct_lab <- c(Luminal = "Luminal", Basal = "Basal", T_cells = "T cell", B_cells = "B cell",
            Macrophage = "Macrophage", Endothelial = "Endothelial", Fibroblast = "Fibroblast",
            Smooth_muscle = "Smooth muscle", Mast_cells = "Mast cell", Cycling = "Cycling")
um[, celltype_marker := factor(celltype_marker, levels = ct_levels)]
PAL10 <- c(Luminal = "#C0504D", Basal = "#E8A33D", T_cells = "#4F81BD", B_cells = "#2E8B8B",
           Macrophage = "#7B5EA7", Endothelial = "#C46BA0", Fibroblast = "#8C9E3F",
           Smooth_muscle = "#6E7B8B", Mast_cells = "#D4A017", Cycling = "#3B7A57")

# ---------------- A: UMAP ----------------
um[, ctn := as.integer(celltype_marker)]
setorder(um, -ctn)
cen <- um[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by = celltype_marker]
pA <- ggplot(um, aes(UMAP1, UMAP2, colour = celltype_marker)) +
  geom_point(size = 0.22, stroke = 0, alpha = 0.75) +
  geom_label(data = cen, aes(label = ct_lab[as.character(celltype_marker)]), family = FONT,
             size = 1.95, label.size = 0, label.padding = unit(0.5, "mm"),
             fill = alpha("white", 0.7), colour = "black") +
  scale_colour_manual(values = PAL10, name = NULL, labels = ct_lab) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1), ncol = 1)) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Cellular composition of the tumour microenvironment",
       subtitle = sprintf("%s cells, %d tumours, 22 clusters", format(nrow(um), big.mark = ","),
                          uniqueN(um$gse_sample))) +
  theme_pub(6.6) +
  theme(legend.position = "right", legend.key.height = unit(3.2, "mm"),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), panel.border = element_rect(fill = NA, linewidth = 0.35))

# ---------------- B: 组成堆叠图 ----------------
comp <- um[, .N, by = .(gse_sample, celltype_marker)]
comp[, frac := N / sum(N), by = gse_sample]
ord <- comp[celltype_marker == "Luminal"][order(-frac), gse_sample]
comp[, gse_sample := factor(gse_sample, levels = ord)]
pB <- ggplot(comp, aes(gse_sample, frac, fill = celltype_marker)) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.1) +
  scale_fill_manual(values = PAL10, name = NULL, labels = ct_lab) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Fraction of cells", title = "Cell-type composition per tumour",
       subtitle = "Tumours ordered by luminal fraction") +
  theme_pub(6.6) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, size = 5.6),
        axis.ticks.x = element_line(linewidth = 0.25))

# ---------------- C: marker dotplot ----------------
setnames(avg, 1, "gene"); cn <- setdiff(names(avg), "gene"); setnames(avg, cn, gsub("-", "_", cn))
d1 <- melt(avg, id.vars = "gene", variable.name = "celltype", value.name = "avg")
d2 <- melt(pct, id.vars = "gene", variable.name = "gene_name", value.name = "pct")
setnames(d2, c("gene", "gene_name"), c("celltype", "gene"))
dot <- merge(d1, d2, by = c("gene", "celltype"))
order_g <- c("KLK3","KLK2","ACPP","AR","MSMB", "KRT5","KRT15","TP63", "PECAM1","VWF",
             "CD3D","CD3E", "CD68","C1QC", "MS4A1","CD79A", "COL1A1","DCN",
             "ACTA2","MYH11", "TPSAB1","CPA3", "MKI67","TOP2A")
dot <- dot[gene %in% order_g]
dot[, gene := factor(gene, levels = order_g)]
dot[, celltype := factor(celltype, levels = ct_levels)]
dot[, z := as.numeric(scale(avg)), by = gene]
dot[, z := pmin(pmax(z, -2), 2)]
pC <- ggplot(dot, aes(gene, celltype)) +
  geom_point(aes(size = pct, colour = z)) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0,
                         name = "Scaled mean", breaks = c(-2, 0, 2), labels = c("-2", "0", "2")) +
  scale_size_continuous(range = c(0.15, 3.0), name = "Fraction\nexpressing",
                        breaks = c(0.25, 0.5, 0.75), limits = c(0, 1)) +
  scale_y_discrete(labels = ct_lab, limits = rev(ct_levels)) +
  labs(x = NULL, y = NULL, title = "Canonical lineage markers",
       subtitle = "Size: fraction of cells expressing; colour: gene-wise scaled mean expression") +
  theme_pub(6.4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 5.6, face = "italic"),
        panel.grid.major = element_line(linewidth = 0.18, colour = "grey93"),
        axis.line = element_blank(), axis.ticks = element_blank())

# ---------------- D: feature plots (single facet plot) ----------------
fd <- rbindlist(lapply(c("KLK3", "AMACR", "KRT5", "PTPRC"), function(g) {
  d <- merge(um[, .(barcode, UMAP1, UMAP2)], ex[, .(barcode, v = get(g))], by = "barcode")
  d[, Gene := factor(g, levels = c("KLK3", "AMACR", "KRT5", "PTPRC"))]
}))
setorder(fd, v)
pD <- ggplot(fd, aes(UMAP1, UMAP2, colour = v)) +
  geom_point(size = 0.16, stroke = 0) +
  facet_wrap(~Gene, nrow = 1) +
  scale_colour_gradientn(colors = c("grey90", "#FDD49E", "#E34A33", "#67000D"),
                         name = "Expression",
                         guide = guide_colourbar(barheight = unit(8, "mm"), barwidth = unit(1.6, "mm"))) +
  labs(x = NULL, y = NULL, title = "Marker gene expression") +
  theme_pub(6.2) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), panel.border = element_rect(fill = NA, linewidth = 0.3),
        strip.text = element_text(size = 6.4, face = "italic"),
        strip.background = element_blank(),
        legend.position = "right", legend.title = element_text(size = 5.4),
        legend.text = element_text(size = 5),
        panel.spacing = unit(1.2, "mm"))

fig3 <- pA / pB / pC / pD +
  plot_layout(heights = c(1.15, 0.95, 1.05, 0.7)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(family = FONT, face = "bold", size = 9)))

save_fig(fig3, file.path(fig_dir, "Fig2_singlecell.png"), 7.09, 9.4)
cat("Fig3 saved\n")
