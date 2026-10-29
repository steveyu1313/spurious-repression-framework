# =============================================================================
# 23_fig4_spatial.R — Figure 4 重绘 (GSE245388; 仅 in-tissue spots, 4 唯一切片)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales); library(ragg)
})
FONT <- "Arial"; fig_dir <- "./08_QC_Revision/figures"
theme_pub <- function(base = 6.4) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.line = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks.length = unit(1, "mm"),
          axis.text = element_text(colour = "black", size = base - 0.6),
          axis.title = element_text(colour = "black", size = base),
          plot.title = element_text(size = base, hjust = 0),
          plot.subtitle = element_text(size = base - 1, colour = "grey25", hjust = 0),
          legend.key.size = unit(2.4, "mm"),
          legend.text = element_text(size = base - 1.2),
          legend.title = element_text(size = base - 0.8),
          strip.text = element_text(size = base - 0.4, face = "bold"),
          strip.background = element_blank(),
          panel.spacing = unit(1.5, "mm"),
          plot.margin = margin(1.2, 1.6, 1.2, 1.6))
}
save_fig <- function(p, file, w, h, dpi = 600) {
  tmp <- file.path(tempdir(), basename(file))
  ggsave(tmp, p, width = w, height = h, dpi = dpi, device = ragg::agg_png, bg = "white")
  file.copy(tmp, file, overwrite = TRUE)
}

md <- fread("./08_QC_Revision/tables/X8_spatial_reanalysed.csv")
md[, cluster := factor(sub("^g", "", cluster), levels = as.character(0:4))]
lab5 <- c(`0` = "D0 fibromuscular stroma", `1` = "D1 luminal epithelium",
          `2` = "D2 luminal (MSMB-mid)", `3` = "D3 luminal (MSMB-high)",
          `4` = "D4 tumour epithelium (AMACR-high)")
PAL5 <- c(`0` = "#6E7B8B", `1` = "#E8A33D", `2` = "#8C9E3F", `3` = "#C46BA0", `4` = "#C0504D")
md[, x := array_col]; md[, y := -array_row]
secs <- c("CDZ", "CFX", "HB", "LGL")
lims <- md[, .(x0 = min(x) - 4, x1 = max(x) + 4, y0 = min(y) - 4, y1 = max(y) + 4), by = section]

# ---------------- A: 空间域图 (逐切片) ----------------
map_dom <- function(s) {
  L <- lims[section == s]
  ggplot(md[section == s], aes(x, y, colour = cluster)) +
    geom_point(shape = 16, size = 1.25, stroke = 0) +
    coord_fixed(ratio = 1, xlim = c(L$x0, L$x1), ylim = c(L$y0, L$y1), expand = FALSE) +
    scale_colour_manual(values = PAL5, name = NULL, labels = lab5, drop = FALSE) +
    labs(title = s, x = NULL, y = NULL) +
    theme_pub(6.2) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.title = element_text(size = 6.6, hjust = 0.5, face = "bold"),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5))
}
pA <- wrap_plots(lapply(secs, map_dom), nrow = 1, guides = "collect") +
  plot_annotation(title = "Spatial domains across four prostate cancer sections",
                  subtitle = sprintf("%d in-tissue spots (64\u2013444 per section); five transcriptional domains",
                                     nrow(md))) &
  theme(legend.position = "right", plot.title = element_text(family = FONT, size = 7.6, face = "bold"),
        plot.subtitle = element_text(family = FONT, size = 6.2, colour = "grey25"))

# ---------------- B: 基因空间表达 ----------------
map_gene <- function(g, s) {
  L <- lims[section == s]
  d <- md[section == s]
  ggplot(d, aes(x, y, colour = get(g))) +
    geom_point(shape = 16, size = 1.25, stroke = 0) +
    coord_fixed(ratio = 1, xlim = c(L$x0, L$x1), ylim = c(L$y0, L$y1), expand = FALSE) +
    scale_colour_gradientn(colors = c("grey93", "#FDD49E", "#E34A33", "#67000D"),
                           name = g, limits = glim[[g]], oob = scales::squish) +
    labs(title = if (g == "KLK3") s else "", x = NULL, y = NULL) +
    theme_pub(6.0) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          plot.title = element_text(size = 6.2, face = "bold", hjust = 0.5),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5))
}
glim <- lapply(setNames(c("KLK3", "AMACR"), c("KLK3", "AMACR")), function(g) range(md[[g]]))
grid_g <- list()
for (g in c("KLK3", "AMACR")) for (s in secs) grid_g[[paste(g, s)]] <- map_gene(g, s)
pB <- wrap_plots(grid_g, nrow = 2, byrow = TRUE, guides = "collect") +
  plot_annotation(title = "Spatial expression of KLK3 and AMACR") &
  theme(legend.position = "right",
        plot.title = element_text(family = FONT, size = 7.4, face = "bold"))

# ---------------- C: UMAP ----------------
pC <- ggplot(md[order(cluster)], aes(UMAP1, UMAP2, colour = cluster)) +
  geom_point(size = 0.5, stroke = 0) +
  scale_colour_manual(values = PAL5, name = NULL, labels = lab5, guide = "none") +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Domain embedding",
       subtitle = "Harmony-integrated by section") +
  theme_pub(6.4) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 0.3))

# ---------------- D: domain marker dotplot ----------------
mk <- c("KLK3", "MSMB", "AMACR", "AR", "KRT5", "ACTA2", "PTPRC")
av <- md[, lapply(.SD, mean), by = cluster, .SDcols = mk]
pct <- md[, lapply(.SD, function(v) mean(v > 0)), by = cluster, .SDcols = mk]
d1 <- melt(av, id.vars = "cluster", variable.name = "gene", value.name = "avg")
d2 <- melt(pct, id.vars = "cluster", variable.name = "gene", value.name = "pct")
dd <- merge(d1, d2, by = c("cluster", "gene"))
dd[, gene := factor(gene, levels = mk)]
dd[, z := as.numeric(scale(avg)), by = gene]
dd[, z := pmin(pmax(z, -2), 2)]
pD <- ggplot(dd, aes(gene, cluster)) +
  geom_point(aes(size = pct, colour = z)) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0,
                         name = "Scaled\nmean", breaks = c(-2, 0, 2), labels = c("-2", "0", "2")) +
  scale_size_continuous(range = c(0.2, 2.9), name = "Fraction\nexpressing", limits = c(0, 1),
                        breaks = c(0.25, 0.5, 0.75)) +
  scale_y_discrete(labels = lab5, limits = rev(levels(dd$cluster))) +
  labs(x = NULL, y = NULL, title = "Domain marker expression",
       subtitle = "Size: fraction of spots; colour: scaled mean") +
  theme_pub(6.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 5.6, face = "italic"),
        panel.grid.major = element_line(linewidth = 0.15, colour = "grey93"),
        axis.line = element_blank(), axis.ticks = element_blank())

fig4 <- wrap_elements(pA) / wrap_elements(pB) / (pC | pD) +
  plot_layout(heights = c(0.55, 0.75, 0.9)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(family = FONT, face = "bold", size = 9)))

save_fig(fig4, file.path(fig_dir, "Fig3_spatial.png"), 7.09, 8.4)
cat("Fig4 saved\n")
