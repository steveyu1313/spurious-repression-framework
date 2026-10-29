# =============================================================================
# 43_fig5_hub.R — Figure 5 关键基因 / PPI 网络
#   A  Top20 hub 棒棒糖图 (degree, 颜色=logFC)
#   B  hub 及其一级邻居的 PPI 子网
#   C  degree vs betweenness (节点重要性双指标)
#   D  hub 基因在单细胞各细胞型的表达热图
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
  library(ragg); library(ggrepel); library(igraph); library(grid)
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
P <- readRDS(file.path(tab_dir, "ppi_objects.rds"))
nd <- as.data.table(P$nodes); g <- P$graph
COL <- c(Up = "#C0392B", Down = "#2E6DA4")

# ---------------- A: top20 hub ----------------
h <- nd[hub == TRUE][order(hub_rank)]
h[, Symbol := factor(Symbol, levels = rev(Symbol))]
# 标签偏移随节点大小线性变化: 小圆点贴边, MYC 这类大圆点让出半径
h[, lab_x := degree + 2.6 + 5.0 * (MCC - min(MCC)) / max(1, diff(range(MCC)))]
pA <- ggplot(h, aes(x = degree, y = Symbol)) +
  geom_segment(aes(x = 0, xend = degree, yend = Symbol, colour = direction),
               linewidth = 0.6) +
  geom_point(aes(size = MCC, colour = direction)) +
  geom_text(aes(x = lab_x, label = sprintf("%.1f", logFC)), hjust = 0,
            size = 2.05, family = FONT, colour = "grey30") +
  scale_colour_manual(values = COL, name = "DEG") +
  scale_size_continuous(range = c(0.8, 3.0), name = "MCC") +
  scale_x_continuous(limits = c(0, 96), breaks = c(0, 20, 40, 60, 80),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(title = "Hub genes ranked by multi-metric consensus",
       subtitle = "bar = degree; label = log2FC", x = "Degree in PPI network", y = NULL) +
  theme_pub() +
  theme(panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey90"),
        axis.text.y = element_text(size = 6.4, face = "italic"),
        legend.position = "right")

# ---------------- B: PPI subnetwork ----------------
hub <- nd[hub == TRUE, Symbol]
nb_all <- unlist(lapply(hub, function(v) {
  i <- which(V(g)$name == v); if (!length(i)) return(character(0))
  as.character(neighbors(g, i)$name)
}))
# 控制节点数: hub + 度>=8 的邻居, 避免网络过密导致标签重叠
keep_deg <- nd[Symbol %in% nb_all & degree >= 8, Symbol]
keep <- unique(c(hub, keep_deg))
sg <- induced_subgraph(g, which(V(g)$name %in% keep))
set.seed(7)
lay <- layout_with_fr(sg, niter = 2000, grid = "nogrid")
xy <- as.data.table(lay); setnames(xy, c("x", "y"))
xy[, Symbol := V(sg)$name]
xy <- merge(xy, nd[, .(Symbol, degree, logFC, direction, hub, hub_rank)], by = "Symbol")
ed <- as.data.table(as_data_frame(sg, what = "edges"))[, .(from, to)]
ed <- merge(ed, xy[, .(name = Symbol, x, y)], by.x = "from", by.y = "name")
ed <- merge(ed, xy[, .(name = Symbol, xend = x, yend = y)], by.x = "to", by.y = "name")

lab <- xy[hub == TRUE][order(hub_rank)][1:min(15, .N)]
pB <- ggplot() +
  geom_segment(data = ed, aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.16, colour = "grey72") +
  geom_point(data = xy[hub == FALSE], aes(x, y, size = degree, fill = logFC),
             shape = 21, stroke = 0.15, colour = "grey35", alpha = 0.85) +
  geom_point(data = xy[hub == TRUE], aes(x, y, size = degree, fill = logFC),
             shape = 21, stroke = 0.55, colour = "black") +
  geom_text_repel(data = lab, aes(x, y, label = Symbol), size = 2.15, family = FONT,
                  fontface = "bold", box.padding = 0.62, point.padding = 0.5,
                  force = 16, force_pull = 0.6, max.iter = 60000,
                  bg.color = "white", bg.r = 0.14,
                  segment.size = 0.22, segment.colour = "grey45",
                  min.segment.length = 0, max.overlaps = Inf) +
  scale_fill_gradient2(low = "#2E6DA4", mid = "grey96", high = "#C0392B",
                       midpoint = 0, name = "log2FC",
                       limits = c(-5.2, 5.2), oob = squish) +
  scale_size_continuous(range = c(0.35, 3.4), name = "Degree") +
  labs(title = "PPI subnetwork of the 20 hub genes",
       subtitle = sprintf("%d nodes, %d edges (STRING v12, score \u2265 400)", vcount(sg), ecount(sg))) +
  coord_fixed() +
  theme_void(base_size = 7, base_family = FONT) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 6.3, colour = "grey25", hjust = 0),
        legend.position = "right", legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6.2), legend.title = element_text(size = 6.5),
        plot.margin = margin(2, 2, 2, 2))

# ---------------- C: hub expression across bulk samples ----------------
Bs <- readRDS(file.path(tab_dir, "bulk_recheck.rds"))
grp <- factor(ifelse(Bs$group == "Tumor", "Tumour", "Benign"), levels = c("Benign", "Tumour"))
ord <- order(grp)
hg <- nd[hub == TRUE][order(hub_rank), Symbol]
sub <- Bs$mat[hg, ord]
z <- t(scale(t(sub))); z[is.na(z)] <- 0
dz <- as.data.table(as.data.frame(z), keep.rownames = "gene")
dl <- melt(dz, id.vars = "gene", variable.name = "sample", value.name = "z")
dl[, `:=`(gene = factor(gene, levels = rev(hg)),
          sample = factor(sample, levels = colnames(sub)))]
nB <- sum(grp == "Benign")
pC <- ggplot(dl, aes(sample, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  annotate("rect", xmin = 0.5, xmax = nB + 0.5, ymin = 0.5, ymax = length(hg) + 0.5,
           fill = NA, colour = "#2E6DA4", linewidth = 0.4) +
  annotate("rect", xmin = nB + 0.5, xmax = length(ord) + 0.5, ymin = 0.5,
           ymax = length(hg) + 0.5, fill = NA, colour = "#C0392B", linewidth = 0.4) +
  annotate("text", x = nB / 2 + 0.5, y = -0.9, label = sprintf("Benign (n=%d)", nB),
           size = 2.2, family = FONT, colour = "#2E6DA4", fontface = "bold") +
  annotate("text", x = nB + (length(ord) - nB) / 2 + 0.5, y = -0.9,
           label = sprintf("Tumour (n=%d)", length(ord) - nB),
           size = 2.2, family = FONT, colour = "#C0392B", fontface = "bold") +
  scale_fill_gradient2(low = "#2E6DA4", mid = "#F7F7F7", high = "#C0392B",
                       midpoint = 0, name = "z-score") +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.2))) +
  coord_cartesian(clip = "off") +
  labs(title = "Hub genes separate benign from tumour tissue",
       subtitle = "row-scaled expression, 14 benign + 36 tumour (GSE46602)",
       x = NULL, y = NULL) +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_text(size = 6.0, face = "italic"),
        plot.margin = margin(3, 3, 9, 2),
        legend.key.height = unit(9, "mm"), legend.key.width = unit(2.4, "mm"),
        legend.title = element_text(size = 6.4))

# ---------------- D: hub x cell type ----------------
sc <- fread(file.path(tab_dir, "T28_hub_scrna_by_celltype.csv"))
ct_cols <- grep("^pct_", names(sc), value = TRUE)
avg_cols <- setdiff(names(sc), c("gene", "is_hub", ct_cols))
sc <- sc[gene %in% nd[hub == TRUE, Symbol]]
m <- as.matrix(sc[, ..avg_cols]); rownames(m) <- sc$gene
z <- t(scale(t(m)))
z[is.na(z)] <- 0
dz <- as.data.table(as.data.frame(z), keep.rownames = "gene")
dl <- melt(dz, id.vars = "gene", variable.name = "CellType", value.name = "z")
dl[, gene := factor(gene, levels = rev(nd[hub == TRUE][order(hub_rank), Symbol]))]
dl <- dl[!is.na(gene)]
ct_order <- c("Luminal", "Basal", "Cycling", "Fibroblast", "Smooth-muscle",
              "Endothelial", "Macrophage", "T-cells", "B-cells", "Mast-cells")
dl[, CellType := factor(CellType, levels = intersect(ct_order, unique(CellType)))]
pD <- ggplot(dl, aes(CellType, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2E6DA4", mid = "#F7F7F7", high = "#C0392B",
                       midpoint = 0, name = "z-score") +
  labs(title = "Hub genes across prostatic cell types",
       subtitle = "mean expression, row-scaled (GSE141445)", x = NULL, y = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 42, hjust = 1, vjust = 1, size = 6.2),
        axis.text.y = element_text(size = 6.0, face = "italic"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right", legend.key.width = unit(2.4, "mm"),
        legend.key.height = unit(9, "mm"),
        legend.title = element_text(size = 6.4),
        plot.margin = margin(3, 2, 3, 2))

# ---------------- E: druggability (DGIdb) ----------------
dr <- fread(file.path(tab_dir, "T36_druggability_DGIdb.csv"))
dmap <- unique(nd[, .(gene = Symbol, direction)])
dr <- merge(dr, dmap, by = "gene", all.x = TRUE)
dr[is.na(direction), direction := "Up"]        # PCA3: panel gene, up-regulated, absent from STRING
dr[, grp := ifelse(in_panel == 1, "Panel", "Hub")]
dr[, lab := as.character(n_drugs)]
dr[, gene := factor(gene, levels = gene[order(n_drugs, gene)])]

pE <- ggplot(dr, aes(n_drugs, gene)) +
  geom_segment(aes(x = 0, xend = n_drugs, yend = gene, colour = direction),
               linewidth = 0.45) +
  geom_point(aes(colour = direction, size = n_sources, shape = grp)) +
  geom_text(aes(x = expm1(log1p(n_drugs) + 0.09), label = lab), hjust = 0, size = 1.95,
            family = FONT, colour = "grey30") +
  scale_x_continuous(trans = "log1p", limits = c(0, 430),
                     breaks = c(0, 1, 5, 20, 50, 200),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_colour_manual(values = c(COL, "Not annotated" = "grey55"), name = NULL) +
  scale_size_continuous(range = c(0.45, 2.1), name = "Sources") +
  scale_shape_manual(values = c(Panel = 17, Hub = 16), name = NULL) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2),
         size = guide_legend(order = 3)) +
  labs(title = "Druggability of the hub and panel genes",
       subtitle = "annotated drugs in DGIdb; only MYC, EZH2 and FASN are both up-regulated and druggable",
       x = "Number of annotated drugs (log scale)", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.0, face = "italic"),
        axis.text.x = element_text(size = 6.6),
        legend.position = "right", legend.box = "vertical",
        legend.key.size = unit(3.0, "mm"),
        plot.margin = margin(3, 3, 3, 2))

design <- "
AAAAAABBBBBB
CCCCCCCDDDDD
EEEEEEEEEEEE
"
p <- pA + pB + pC + pD + pE +
  plot_layout(design = design, heights = c(1.00, 1.28, 1.12)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9.5, face = "bold", family = FONT))

save_fig(p, file.path(fig_dir, "Fig4_hub_PPI.png"), 7.2, 9.25)
ggsave(file.path(tempdir(), "Fig4.pdf"), p, width = 7.2, height = 9.25, device = cairo_pdf)
file.copy(file.path(tempdir(), "Fig4.pdf"), file.path(fig_dir, "Fig4_hub_PPI.pdf"), overwrite = TRUE)
cat("Fig4 saved\n")
