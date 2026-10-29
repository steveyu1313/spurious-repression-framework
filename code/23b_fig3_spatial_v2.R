# =============================================================================
# 23b_fig3_spatial_v2.R — Figure 3 重绘 v2
#   A/B 换成等大正方形拼版: H&E 组织背景(漂白玻片灰底) + 1 mm 标尺,
#   四个唯一切片 (CDZ/CFX/HB/LGL) 同尺寸; A = 五个转录域, B = KLK3/AMACR。
#   C/D (UMAP, dotplot) 沿用 23_fig4_spatial.R。
#   数据: X8_spatial_reanalysed.csv (873 spots) + Visium hires 图像。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
  library(ragg); library(png); library(jsonlite)
})
FONT <- "Arial"; fig_dir <- "./08_QC_Revision/figures"
base_dir <- "04_Spatial/data/Visium"
SECS4 <- c("CDZ", "CFX", "HB", "LGL")
SPOT_UM <- 55; PT_MM <- 0.50

theme_pub <- function(base = 6.4) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.line = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks.length = unit(1, "mm"),
          axis.text = element_text(colour = "black", size = base - 0.6),
          axis.title = element_text(colour = "black", size = base),
          plot.title = element_text(size = base, hjust = 0, face = "bold"),
          legend.key.size = unit(2.4, "mm"),
          legend.text = element_text(size = base - 1.2),
          legend.title = element_text(size = base - 0.8),
          panel.spacing = unit(1.5, "mm"),
          plot.margin = margin(1.2, 1.6, 1.2, 1.6))
}
lab5 <- c(`0` = "D0 fibromuscular stroma", `1` = "D1 luminal epithelium",
          `2` = "D2 luminal (MSMB-mid)", `3` = "D3 luminal (MSMB-high)",
          `4` = "D4 tumour epithelium (AMACR-high)")
PAL5 <- c(`0` = "#6E7B8B", `1` = "#E8A33D", `2` = "#8C9E3F", `3` = "#C46BA0", `4` = "#C0504D")

md <- fread("./08_QC_Revision/tables/X8_spatial_reanalysed.csv")
md[, bc := sub("_[0-9]+$", "", barcode)]
md[, cluster := factor(cluster, levels = as.character(0:4))]

# ---------- 读取 hires 图像与像素坐标 ----------
read_section <- function(s) {
  d  <- file.path(base_dir, s, "outs", "spatial")
  sf <- fromJSON(file.path(d, "scalefactors_json.json"))
  pos <- fread(file.path(d, "tissue_positions_list.csv"), header = FALSE,
               col.names = c("barcode", "in_tissue", "array_row", "array_col",
                             "pxl_row", "pxl_col"))
  h <- sf$tissue_hires_scalef
  pos[, `:=`(x = pxl_col * h, y = pxl_row * h)]
  img <- readPNG(file.path(d, "tissue_hires_image.png"))
  if (dim(img)[3] == 4) {                    # 压平 alpha, 防止导出棋盘格
    a <- img[, , 4]
    for (k in 1:3) img[, , k] <- img[, , k] * a + (1 - a)
    img <- img[, , 1:3, drop = FALSE]
  }
  list(pos = pos[in_tissue == 1], img = img,
       um_px = SPOT_UM / (sf$spot_diameter_fullres * h),
       img_w = dim(img)[2], img_h = dim(img)[1])
}
secs <- lapply(SECS4, read_section); names(secs) <- SECS4

for (s in SECS4) {
  p <- secs[[s]]$pos
  w <- max(p$x) - min(p$x); hh <- max(p$y) - min(p$y)
  px_pad <- w * 0.06; py_pad <- hh * 0.06
  x0 <- max(0, min(p$x) - px_pad); x1 <- min(secs[[s]]$img_w, max(p$x) + px_pad)
  y0 <- max(0, min(p$y) - py_pad); y1 <- min(secs[[s]]$img_h, max(p$y) + py_pad)
  S <- max(x1 - x0, y1 - y0)
  secs[[s]]$crop <- c(x0 = x0, x1 = x1, y0 = y0, y1 = y1, S = S)
}

# ---------- 等大正方形细胞面板 ----------
theme_cell <- theme_void(base_size = 7, base_family = FONT) +
  theme(plot.margin = margin(1.2, 1.2, 1.2, 1.2),
        legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "grey85", linewidth = 0.3),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(size = 6.0, face = "italic", hjust = 0.5,
                                  colour = "grey15", family = FONT))

make_cell <- function(s, mode = c("cluster", "gene"), gene = NULL, lim = NULL,
                      show_sample = TRUE, title = NULL) {
  mode <- match.arg(mode)
  o <- secs[[s]]; b <- o$crop
  ri <- max(1, floor(b["y0"])):min(o$img_h, floor(b["y1"]))
  ci <- max(1, floor(b["x0"])):min(o$img_w, floor(b["x1"]))
  img <- o$img[ri, ci, , drop = FALSE]
  # H&E 提亮: 压暗部 ~48%, 组织结构保留但不再压过表达点
  img <- 1 - (1 - img) * 0.52
  # 载玻片灰底漂白: 中性灰像素置白, 组织(粉/紫)保留
  r <- img[, , 1]; gg <- img[, , 2]; bb <- img[, , 3]
  mx <- pmax(r, gg, bb); mn <- pmin(r, gg, bb)
  bg <- (mx - mn) < 0.055
  for (k in 1:3) { ch <- img[, , k]; ch[bg] <- 1; img[, , k] <- ch }
  w <- length(ci); h <- length(ri); S <- b["S"]
  xo <- (S - w) / 2; yo <- (S - h) / 2

  d <- merge(md[section == s, .(bc, cluster, KLK3, AMACR)],
             o$pos[, .(barcode, x, y)], by.x = "bc", by.y = "barcode")
  d[, `:=`(px = x - ci[1] + 0.5 + xo, py = (ri[length(ri)] - y) + 0.5 + yo)]

  g <- ggplot() +
    annotation_raster(img, xmin = xo, xmax = xo + w, ymin = yo, ymax = yo + h,
                      interpolate = TRUE)
  if (mode == "cluster") {
    g <- g + geom_point(data = d, aes(px, py, colour = cluster),
                        size = PT_MM, shape = 16, stroke = 0, alpha = 0.92) +
      scale_colour_manual(values = PAL5, name = NULL, labels = lab5, drop = FALSE,
                          guide = guide_legend(override.aes = list(size = 2.1, alpha = 1)))
  } else {
    g <- g + geom_point(data = d, aes(px, py, colour = d[[gene]]),
                        size = PT_MM, shape = 16, stroke = 0, alpha = 0.92) +
      scale_colour_gradientn(colors = c("grey93", "#FDD49E", "#E34A33", "#67000D"),
                             name = gene, limits = lim, oob = scales::squish,
                             guide = guide_colourbar(barheight = unit(16, "mm"),
                                                     barwidth = unit(2.2, "mm")))
  }
  bar <- 1000 / o$um_px
  g <- g +
    annotate("segment", x = S * 0.965 - bar, xend = S * 0.965, y = S * 0.055,
             yend = S * 0.055, linewidth = 0.45, colour = "grey15") +
    annotate("text", x = S * 0.965 - bar / 2, y = S * 0.095, label = "1 mm",
             size = 1.85, family = FONT, colour = "grey15", vjust = 0)
  if (show_sample)
    g <- g + annotate("text", x = S * 0.03, y = S * 0.975, label = s,
                      hjust = 0, vjust = 1, size = 2.5, family = FONT,
                      fontface = "bold", colour = "black")
  if (!is.null(title))
    g <- g + labs(title = title)
  g + coord_fixed(xlim = c(0, S), ylim = c(0, S), expand = FALSE, clip = "off") +
    theme_cell
}

glim <- lapply(c("KLK3", "AMACR"), function(g) range(md[[g]]))
names(glim) <- c("KLK3", "AMACR")

# A: 转录域 (1 行 4 列, 图例收集到右侧)
pA <- wrap_plots(lapply(SECS4, make_cell, mode = "cluster"), nrow = 1,
                 guides = "collect") &
  theme(legend.position = "right")

# B: KLK3 / AMACR (2 行 4 列; 基因名标在每行第一格上方)
cellsB <- list()
for (g in c("KLK3", "AMACR"))
  for (i in seq_along(SECS4)) {
    s <- SECS4[i]
    cellsB[[paste(g, s)]] <- make_cell(s, mode = "gene", gene = g, lim = glim[[g]],
                                       show_sample = FALSE,
                                       title = if (i == 1) g else NULL)
  }
pB <- wrap_plots(cellsB, nrow = 2, byrow = TRUE, guides = "collect") &
  theme(legend.position = "right")

# ---------- C: UMAP (沿用 23) ----------
pC <- ggplot(md[order(cluster)], aes(UMAP1, UMAP2, colour = cluster)) +
  geom_point(size = 0.5, stroke = 0) +
  scale_colour_manual(values = PAL5, name = NULL, labels = lab5, guide = "none") +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Domain embedding",
       subtitle = "Harmony-integrated by section") +
  theme_pub(6.4) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 0.3))

# ---------- D: domain marker dotplot (沿用 23) ----------
mkn <- c("KLK3", "MSMB", "AMACR", "AR", "KRT5", "ACTA2", "PTPRC")
av  <- md[, lapply(.SD, mean), by = cluster, .SDcols = mkn]
pct <- md[, lapply(.SD, function(v) mean(v > 0)), by = cluster, .SDcols = mkn]
d1 <- melt(av,  id.vars = "cluster", variable.name = "gene", value.name = "avg")
d2 <- melt(pct, id.vars = "cluster", variable.name = "gene", value.name = "pct")
dd <- merge(d1, d2, by = c("cluster", "gene"))
dd[, gene := factor(gene, levels = mkn)]
dd[, z := as.numeric(scale(avg)), by = gene]
dd[, z := pmin(pmax(z, -2), 2)]
pD <- ggplot(dd, aes(gene, cluster)) +
  geom_point(aes(size = pct, colour = z)) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0,
                         name = "Scaled\nmean", breaks = c(-2, 0, 2),
                         labels = c("-2", "0", "2")) +
  scale_size_continuous(range = c(0.2, 2.9), name = "Fraction\nexpressing",
                        limits = c(0, 1), breaks = c(0.25, 0.5, 0.75)) +
  scale_y_discrete(labels = lab5, limits = rev(levels(dd$cluster))) +
  labs(x = NULL, y = NULL, title = "Domain marker expression",
       subtitle = "Size: fraction of spots; colour: scaled mean") +
  theme_pub(6.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 5.6,
                                   face = "italic"),
        panel.grid.major = element_line(linewidth = 0.15, colour = "grey93"),
        axis.line = element_blank(), axis.ticks = element_blank())

fig3 <- wrap_elements(pA) / wrap_elements(pB) / (pC | pD) +
  plot_layout(heights = c(0.60, 1.18, 0.92)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(family = FONT,
                                                        face = "bold", size = 9)))

dir.create(file.path(fig_dir, "archive"), showWarnings = FALSE)
if (file.exists(file.path(fig_dir, "Fig3_spatial.png")))
  file.copy(file.path(fig_dir, "Fig3_spatial.png"),
            file.path(fig_dir, "archive", "Fig3_spatial_spotonly.png"), overwrite = TRUE)

tmp <- file.path(tempdir(), "Fig3_spatial.png")
ggsave(tmp, fig3, width = 7.09, height = 8.0, dpi = 600, device = ragg::agg_png,
       bg = "white")
file.copy(tmp, file.path(fig_dir, "Fig3_spatial.png"), overwrite = TRUE)
cat("Fig3 v2 saved\n")
