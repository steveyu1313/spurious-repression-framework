# =============================================================================
# 44_fig6_pathway.R — Figure 6 通路富集与上下游调控
#   A  GO BP/CC/MF 分面气泡图
#   B  KEGG 通路
#   C  Hallmark GSEA (NES)
#   D  上游转录因子 (ChEA3, 全部 DEG)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
  library(ragg); library(stringr); library(grid)
})
FONT <- "Arial"
fig_dir <- "./08_QC_Revision/figures"; tab_dir <- "./08_QC_Revision/tables"
theme_pub <- function(base = 7) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks.length = unit(1.2, "mm"),
          axis.text = element_text(colour = "black", size = base),
          axis.title = element_text(colour = "black", size = base + 0.5),
          plot.title = element_text(size = base + 1, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = base - 0.8, colour = "grey25", hjust = 0),
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = base - 0.9),
          legend.title = element_text(size = base - 0.4),
          plot.margin = margin(2, 3, 2, 2))
}
save_fig <- function(p, file, w, h, dpi = 600) {
  tmp <- file.path(tempdir(), basename(file))
  ggsave(tmp, p, width = w, height = h, dpi = dpi, device = ragg::agg_png, bg = "white")
  file.copy(tmp, file, overwrite = TRUE)
}
wrap <- function(x, n = 42) str_wrap(x, width = n)

# ---------------- A: GO ----------------
go <- fread(file.path(tab_dir, "T19_GO_enrichment.csv"))
go <- go[gene_set == "All"]
go[, ontology := factor(ontology, levels = c("BP", "CC", "MF"),
                        labels = c("BP", "CC", "MF"))]
top <- go[order(ontology, p.adjust)][, head(.SD, 7), by = ontology]
top[, Description := wrap(tolower(Description), 44)]
top[, Description := factor(Description, levels = rev(unique(Description[order(ontology, p.adjust)])))]
top[, logp := -log10(p.adjust)]
top[, gr := as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))]
pA <- ggplot(top, aes(gr, Description)) +
  geom_point(aes(size = Count, colour = logp)) +
  facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
  scale_colour_gradientn(colours = c("#FDD0A2", "#E6550D", "#7F2704"),
                         name = expression(-log[10]~italic(P)[adj])) +
  scale_size_continuous(range = c(0.7, 3.2), name = "Count") +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.12, 0.18))) +
  labs(title = "GO enrichment (704 DEGs)",
       subtitle = "BP / CC / MF, BH-adjusted; universe = 21,358 genes",
       x = "Gene ratio", y = NULL) +
  theme_pub() +
  theme(strip.text.y = element_text(size = 6, angle = 0, face = "bold", hjust = 0),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.placement = "outside",
        panel.spacing.y = unit(1.8, "mm"),
        axis.text.y = element_text(size = 5.3),
        legend.position = "right",
        legend.box = "vertical")

# ---------------- B: KEGG ----------------
kk <- fread(file.path(tab_dir, "T20_KEGG_enrichment.csv"))
kk[, logp := -log10(p.adjust)]
kk[, Description := wrap(tolower(Description), 40)]
setorder(kk, logp)
kk[, Description := factor(Description, levels = Description)]
pB <- ggplot(kk, aes(logp, Description)) +
  geom_segment(aes(x = 0, xend = logp, yend = Description), linewidth = 0.5, colour = "grey60") +
  geom_point(aes(size = as.numeric(sub("/.*", "", GeneRatio))), colour = "#3B6EA5") +
  scale_size_continuous(range = c(0.9, 3.0), name = "Count") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "KEGG pathways",
       subtitle = sprintf("%d pathways, FDR < 0.05", nrow(kk)),
       x = expression(-log[10]~italic(P)[adj]), y = NULL) +
  theme_pub() + theme(axis.text.y = element_text(size = 5.5))

# ---------------- C: Hallmark ----------------
hl <- fread(file.path(tab_dir, "T21_GSEA_hallmark.csv"))
hl <- hl[padj < 0.05]
sel <- rbind(hl[NES > 0][order(padj)][1:min(9, .N)], hl[NES < 0][order(padj)][1:min(11, .N)])
sel[, name := gsub("^HALLMARK_", "", pathway)]
sel[, name := wrap(str_to_sentence(gsub("_", " ", tolower(name))), 34)]
sel[, dir := ifelse(NES > 0, "Up in tumour", "Down in tumour")]
setorder(sel, NES)
sel[, name := factor(name, levels = name)]
pC <- ggplot(sel, aes(NES, name, fill = dir)) +
  geom_col(width = 0.68) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("Up in tumour" = "#C0392B", "Down in tumour" = "#2E6DA4"), name = NULL) +
  labs(title = "Hallmark gene sets (GSEA)",
       subtitle = sprintf("%d of 50 sets at FDR < 0.05", nrow(hl)),
       x = "Normalised enrichment score", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 5.4),
        legend.position = "none")

# ---------------- D: upstream TF (ChEA3) ----------------
ch <- fread(file.path(tab_dir, "T23_upstream_TF_ChEA3.csv"))
ch <- ch[1:min(12, .N)]
ch[, n_overlap := length(strsplit(Overlapping_Genes, ",")[[1]]), by = Rank]
ch[, Score := as.numeric(Score)]
ch[, TF := gsub("^NKX31$", "NKX3-1", TF)]
ch[, TF := factor(TF, levels = rev(TF))]
ch[, hl := ifelse(Rank <= 3, "top 3", "other")]
pD <- ggplot(ch, aes(Score, TF)) +
  geom_segment(aes(x = 0, xend = Score, yend = TF), linewidth = 0.5, colour = "grey65") +
  geom_point(aes(size = n_overlap, colour = hl)) +
  geom_text(aes(label = n_overlap), hjust = -0.9, size = 1.85, family = FONT, colour = "grey30") +
  scale_colour_manual(values = c("top 3" = "#7B3294", "other" = "grey45"), guide = "none") +
  scale_size_continuous(range = c(1.0, 3.2), name = "Overlap") +
  scale_x_continuous(limits = c(0, max(ch$Score) * 1.32),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Upstream TFs (ChEA3)",
       subtitle = "lower rank = stronger evidence",
       x = "ChEA3 mean rank score", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.4, face = "bold"),
        legend.position = "none")

design <- "
AABB
CCDD
"
p <- pA + pB + pC + pD + plot_layout(design = design, widths = c(1.15, 1), heights = c(1.05, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))
save_fig(p, file.path(fig_dir, "Fig5_pathway_TF.png"), 7.4, 7.4)
ggsave(file.path(tempdir(), "Fig5.pdf"), p, width = 7.4, height = 7.4, device = cairo_pdf)
file.copy(file.path(tempdir(), "Fig5.pdf"), file.path(fig_dir, "Fig5_pathway_TF.pdf"), overwrite = TRUE)
cat("Fig5 saved\n")
