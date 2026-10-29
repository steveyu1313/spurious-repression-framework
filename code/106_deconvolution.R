#!/usr/bin/env Rscript
# =============================================================
# 106_deconvolution.R  (W7: established deconvolution tools)
# 目的: 用现成工具量化 bulk 队列的区室构成变化, 检验
#   "down-regulated hubs 反映 compartment loss" 的论断,
#   并确认 panel 基因对比并非组成效应。
# 工具: MCP-counter (Becht et al., Genome Biol 2016);
#       EPIC (Racle et al., eLife 2017, 输入 TPM)
# 队列: GSE46602 (发现), GSE70770 Cambridge (外部 Illumina),
#       GSE104131 (RNA-seq, 配对)
# 输出: tables/T46k_mcp_scores.csv / T46k_mcp_stats.csv (Table S21)
#       tables/T46k_epic_scores.csv / T46k_epic_stats.csv
#       figures/FigS5_deconvolution.pdf / .jpg
# =============================================================
suppressPackageStartupMessages({
  library(data.table); library(limma); library(ggplot2); library(patchwork)
})
root    <- "."
tab_dir <- file.path(root, "tables")
fig_dir <- file.path(root, "figures")
ext_d   <- file.path(root, "data/external")
source(file.path(ext_d, "MCPcounter.R"))

get_field <- function(lines, key) {
  l <- grep(paste0("^", key, "\t"), lines, value = TRUE)
  if (!length(l)) return(NULL)
  gsub('^"|"$', "", strsplit(l[1], "\t")[[1]][-1])
}
read_sm <- function(path) {
  con <- gzfile(path, "rt"); lines <- readLines(con, warn = FALSE); close(con)
  i1 <- grep("^!series_matrix_table_begin", lines)
  i2 <- grep("^!series_matrix_table_end", lines)
  expr <- as.matrix(read.delim(text = lines[(i1 + 1):(i2 - 1)], header = TRUE,
                               row.names = 1, check.names = FALSE, quote = "\""))
  list(lines = lines, expr = expr)
}
ch_line_vals <- function(lines, pat) {
  ch <- grep("^!Sample_characteristics_ch1\t", lines, value = TRUE)
  out <- NULL
  for (l in ch) {
    v <- gsub('^"|"$', "", strsplit(l, "\t")[[1]][-1])
    if (any(grepl(pat, v, ignore.case = TRUE))) {
      out <- rep(NA_character_, length(v))
      hit <- grepl(pat, v, ignore.case = TRUE)
      out[hit] <- sub(".*:\\s*", "", v[hit])
    }
  }
  out
}

# ---------------- 队列 A: GSE46602 (基因级 RMA) ----------------
cat("== GSE46602 ==\n")
br <- readRDS(file.path(tab_dir, "bulk_recheck.rds"))
m466 <- br$mat; grp466 <- factor(br$group, levels = c("Benign", "Tumor"))
mcp466 <- MCPcounter.estimate(m466, featuresType = "HUGO_symbols")
cat("  populations:", nrow(mcp466), "| samples:", ncol(m466), "\n")

# ---------------- 队列 B: GSE70770 Cambridge (Illumina HumanHT-12 V4) ----------------
cat("== GSE70770 ==\n")
suppressPackageStartupMessages({ library(illuminaHumanv4.db); library(AnnotationDbi) })
sm <- read_sm(file.path(ext_d, "GSE70770-GPL10558_series_matrix.txt.gz"))
gsm <- get_field(sm$lines, "!Sample_geo_accession")
stype <- ch_line_vals(sm$lines, "sample type")
expr707 <- sm$expr
sym <- suppressMessages(AnnotationDbi::mapIds(illuminaHumanv4.db, keys = rownames(expr707),
                                              column = "SYMBOL", keytype = "PROBEID", multiVals = "first"))
keep <- !is.na(sym); expr707 <- expr707[keep, , drop = FALSE]; sym <- sym[keep]
o <- order(rowMeans(expr707, na.rm = TRUE), decreasing = TRUE)
expr707 <- expr707[o, , drop = FALSE]; sym <- sym[o]
expr707 <- expr707[!duplicated(sym), , drop = FALSE]; rownames(expr707) <- sym[!duplicated(sym)]
if (max(expr707, na.rm = TRUE) > 50) expr707 <- log2(expr707 + 1)
expr707 <- normalizeBetweenArrays(expr707, method = "quantile")
colnames(expr707) <- gsm
grp707 <- factor(ifelse(stype == "Tumour", "Tumor", "Benign"), levels = c("Benign", "Tumor"))
ok707 <- !is.na(grp707)
expr707 <- expr707[, ok707, drop = FALSE]; grp707 <- grp707[ok707]
cat("  tumour:", sum(grp707 == "Tumor"), "/ benign:", sum(grp707 == "Benign"),
    "| unlabelled (Stockholm) dropped:", sum(!ok707), "\n")
mcp707 <- MCPcounter.estimate(expr707, featuresType = "HUGO_symbols")

# ---------------- 队列 C: GSE104131 (RNA-seq FPKM, 配对) ----------------
cat("== GSE104131 ==\n")
fp <- fread(file.path(ext_d, "GSE104131_batch1.fpkms.txt.gz"))
setnames(fp, "GENESYMBOL", "gene")
fp <- fp[!duplicated(gene)]
sm104 <- read_sm(file.path(ext_d, "GSE104131_series_matrix.txt.gz"))
title <- sub("^Sample_", "", get_field(sm104$lines, "!Sample_title"))
tissue <- ch_line_vals(sm104$lines, "^tissue:")
meta104 <- data.table(sample = title, tissue = tissue)
meta104[, grp := fcase(tissue == "prostate tumor", "Tumor",
                       tissue == "normal prostate", "Benign",
                       default = NA_character_)]
meta104 <- meta104[!is.na(grp)]
common <- intersect(meta104$sample, setdiff(names(fp), c("GID", "Locus")))
meta104 <- meta104[sample %in% common]
fp104 <- as.matrix(fp[, ..common]); rownames(fp104) <- fp$gene[seq_len(nrow(fp))]
meta104 <- meta104[match(common, sample)]
grp104 <- factor(meta104$grp, levels = c("Benign", "Tumor"))
cat("  paired samples:", ncol(fp104), "| tumor:", sum(grp104 == "Tumor"),
    "/ benign:", sum(grp104 == "Benign"), "\n")
mcp104 <- MCPcounter.estimate(log2(fp104 + 1), featuresType = "HUGO_symbols")
# EPIC 输入: FPKM -> TPM (每样本列归一到 1e6, 不取 log)
tpm104 <- sweep(fp104, 2, colSums(fp104, na.rm = TRUE), "/") * 1e6
epic104 <- EPIC::EPIC(tpm104)
epic_frac <- as.data.table(epic104$cellFractions, keep.rownames = "sample")
cat("  EPIC cell types:", paste(setdiff(names(epic_frac), "sample"), collapse = ", "), "\n")

# ---------------- 汇总: MCP-counter 长表 + 统计 ----------------
pop_order <- rownames(mcp466)
longify <- function(m, cohort, grp) {
  d <- as.data.table(m, keep.rownames = "population")
  d <- melt(d, id.vars = "population", variable.name = "gsm", value.name = "score")
  d[, cohort := cohort][, group := as.character(grp)[match(gsm, colnames(m))]]
  d
}
mcp_all <- rbindlist(list(
  longify(mcp466, "GSE46602", grp466),
  longify(mcp707, "GSE70770", grp707),
  longify(mcp104, "GSE104131", grp104)))
mcp_all[, population := factor(population, levels = pop_order)]
fwrite(mcp_all, file.path(tab_dir, "T46k_mcp_scores.csv"))

stat_one <- function(d) {
  t <- d[group == "Tumor"]$score; b <- d[group == "Benign"]$score
  p <- suppressWarnings(wilcox.test(t, b)$p.value)
  data.table(median_benign = round(median(b, na.rm = TRUE), 3),
             median_tumor  = round(median(t, na.rm = TRUE), 3),
             log2FC_median = round(log2(median(t, na.rm = TRUE) / median(b, na.rm = TRUE)), 2),
             p = signif(p, 3))
}
mcp_stats <- mcp_all[, stat_one(.SD), by = .(cohort, population)]
mcp_stats[, p_adj := signif(p.adjust(p, method = "BH"), 3), by = cohort]
setorder(mcp_stats, cohort, population)
fwrite(mcp_stats, file.path(tab_dir, "T46k_mcp_stats.csv"))
cat("\n=== MCP-counter: tumor vs benign ===\n"); print(mcp_stats)

# EPIC 统计
epic_l <- melt(epic_frac, id.vars = "sample", variable.name = "population", value.name = "fraction")
epic_l <- merge(epic_l, meta104[, .(sample, grp)], by = "sample")
epic_stats <- epic_l[, {
  t <- fraction[grp == "Tumor"]; b <- fraction[grp == "Benign"]
  .(median_benign = round(median(b), 4), median_tumor = round(median(t), 4),
    p = signif(suppressWarnings(wilcox.test(t, b)$p.value), 3))
}, by = population]
epic_stats[, p_adj := signif(p.adjust(p, method = "BH"), 3)]
fwrite(epic_frac, file.path(tab_dir, "T46k_epic_scores.csv"))
fwrite(epic_stats, file.path(tab_dir, "T46k_epic_stats.csv"))
cat("\n=== EPIC (GSE104131) ===\n"); print(epic_stats)

# ---------------- Fig. S5 ----------------
theme_pub <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 8, hjust = 0),
        plot.subtitle = element_text(size = 7, hjust = 0),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6.5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8))
mcp_all[, cohort := factor(cohort, levels = c("GSE46602", "GSE70770", "GSE104131"))]
pA <- ggplot(mcp_all, aes(population, score + 0.01, fill = group)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = c(Benign = "#BDD7E7", Tumor = "#B2182B"), name = NULL) +
  scale_y_log10() +
  facet_wrap(~cohort, nrow = 1) +
  labs(title = "MCP-counter population scores",
       subtitle = "Tumour versus benign tissue; identical gene set per cohort",
       x = NULL, y = "Population score (log10)") +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

epic_l[, population := sub("x", "", as.character(population))]
pB <- ggplot(epic_l[population != "otherCells"],
             aes(population, fraction, fill = grp)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = c(Benign = "#BDD7E7", Tumor = "#B2182B"), name = NULL) +
  labs(title = "EPIC cell fractions (GSE104131)",
       subtitle = "RNA-seq cohort, 15 paired specimens; uncharacterised cells omitted",
       x = NULL, y = "Estimated fraction") +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5))

pw <- (pA / pB) + patchwork::plot_layout(heights = c(1.15, 1)) +
  patchwork::plot_annotation(tag_levels = "A")
tmp <- tempdir()
ggsave(file.path(tmp, "FigS5_deconvolution.pdf"), pw, width = 7.2, height = 5.6)
ragg::agg_jpeg(file.path(tmp, "FigS5_deconvolution.jpg"),
               width = 7.2 * 300, height = 5.6 * 300, res = 300)
print(pw); dev.off()
ok1 <- file.copy(file.path(tmp, "FigS5_deconvolution.pdf"), fig_dir, overwrite = TRUE)
ok2 <- file.copy(file.path(tmp, "FigS5_deconvolution.jpg"), fig_dir, overwrite = TRUE)
cat("figures copied:", ok1, ok2, "\n")
cat("DONE\n")
