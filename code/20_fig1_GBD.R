# =============================================================================
# 20_fig1_GBD.R — Figure 1 重绘: 全球疾病负担 (面板合并 + 去冗余)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
  library(ragg); library(ggrepel)
})
set.seed(1)

fig_dir <- "./08_QC_Revision/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

FONT <- "Arial"
theme_pub <- function(base = 7) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks.length = unit(1.2, "mm"),
          axis.text = element_text(colour = "black", size = base),
          axis.title = element_text(colour = "black", size = base + 0.5),
          plot.title = element_text(size = base + 1, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = base - 0.5, colour = "grey25", hjust = 0),
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = base - 0.5),
          legend.title = element_text(size = base - 0.3),
          plot.margin = margin(2, 3, 2, 2),
          strip.background = element_blank())
}
PAL <- c(up = "#C0392B", down = "#2E6DA4", ns = "grey62")

rd <- function(f, m) {
  d <- fread(file.path("./01_GBD_Analysis/data", f))
  d[metric == m]
}
inc <- rd("GBD2023_prostate_incidence.csv", "Rate")
dea <- rd("GBD2023_prostate_deaths.csv", "Rate")
dal <- rd("GBD2023_prostate_DALYs.csv", "Rate")

# ---------------- A: 全球 ASIR / ASMR 趋势 ----------------
g <- rbind(inc[location == "Global" & age == "Age-standardized", .(year, val, lower, upper, metric = "ASIR")],
           dea[location == "Global" & age == "Age-standardized", .(year, val, lower, upper, metric = "ASMR")])
pA <- ggplot(g, aes(year, val, colour = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.0, stroke = 0) +
  annotate("text", x = 2012.2, y = 15.37, label = "ASIR 15.37", family = FONT, size = 2.2,
           hjust = 0, vjust = -0.6, colour = "#C0392B") +
  annotate("text", x = 2022.8, y = 15.33, label = "15.33", family = FONT, size = 2.2,
           hjust = 1, vjust = -0.7, colour = "#C0392B") +
  annotate("text", x = 2022.8, y = 5.29, label = "5.29", family = FONT, size = 2.2,
           hjust = 1, vjust = -0.7, colour = "#2E6DA4") +
  scale_colour_manual(values = c(ASIR = "#C0392B", ASMR = "#2E6DA4"), guide = "none") +
  scale_fill_manual(values = c(ASIR = "#C0392B", ASMR = "#2E6DA4"), guide = "none") +
  scale_x_continuous(breaks = seq(2012, 2022, 2), expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(limits = c(0, 19), breaks = seq(0, 15, 5), expand = c(0, 0)) +
  labs(x = NULL, y = "Age-standardised rate (per 100,000)",
       title = "Global temporal trends, 2012\u20132023",
       subtitle = "Shaded bands: 95% uncertainty intervals") +
  theme_pub() +
  annotate("text", x = 2012.2, y = 17.2, label = "ASIR", family = FONT, size = 2.6,
           hjust = 0, colour = "#C0392B", fontface = "bold") +
  annotate("text", x = 2012.2, y = 6.4, label = "ASMR", family = FONT, size = 2.6,
           hjust = 0, colour = "#2E6DA4", fontface = "bold")

# ---------------- B: 全球 ASDR 趋势 ----------------
d <- dal[location == "Global" & age == "Age-standardized"]
pB <- ggplot(d, aes(year, val)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#7B5EA7", alpha = 0.18) +
  geom_line(colour = "#5B3E8E", linewidth = 0.7) +
  geom_point(colour = "#5B3E8E", size = 1.0) +
  annotate("text", x = 2022.8, y = 97.43, label = "97.4", family = FONT, size = 2.2,
           hjust = 1, vjust = -0.7, colour = "#5B3E8E") +
  scale_x_continuous(breaks = seq(2012, 2022, 2), expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(limits = c(80, 115), breaks = seq(80, 110, 10)) +
  labs(x = NULL, y = "Age-standardised DALY rate\n(per 100,000)",
       title = "Global DALY rate", subtitle = " ") +
  theme_pub()

# ---------------- C: EAPC forest ----------------
t7 <- fread("./08_QC_Revision/tables/T7_numbers_table.csv")
grp_lab <- c(Country = "Country", `GBD region` = "GBD region", `SDI quintile` = "SDI quintile")
f <- copy(t7)
f[, short := location]
f[short == "Central Europe, Eastern Europe, and Central Asia", short := "Central/Eastern Europe & C. Asia"]
f[short == "Southeast Asia, East Asia, and Oceania", short := "SE/East Asia & Oceania"]
f[short == "North Africa and Middle East", short := "North Africa & Middle East"]
f[short == "Latin America and Caribbean", short := "Latin America & Caribbean"]
f[short == "United States of America", short := "United States"]
f[short == "High-middle SDI", short := "High-middle SDI"]
f[, dir := fcase(EAPC > 0 & p < 0.05, "up", EAPC < 0 & p < 0.05, "down", default = "ns")]
setorder(f, -EAPC)
f[, y := .N:1]
f[, ylab := sprintf("%s  (%.2f)", short, ASR2023)]

pC <- ggplot(f, aes(EAPC, y, colour = dir)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.45) +
  geom_point(size = 1.5, stroke = 0) +
  geom_text(aes(x = 5.35, label = sprintf("%.2f", EAPC)), family = FONT, size = 2.15,
            hjust = 1, colour = "black") +
  scale_colour_manual(values = PAL, guide = "none") +
  scale_y_continuous(breaks = f$y, labels = f$ylab, expand = expansion(add = 0.6)) +
  scale_x_continuous(limits = c(-2.0, 5.5), breaks = seq(-2, 5, 2)) +
  labs(x = "EAPC of age-standardised incidence (% per year, 95% CI)", y = NULL,
       title = "Estimated annual percentage change by location",
       subtitle = "Red: significant increase; blue: significant decrease; grey: not significant") +
  theme_pub(6.4) +
  theme(panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"),
        axis.text.y = element_text(size = 5.9))

# ---------------- D: 2023 ASIR 排序 ----------------
f2 <- copy(t7)[group != "Global"]
f2[, grp2 := factor(group, levels = c("Country", "GBD region", "SDI quintile"))]
f2[, short := f$short[match(location, f$location)]]
setorder(f2, ASR2023)
f2[, y := 1:.N]
f2[, dir := fcase(EAPC > 0 & p < 0.05, "up", EAPC < 0 & p < 0.05, "down", default = "ns")]

pD <- ggplot(f2, aes(ASR2023, y)) +
  geom_segment(aes(x = 0, xend = ASR2023, yend = y), linewidth = 0.3, colour = "grey75") +
  geom_point(aes(colour = dir, shape = grp2), size = 1.7, stroke = 0) +
  geom_text(aes(label = sprintf("%.1f", ASR2023)), family = FONT, size = 2.1,
            hjust = -0.35, vjust = 0.4) +
  scale_colour_manual(values = PAL, guide = "none") +
  scale_shape_manual(values = c(Country = 16, `GBD region` = 17, `SDI quintile` = 15), name = NULL) +
  scale_y_continuous(breaks = f2$y, labels = f2$short, expand = expansion(add = 0.6)) +
  scale_x_continuous(limits = c(0, 64), breaks = seq(0, 60, 20),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = "ASIR in 2023 (per 100,000)", y = NULL,
       title = "Incidence rate by location, 2023",
       subtitle = "\u25CF country   \u25B2 GBD region   \u25A0 SDI quintile") +
  theme_pub(6.4) +
  theme(panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"),
        axis.text.y = element_text(size = 5.9),
        legend.position = "none")

# ---------------- 组装 ----------------
fig1 <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1.55), widths = c(1, 1)) +
  plot_annotation(tag_levels = "A", tag_prefix = "",
                  theme = theme(plot.tag = element_text(family = FONT, face = "bold", size = 9)))

W <- 7.09; H <- 7.7
tmp <- file.path(tempdir(), "Fig1_global_burden.png")
ggsave(tmp, fig1, width = W, height = H, dpi = 600, device = ragg::agg_png, bg = "white")
file.copy(tmp, file.path(fig_dir, "Fig1_global_burden.png"), overwrite = TRUE)
tmp2 <- file.path(tempdir(), "Fig1_global_burden.pdf")
ggsave(tmp2, fig1, width = W, height = H, device = cairo_pdf)
file.copy(tmp2, file.path(fig_dir, "Fig1_global_burden.pdf"), overwrite = TRUE)
cat("Fig1 saved\n")
