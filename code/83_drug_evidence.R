# 83_drug_evidence.R  -- GDSC/CTRP expression-vs-IC50 evidence for MYC-EZH2-FASN
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ext_d <- "08_QC_Revision/data/external"
fig_d <- "08_QC_Revision/figures"
tab_d <- "08_QC_Revision/tables"

# GDSC2: EZH2 inhibitors -------------------------------------------------
e2 <- readRDS(file.path(ext_d, "GDSC2_Expr (RMA Normalized and Log Transformed).rds"))
r2 <- readRDS(file.path(ext_d, "GDSC2_Res.rds")); r2 <- exp(r2)  # IC50 uM
ez2_drugs <- grep("GSK343|EPZ", colnames(r2), value = TRUE, ignore.case = TRUE)
cat("GDSC2 EZH2-inhibitor columns:", paste(ez2_drugs, collapse = ", "), "\n")

# CTRP2: FASN inhibitor (cerulenin) --------------------------------------
ec <- readRDS(file.path(ext_d, "CTRP2_Expr (TPM, not log transformed).rds"))
rc <- readRDS(file.path(ext_d, "CTRP2_Res.rds")); rc <- exp(rc)
cer <- grep("cerulenin", colnames(rc), value = TRUE, ignore.case = TRUE)
cat("CTRP2 cerulenin columns:", paste(cer, collapse = ", "), "\n")

# GDSC1: docetaxel (clinical context) -------------------------------------
e1 <- readRDS(file.path(ext_d, "GDSC1_Expr (RMA Normalized and Log Transformed).rds"))
r1 <- readRDS(file.path(ext_d, "GDSC1_Res.rds")); r1 <- exp(r1)
dtx <- grep("^Docetaxel", colnames(r1), value = TRUE)

get_expr <- function(E, gene) {   # oncoPredict Expr: genes (rows) x cell lines (cols)
  i <- grep(paste0("^", gene, "$"), rownames(E))
  if (!length(i)) return(NULL)
  v <- E[i[1], ]
  if (max(v, na.rm = TRUE) > 200) v <- log2(v + 1)   # TPM -> log2
  v
}
common <- function(a, b) intersect(names(a), names(b))

cor_pair <- function(E, gene, R, drug, dataset) {
  g <- get_expr(E, gene); d <- R[, drug]
  names(d) <- rownames(R); if (is.null(g)) return(NULL)
  k <- common(g, d); g <- g[k]; d <- d[k]
  ok <- !is.na(g) & !is.na(d)
  ct <- cor.test(g[ok], d[ok], method = "spearman")
  data.table(dataset = dataset, gene = gene, drug = drug, n = sum(ok),
             rho = unname(ct$estimate), p = ct$p.value)
}

rows <- list()
for (dg in ez2_drugs) rows[[length(rows)+1]] <- cor_pair(e2, "EZH2", r2, dg, "GDSC2")
for (dg in cer)      rows[[length(rows)+1]] <- cor_pair(ec, "FASN", rc, dg, "CTRP2")
for (dg in dtx) for (g in c("MYC","EZH2","FASN"))
  rows[[length(rows)+1]] <- cor_pair(e1, g, r1, dg, "GDSC1")
res <- rbindlist(rows)
res[, padj := p.adjust(p, "BH")]
res[, q_label := ifelse(p < 1e-4, "P < 1e-4", paste0("P = ", signif(p, 2)))]
fwrite(res, file.path(tab_d, "T45_drug_expression_corr.csv"))
cat("\n== T45 ==\n"); print(res)

# scatter figure (FigS4) --------------------------------------------------
mk_dt <- function(E, gene, R, drug, dataset) {
  g <- get_expr(E, gene); d <- R[, drug]; names(d) <- rownames(R)
  k <- common(g, d); ok <- !is.na(g[k]) & !is.na(d[k])
  data.table(dataset = dataset, drug = drug, gene = gene, expr = g[k][ok], ic50 = d[k][ok])
}
pl <- list()
for (dg in ez2_drugs) pl[[length(pl)+1]] <- mk_dt(e2, "EZH2", r2, dg, "GDSC2")
for (dg in cer)      pl[[length(pl)+1]] <- mk_dt(ec, "FASN", rc, dg, "CTRP2")
pl <- rbindlist(pl)
pl[, facet := paste0(dataset, " | ", gene, " vs ", sub("_[0-9]+$", "", drug))]
gg <- ggplot(pl, aes(expr, ic50)) +
  geom_point(size = 0.7, alpha = 0.45, colour = "#3D8B8B") +
  geom_smooth(method = "lm", colour = "#C94F4F", linewidth = 0.5, se = FALSE) +
  facet_wrap(~facet, scales = "free") +
  labs(x = "Cell-line expression (log2)", y = expression(IC50~" ("*mu*"M)"),
       title = "Expression of core genes vs inhibitor sensitivity (Spearman)",
       subtitle = paste0("EZH2 vs GSK343: rho = ", signif(res[drug=="GSK343_1627"]$rho, 2),
                         ", ", res[drug=="GSK343_1627"]$q_label,
                         "; FASN vs cerulenin: rho = ", signif(res[drug=="cerulenin"]$rho, 2),
                         ", ", res[drug=="cerulenin"]$q_label)) +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face="bold", size=10))
ggsave(file.path(fig_d, "FigS4_drug_evidence.png"), gg, width = 7.2, height = 3.2, dpi = 600)

# prostate-line subset (context) ------------------------------------------
pr <- c("22RV1","DU-145","LNCAP","VCAP","PC3","PC-3","DU145","LNCaP","VDJ")
res_pr <- copy(res); res_pr[, n_prostate := NA_integer_]; res_pr[, rho_prostate := NA_real_]
for (i in seq_len(nrow(res))) {
  R <- switch(res$dataset[i], GDSC2 = r2, CTRP2 = rc, GDSC1 = r1)
  E <- switch(res$dataset[i], GDSC2 = e2, CTRP2 = ec, GDSC1 = e1)
  g <- get_expr(E, res$gene[i]); d <- R[, res$drug[i]]; names(d) <- rownames(R)
  k <- common(g, d); ok <- !is.na(g[k]) & !is.na(d[k]) & (k %in% pr)
  if (sum(ok) >= 6) {
    res_pr[i, `:=`(n_prostate = sum(ok), rho_prostate = unname(cor(g[k][ok], d[k][ok], method = "spearman")))]
  }
}
fwrite(res_pr, file.path(tab_d, "T45_drug_expression_corr.csv"))
cat("\n== T45 with prostate subset ==\n"); print(res_pr)
cat("FigS4 + T45 written\n")
