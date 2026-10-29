# ============================================================
# 40_hub_ppi.R
# 目标:
#   1) 从 GSE46602 复算 DEG 中筛出显著差异基因
#   2) 调 STRING v12 API 构建 PPI 网络 (medium confidence, score>=400)
#   3) igraph 计算拓扑指标 -> 识别 hub gene
#   4) 输出 节点表 / 边表 / hub 排名 / MCODE 模块
# 运行: cd 项目根 && Rscript 08_QC_Revision/scripts/40_hub_ppi.R
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(igraph)
})

set.seed(20261020)
fig_dir <- "./08_QC_Revision/figures"
tab_dir <- "./08_QC_Revision/tables"
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. DEG ----------
deg <- fread(file.path(tab_dir, "T6_GSE46602_DEG_rechecked.csv"))
deg <- deg[!is.na(Symbol) & Symbol != "" & !duplicated(Symbol)]
deg[, sig := abs(logFC) > 1 & adj.P.Val < 0.05]
sig <- deg[sig == TRUE]
sig[, direction := ifelse(logFC > 0, "Up", "Down")]
cat(sprintf("全部可注释基因 %d ; 显著 DEG %d (上调 %d / 下调 %d)\n",
            nrow(deg), nrow(sig), sum(sig$direction == "Up"), sum(sig$direction == "Down")))
fwrite(sig, file.path(tab_dir, "T12_DEG_significant.csv"))

# ---------- 2. STRING PPI ----------
# 注意: 必须一次性提交全部基因。分批提交时 STRING 只返回批内互作,
#       跨批次的边会被静默丢弃(实测 8 批 258 条 vs 一次性 1736 条)。
fetch_string <- function(genes, score) {
  body <- paste0("identifiers=", paste(genes, collapse = "%0d"),
                 "&species=9606&required_score=", score,
                 "&caller_identity=prostate_qc_revision")
  tmp <- tempfile(fileext = ".tsv")
  system2("curl", c("-s", "-m", "300", "-X", "POST", "--data", shQuote(body),
                    "-o", shQuote(tmp),
                    shQuote("https://string-db.org/api/tsv/network")),
          stdout = FALSE, stderr = FALSE)
  if (!file.exists(tmp) || file.info(tmp)$size == 0) return(NULL)
  d <- try(fread(tmp, sep = "\t", header = TRUE, quote = ""), silent = TRUE)
  if (inherits(d, "try-error") || !"preferredName_A" %in% names(d)) return(NULL)
  unique(d, by = c("preferredName_A", "preferredName_B"))
}

cache <- file.path(tab_dir, "T13_STRING_edges_raw.tsv")
cache_low <- file.path(tab_dir, "T13b_STRING_edges_lowconf.tsv")
if (!file.exists(cache)) {
  edges <- fetch_string(sig$Symbol, 400)
  if (is.null(edges)) stop("STRING 请求失败")
  fwrite(edges, cache)
} else edges <- fread(cache)
if (!file.exists(cache_low)) {
  e_low <- fetch_string(sig$Symbol, 150)
  if (!is.null(e_low)) fwrite(e_low, cache_low)
} else e_low <- fread(cache_low)

# 只保留本次 DEG 内部的互作, 并统一无向去重
edges <- edges[preferredName_A %in% sig$Symbol & preferredName_B %in% sig$Symbol]
edges[, `:=`(a = pmin(preferredName_A, preferredName_B),
             b = pmax(preferredName_A, preferredName_B))]
edges <- unique(edges, by = c("a", "b"))
cat(sprintf("PPI 边 %d 条, 涉及节点 %d 个\n", nrow(edges),
            length(unique(c(edges$a, edges$b)))))

# ---------- 3. igraph 拓扑 ----------
g <- graph_from_data_frame(edges[, .(a, b, score)], directed = FALSE,
                           vertices = sig[, .(Symbol, logFC, adj.P.Val, direction)])
g <- simplify(g, edge.attr.comb = list(score = "max", "ignore"))

comp <- components(g)
g <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
cat(sprintf("最大连通分量: %d 节点 / %d 边\n", vcount(g), ecount(g)))

V(g)$degree      <- degree(g)
V(g)$betweenness <- betweenness(g, normalized = TRUE)
V(g)$closeness   <- closeness(g, normalized = TRUE)
V(g)$eigen       <- eigen_centrality(g)$vector
V(g)$transitivity <- transitivity(g, type = "local", isolates = "zero")

# MCC (CytoHubba 近似): 邻居对的共同邻居数之和
# 邻接矩阵向量化, 避免逐对 intersect: MCC(v) = (sum(C[N,N]) - sum(diag(C)[N])) / 2
A <- as_adjacency_matrix(g, sparse = TRUE)
A <- as.matrix(A)
CN <- A %*% A
diagCN <- diag(CN)
V(g)$MCC <- sapply(seq_len(vcount(g)), function(i) {
  N <- which(A[i, ] > 0)
  if (length(N) < 2) return(0)
  (sum(CN[N, N]) - sum(diagCN[N])) / 2
})

md <- as.data.table(igraph::as_data_frame(g, what = "vertices"))
setnames(md, "name", "Symbol")
for (m in c("degree", "betweenness", "closeness", "eigen", "MCC")) {
  md[, paste0("rank_", m) := frank(-get(m), ties.method = "min")]
}
md[, rank_mean := rowMeans(.SD), .SDcols = paste0("rank_", c("degree", "betweenness", "closeness", "eigen", "MCC"))]
setorder(md, rank_mean)
md[, hub_rank := 1:.N]
md[, hub := hub_rank <= 20]
setcolorder(md, c("hub_rank", "Symbol", "hub", "degree", "betweenness", "closeness",
                  "eigen", "MCC", "transitivity", "logFC", "adj.P.Val", "direction"))
fwrite(md, file.path(tab_dir, "T14_PPI_node_topology.csv"))
fwrite(edges[, .(preferredName_A, preferredName_B, score)], file.path(tab_dir, "T15_PPI_edges.csv"))

cat("\n=== Top 20 hub genes ===\n")
print(md[1:20, .(hub_rank, Symbol, degree, betweenness = round(betweenness, 4),
                 MCC, logFC = round(logFC, 2), direction)])

# ---------- 4. MCODE 式模块 ----------
# 以 degree 加权取 seed, 扩展至邻居中权重>=0.2*seed 的子图, 要求 >=5 节点
vnames <- as.character(V(g)$name)
nbr <- function(v) {
  idx <- which(vnames == v)[1]
  if (is.na(idx)) return(character(0))
  as.character(neighbors(g, idx)$name)
}
cores <- list()
for (v in md[1:min(60, .N), Symbol]) {
  nb <- nbr(v)
  if (length(nb) < 4) next
  sub <- induced_subgraph(g, unique(c(v, nb)))
  sub <- induced_subgraph(sub, which(degree(sub) >= 0.2 * max(degree(sub))))
  if (vcount(sub) >= 5) {
    key <- paste(sort(as.character(V(sub)$name)), collapse = "|")
    cores[[key]] <- as.character(V(sub)$name)
  }
}
cores <- cores[!duplicated(names(cores))]
mods <- data.table(module = paste0("M", seq_along(cores)),
                   size = sapply(cores, length),
                   genes = sapply(cores, function(x) paste(x, collapse = ", ")),
                   seed = sapply(cores, function(x) x[1]))
setorder(mods, -size)
fwrite(mods, file.path(tab_dir, "T16_PPI_modules.csv"))
cat(sprintf("\n模块数 %d ; 最大模块 %d 个基因\n", nrow(mods), mods$size[1]))

# ---------- 5. hub 的一级邻居(下游效应模块) ----------
hub_genes <- md[hub == TRUE, Symbol]
nb_down <- unique(unlist(lapply(hub_genes, nbr)))
nb_down <- setdiff(nb_down, hub_genes)
fwrite(data.table(hub = hub_genes), file.path(tab_dir, "T17_hub_genes.csv"))
fwrite(data.table(downstream_neighbour = sort(nb_down)),
       file.path(tab_dir, "T18_hub_first_neighbours.csv"))
cat(sprintf("hub 一级邻居 %d 个\n", length(nb_down)))

saveRDS(list(graph = g, nodes = md, edges = edges, modules = mods,
             hub = hub_genes, neighbours = nb_down, deg = deg, sig = sig),
        file.path(tab_dir, "ppi_objects.rds"))
cat("\n完成 -> T12..T18\n")
