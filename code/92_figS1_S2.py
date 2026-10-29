# -*- coding: utf-8 -*-
"""92_figS1_S2.py — generate Supplementary FigS1 (scRNA QC) and FigS2 (spatial QC).
Run with miniconda python (pandas/matplotlib).
Outputs: supplementary/FigS1_scRNA_QC.png, supplementary/FigS2_spatial_QC.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = r"D:/yuproject/前列腺癌/08_QC_Revision"
TAB = BASE + r"/tables"
SUP = BASE + r"/supplementary"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "savefig.dpi": 300,
})

# ---------------- FigS1 : scRNA QC ----------------
x1 = pd.read_csv(TAB + "/X1_scrna_umap.csv")
order = (x1.groupby("gse_sample").size().sort_values(ascending=False).index.tolist())

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
ax = axes[0, 0]
data = [x1.loc[x1.gse_sample == s, "nFeature_RNA"].values for s in order]
vp = ax.violinplot(data, showmedians=True, widths=0.8)
for b in vp["bodies"]:
    b.set_facecolor("#4C78A8"); b.set_alpha(0.6)
ax.set_xticks(range(1, len(order) + 1)); ax.set_xticklabels(order, rotation=45, ha="right")
ax.set_ylabel("Detected genes per cell")
ax.set_title("A", loc="left", fontweight="bold")

ax = axes[0, 1]
data = [x1.loc[x1.gse_sample == s, "nCount_RNA"].values for s in order]
vp = ax.violinplot(data, showmedians=True, widths=0.8)
for b in vp["bodies"]:
    b.set_facecolor("#4C78A8"); b.set_alpha(0.6)
ax.set_xticks(range(1, len(order) + 1)); ax.set_xticklabels(order, rotation=45, ha="right")
ax.set_ylabel("UMI counts per cell")
ax.set_title("B", loc="left", fontweight="bold")

ax = axes[1, 0]
cnt = x1["gse_sample"].value_counts().reindex(order)
ax.bar(range(len(order)), cnt.values, color="#4C78A8")
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=45, ha="right")
ax.set_ylabel("Cells retained after filtering")
ax.set_title("C", loc="left", fontweight="bold")
for i, v in enumerate(cnt.values):
    ax.text(i, v + 60, str(v), ha="center", fontsize=7)

ax = axes[1, 1]
sub = x1.sample(min(8000, len(x1)), random_state=0)
types = sub["celltype_marker"].value_counts().index.tolist()
cmap = plt.get_cmap("tab20")
for k, ct in enumerate(types):
    d = sub[sub.celltype_marker == ct]
    ax.scatter(d["nCount_RNA"], d["nFeature_RNA"], s=3, alpha=0.35,
               color=cmap(k % 20), label=ct, rasterized=True)
ax.set_xlabel("UMI counts per cell"); ax.set_ylabel("Detected genes per cell")
ax.set_xlim(0, np.percentile(x1["nCount_RNA"], 99.5))
ax.set_ylim(0, np.percentile(x1["nFeature_RNA"], 99.5))
ax.legend(fontsize=5.5, ncol=2, loc="upper left", frameon=False, markerscale=3)
ax.set_title("D", loc="left", fontweight="bold")

fig.suptitle("Single-cell quality metrics after filtering (GSE141445; 36,424 cells, 13 tumours)",
             fontsize=11, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.98))
fig.savefig(SUP + "/FigS1_scRNA_QC.png")
plt.close(fig)
print("FigS1 done")

# ---------------- FigS2 : spatial QC before/after ----------------
x7 = pd.read_csv(TAB + "/X7_spatial_meta.csv")     # all captured spots
x8 = pd.read_csv(TAB + "/X8_spatial_reanalysed.csv")  # in-tissue spots used
secs = ["CDZ", "CFX", "HB", "LGL"]

fig, axes = plt.subplots(2, 2, figsize=(10, 8))

def box_before_after(ax, col_before, col_after, ylab):
    pos, data, colors = [], [], []
    for i, s in enumerate(secs):
        b = x7.loc[x7["sample"] == s, col_before].values
        a = x8.loc[x8["section"] == s, col_after].values
        pos += [i - 0.18, i + 0.18]
        data += [b, a]
        colors += ["#BBBBBB", "#E45756"]
    bp = ax.boxplot(data, positions=pos, widths=0.3, showfliers=False,
                    patch_artist=True, medianprops=dict(color="black", lw=1))
    for p, c in zip(bp["boxes"], colors):
        p.set_facecolor(c); p.set_alpha(0.75); p.set_edgecolor("none")
    ax.set_xticks(range(len(secs))); ax.set_xticklabels(secs)
    ax.set_ylabel(ylab)
    ax.plot([], [], "s", color="#BBBBBB", label="All captured spots (4,992/section)")
    ax.plot([], [], "s", color="#E45756", label="In-tissue spots used")
    ax.legend(fontsize=7, frameon=False)

ax = axes[0, 0]
box_before_after(ax, "nFeature", "nFeature", "Detected genes per spot")
ax.set_title("A", loc="left", fontweight="bold")
ax = axes[0, 1]
box_before_after(ax, "nCount", "nCount", "UMI counts per spot")
ax.set_title("B", loc="left", fontweight="bold")
ax = axes[1, 0]
box_before_after(ax, "mt", "mt", "Mitochondrial fraction (%)")
ax.set_title("C", loc="left", fontweight="bold")

ax = axes[1, 1]
before = [4992] * len(secs)
after = [int((x8["section"] == s).sum()) for s in secs]
xp = np.arange(len(secs))
ax.bar(xp - 0.18, before, width=0.36, color="#BBBBBB", label="All captured spots")
ax.bar(xp + 0.18, after, width=0.36, color="#E45756", label="In-tissue spots used")
for i, v in enumerate(after):
    ax.text(i + 0.18, v + 60, str(v), ha="center", fontsize=8)
ax.set_xticks(xp); ax.set_xticklabels(secs)
ax.set_ylabel("Spots per section")
ax.legend(fontsize=7, frameon=False)
ax.set_title("D", loc="left", fontweight="bold")

fig.suptitle("Spatial spot-level quality metrics before and after restriction to in-tissue spots (GSE245388)",
             fontsize=11, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.98))
fig.savefig(SUP + "/FigS2_spatial_QC.png")
plt.close(fig)
print("FigS2 done")

# also derive Table S1 (scRNA QC by tumour) while we are here
t_s1 = x1.groupby("gse_sample").agg(
    n_cells=("barcode", "size"),
    median_genes=("nFeature_RNA", "median"),
    median_UMI=("nCount_RNA", "median"),
    min_genes=("nFeature_RNA", "min"),
    max_genes=("nFeature_RNA", "max"),
).reset_index().rename(columns={"gse_sample": "sample"})
t_s1.to_csv(SUP + "/Table_S1_scRNA_QC_by_tumour.csv", index=False)
print("Table S1 derived:", t_s1.shape)
