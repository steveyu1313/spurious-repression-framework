#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================
# 103a_fetch_tcga_prad.py
# 经 cBioPortal API 抓取 TCGA-PRAD (pan-cancer atlas 2018):
#   - 4 个 panel 基因 (AMACR/PCA3/HOXC6/HPN) 的 RSEM RNA-seq 表达
#   - 全部 panel 基因表达 (用于抓 panel 邻近基因? 不需要, 只取 4 个)
#   - 患者级临床数据 (DFS_STATUS/DFS_MONTHS = 生化/临床复发)
#   - 样本级临床数据 (样本类型 primary tumor 等)
# 输出:
#   data/external/tcga/T45_tcga_prad_expression.csv   (样本 x 4 基因 RSEM)
#   data/external/tcga/T45b_tcga_prad_clinical.csv
# =============================================================
import json
import urllib.request
import os
import csv

BASE = "https://www.cbioportal.org/api"
STUDY = "prad_tcga_pan_can_atlas_2018"
PROFILE = STUDY + "_rna_seq_v2_mrna"
PANEL = ["AMACR", "PCA3", "HOXC6", "HPN"]
OUT_D = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     "data", "external", "tcga")
os.makedirs(OUT_D, exist_ok=True)


def api(path, body=None, params=None):
    url = BASE + path
    if params:
        parts = []
        for k, v in params.items():
            if isinstance(v, (list, tuple)):
                parts += [f"{k}={x}" for x in v]
            else:
                parts.append(f"{k}={v}")
        url += "?" + "&".join(parts)
    data = None
    headers = {"Accept": "application/json"}
    if body is not None:
        data = json.dumps(body).encode()
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers, method="POST" if data else "GET")
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print(f"  retry {attempt+1} for {url}: {e}")
    raise RuntimeError(f"API failed: {url}")


# ---------- 1. 分子数据 (GET /molecular-data, 每基因一次请求) ----------
ENTREZ = {"AMACR": 23600, "PCA3": 50652, "HOXC6": 3223, "HPN": 3249}
SYM_OF = {v: k for k, v in ENTREZ.items()}
print("fetching molecular data ...")
mol = []
for gid in ENTREZ.values():
    m = api(f"/molecular-profiles/{PROFILE}/molecular-data",
            params={"sampleListId": STUDY + "_all",
                    "entrezGeneId": gid, "projection": "SUMMARY"})
    print(f"  entrez {gid} ({SYM_OF[gid]}): rows {len(m)}")
    mol += m

# 透传为 样本 x 基因 矩阵
by_sample = {}
for row in mol:
    s = row["sampleId"]
    by_sample.setdefault(s, {})[SYM_OF[row["entrezGeneId"]]] = row["value"]
expr_path = os.path.join(OUT_D, "T45_tcga_prad_expression.csv")
with open(expr_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sampleId", "patientId"] + PANEL)
    for s in sorted(by_sample):
        vals = [by_sample[s].get(g, "") for g in PANEL]
        if all(v != "" for v in vals):          # 4 基因齐全才保留
            w.writerow([s, s[:12]] + vals)
print("  expression ->", expr_path, "| samples:", sum(1 for _ in open(expr_path)) - 1)

# ---------- 2. 患者级临床 (DFS/PFS/分级/分期在 PATIENT 级) ----------
print("fetching patient clinical ...")
OUT_ATTRS = ["DFS_STATUS", "DFS_MONTHS", "PFS_STATUS", "PFS_MONTHS",
             "DSS_STATUS", "DSS_MONTHS", "OS_STATUS", "OS_MONTHS",
             "GRADE", "AJCC_PATHOLOGIC_TUMOR_STAGE", "PATH_T_STAGE",
             "PATH_N_STAGE", "AGE_AT_DIAGNOSIS", "CANCER_TYPE"]
clin = api("/studies/" + STUDY + "/clinical-data",
           params={"clinicalDataType": "PATIENT", "projection": "DETAILED"})
pt = {}
for c in clin:
    if c["clinicalAttributeId"] in OUT_ATTRS:
        pt.setdefault(c["patientId"], {})[c["clinicalAttributeId"]] = c["value"]
samp_path = os.path.join(OUT_D, "T45b_tcga_prad_clinical.csv")
with open(samp_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["patientId"] + OUT_ATTRS)
    for s in sorted(pt):
        w.writerow([s] + [pt[s].get(k, "") for k in OUT_ATTRS])
print("  clinical ->", samp_path, "| patients:", len(pt))
print("DONE")
