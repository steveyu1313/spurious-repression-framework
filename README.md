# A multi-scale de-biasing framework for bulk-tissue biomarker and target discovery, with a zero-parameter transportable classifier (prostate cancer case study)

Code, result tables and figures for the manuscript:

> A multi-scale de-biasing framework for bulk-tissue biomarker and target discovery, with a zero-parameter transportable classifier: prostate cancer as a case study.

All analyses are secondary analyses of public, de-identified datasets. No wet-laboratory experiment was performed. The outputs are hypothesis-generating: a diagnostic panel candidate and an ordered target agenda, not validated targets or assays.

## What this repository contains

1. **Multi-scale de-biasing operator** — resolves *spurious target repression* in bulk transcriptomics: down-regulated "hub" genes that merely reflect stromal/immune compartment loss rather than tumour-cell repression. Direction, cell-of-origin (single-cell), spatial and druggability filters are combined into a corrected target shortlist (MYC–EZH2–FASN in prostate cancer).
2. **Zero-parameter transportable classifier** — a fixed-weight four-gene epithelial panel scored as the mean of within-cohort z-scores, benchmarked under an identical protocol against logistic regression, LASSO, random forest and radial-SVM across Affymetrix → Illumina → RNA-seq transfers, with paired DeLong tests, calibration and decision-curve analysis at every step.
3. **Established-tool deconvolution check** — MCP-counter and EPIC applied as published to bound the composition argument (Fig S5, Tables S21–S22).

## Repository layout

```
GitHub_release/
├── code/                 # full analysis pipeline, run order = file number
├── results/tables/       # all intermediate/final tables (T*/X* CSV/TSV)
├── figures/              # final manuscript and supplementary figures
├── LICENSE               # MIT
└── README.md
```

`results/tables/` contains every table produced by the pipeline, so all figures can be regenerated from these CSV/TSV files alone without re-downloading or re-running the upstream cohorts.

## Data sources (all public)

| Source | Accession / URL | Used by |
|---|---|---|
| GSE46602 (discovery, Affymetrix, laser-capture) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46602 | 05, 16–18, 100 |
| GSE69223 | GEO | bulk cohorts |
| GSE55945 | GEO | bulk cohorts |
| GSE70770 (external, Illumina HumanHT-12) | GEO | 82, 100, 106 |
| GSE104131 (RNA-seq, AA/EA paired) | GEO | 103b, 106 |
| GSE141445 (single-cell) | GEO | 01–03, 12, 52_panel |
| GSE245388 (spatial) | GEO | 04, 06, 13, 15 |
| TCGA-PRAD pan-cancer atlas 2018 | https://www.cbioportal.org | 103a (fetch), 103b |
| STRING PPI | https://string-db.org | 40 |
| ChEA3 / KEA3 upstream analysis | https://maayanlab.cloud | 41 |
| DGIdb druggability | https://www.dgidb.org | 61 |
| GDSC1/GDSC2/CTRP2 via oncoPredict training files | https://osf.io/c6tfx/ | 83 |
| Human Protein Atlas v24 | https://www.proteinatlas.org (CC BY-SA 3.0) | 84, 85, 86 |

GEO expression matrices are downloaded manually from the accessions above into the project data folder; `103a_fetch_tcga_prad.py` (cBioPortal API) and `84/85_fetch_hpa*.py` (HPA) fetch their data programmatically.

## Run order

| Stage | Scripts | Purpose |
|---|---|---|
| 0. Data audit | `00`–`07` | Structural and numeric consistency checks of GBD, scRNA, spatial and bulk inputs |
| 1. Bulk reanalysis | `05`, `16`, `17`, `18` | DEG recheck, outlier and sensitivity analyses, baseline contrast |
| 2. sc/spatial export & reanalysis | `12`, `13`, `15` | Cell-type expression matrices, spatial reanalysis |
| 3. Main figures 1–4 | `20`, `21`, `22`, `23`, `23b` | GBD epidemiology, bulk, single-cell, spatial figures |
| 4. Network / pathway / druggability | `40`–`44`, `61` | PPI hubs, GO/KEGG/GSEA, cross-layer hub table, DGIdb druggability |
| 5. Diagnostic panel & validation | `50`, `52`, `52_panel_scrna`, `51`, `80`, `82`, `83` | Panel construction, incremental value, external validation, drug-evidence correlation |
| 6. HPA protein evidence | `84`, `85`, `86` | Fetch and compose IHC figure |
| 7. Cross-platform benchmark | `100` | ZeroParam vs LR/LASSO/RF/SVM, GSE46602 → GSE70770 |
| 8. External RNA-seq validation | `103a`, `103b` | TCGA-PRAD fetch; GSE104131 + TCGA panel validation, PFS |
| 9. Deconvolution | `106` | MCP-counter + EPIC on three bulk cohorts |
| 10. Supplementary figures | `92` | Figures S1–S2 |

Figure scripts (`20`–`23b`, `43`, `44`, `51`, `52`, `86`, `92`) can be run directly from `results/tables/`.

## Software environment

- **R ≥ 4.2**: `data.table, limma, Seurat, harmony, clusterProfiler, fgsea, msigdbr, igraph, ggplot2, ggrepel, patchwork, scales, pROC, glmnet, randomForest, e1071, survival, stringr, jsonlite, png, ragg, AnnotationDbi, org.Hs.eg.db, hgu133plus2.db, illuminaHumanv4.db`
- **Python ≥ 3.10**: `numpy, pandas, matplotlib, requests/urllib`
- `106_deconvolution.R` sources MCP-counter directly from the authors' repository (https://github.com/ebecht/MCPcounter) at runtime; EPIC is installed from https://github.com/dtsoucas/EPIC (also on CRAN as `EPIC`).

## Citation

If you use this code, please cite the manuscript (citation to be updated upon publication). The repository will be archived on Zenodo with a DOI upon acceptance.

## License

MIT License — see [LICENSE](LICENSE).
