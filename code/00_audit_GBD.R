# =============================================================================
# 00_audit_GBD.R  —  GBD 2023 前列腺癌数据合理性核查
# 目的: 复算正文中所有 GBD 数字, 检查内部一致性、非单调性、量级合理性
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

root <- "."  # 相对路径, 规避 Windows 中文路径编码问题
out  <- file.path(root, "08_QC_Revision")
dir.create(file.path(out, "tables"), showWarnings = FALSE, recursive = TRUE)

rd <- function(f) fread(file.path(root, "01_GBD_Analysis/data", f), encoding = "UTF-8")

inc <- rd("GBD2023_prostate_incidence.csv")
dea <- rd("GBD2023_prostate_deaths.csv")
dal <- rd("GBD2023_prostate_DALYs.csv")

cat("================ 1. 数据结构 ================\n")
cat("incidence rows:", nrow(inc), "\n")
cat("metrics:", paste(unique(inc$metric), collapse = " | "), "\n")
cat("years:", min(inc$year), "-", max(inc$year), " (n =", length(unique(inc$year)), ")\n")
cat("sex:", paste(unique(inc$sex), collapse = " | "), "\n")
cat("n locations:", length(unique(inc$location)), "\n")
cat("locations:\n"); print(sort(unique(inc$location)))
cat("age groups:\n"); print(unique(inc$age))

cat("\n================ 2. 指标齐全性 ================\n")
print(table(inc$metric))
print(table(inc$age))

cat("\n================ 3. 2023 年全球关键值 ================\n")
key <- inc[location == "Global" & year == 2023 & age == "Age-standardized"]
print(key[, .(measure, metric, age, year, val, lower, upper)])

cat("\n--- 2013 全球(核对正文 2013:15.37) ---\n")
print(inc[location == "Global" & year %in% c(2013, 2023) & age == "Age-standardized" & metric == "Rate"])

cat("\n--- 全球死亡 2023 ---\n")
print(dea[location == "Global" & year == 2023 & age == "Age-standardized"])
cat("\n--- 全球 DALY 2023 ---\n")
print(dal[location == "Global" & year == 2023 & age == "Age-standardized"])

cat("\n================ 4. 全球粗率/年龄标化率换算校验 ================\n")
# 用 Number 与 Rate 反推分母人口, 检验 Rate 是否以全人群为分母
chk <- inc[location == "Global" & age == "All ages" & year == 2023]
print(chk[, .(metric, val)])
num <- chk[metric == "Number", val]
crd <- chk[metric == "Rate", val]
if (length(num) && length(crd)) {
  pop_implied <- num / crd * 1e5
  cat(sprintf("Number = %.0f, All-ages Rate = %.3f -> implied denominator population = %.3f billion\n",
              num, crd, pop_implied / 1e9))
}

cat("\n================ 5. SDI 五分位 2023 年龄标化发病率 ================\n")
sdi <- inc[age == "Age-standardized" & metric == "Rate" &
             location %in% c("High SDI", "High-middle SDI", "Middle SDI", "Low-middle SDI", "Low SDI")]
w <- dcast(sdi[year == 2023], location ~ ., value.var = "val")
setnames(w, ".", "ASR_2023")
sdi_order <- c("High SDI", "High-middle SDI", "Middle SDI", "Low-middle SDI", "Low SDI")
w[, order := match(location, sdi_order)]
setorder(w, order)
print(w)
cat("\n单调性检验 (应随 SDI 单调下降):\n")
v <- w$ASR_2023
cat("ASR 序列:", paste(round(v, 2), collapse = " > "), "\n")
cat("是否严格单调递减:", all(diff(v) < 0), "\n")
cat(sprintf("High SDI / Low-middle SDI = %.2f 倍\n", v[1] / v[4]))
cat(sprintf("High SDI / Low SDI = %.2f 倍\n", v[1] / v[5]))

cat("\n--- SDI 各年份 ASR 明细 ---\n")
print(dcast(sdi, location ~ year, value.var = "val")[, lapply(.SD, function(x) round(x, 2)),
      .SDcols = as.character(sort(unique(sdi$year)))])

cat("\n================ 6. EAPC 重算 (log-linear) ================\n")
asr <- inc[age == "Age-standardized" & metric == "Rate"]
eapc_re <- asr[, {
  if (.N >= 3 && all(val > 0)) {
    fit <- lm(log(val) ~ year)
    b <- coef(fit)[2]; se <- summary(fit)$coefficients[2, 2]
    .(EAPC = 100 * (exp(b) - 1),
      lo = 100 * (exp(b - 1.96 * se) - 1),
      hi = 100 * (exp(b + 1.96 * se) - 1),
      p = summary(fit)$coefficients[2, 4], n_yr = .N)
  } else .(EAPC = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, n_yr = .N)
}, by = location]
setorder(eapc_re, -EAPC)
print(eapc_re[!is.na(EAPC)][1:26])
fwrite(eapc_re, file.path(out, "tables/T1_EAPC_recomputed.csv"))

cat("\n--- 与既有 eapc_incidence.csv 对照 ---\n")
old <- fread(file.path(root, "01_GBD_Analysis/results/tables/eapc_incidence.csv"))
cmp <- merge(old[, .(location, EAPC_old = EAPC, p_old = pval)],
             eapc_re[, .(location, EAPC_new = EAPC, p_new = p)],
             by = "location", all = TRUE)
cmp[, diff := round(EAPC_new - EAPC_old, 4)]
print(cmp[order(-abs(diff))])
cat("\n最大绝对差:", max(abs(cmp$diff), na.rm = TRUE), "\n")
fwrite(cmp, file.path(out, "tables/T2_EAPC_recompute_vs_original.csv"))

cat("\n================ 7. 正文数字逐项核对 ================\n")
chk_num <- function(label, claimed, actual, tol = 0.02) {
  ok <- !is.na(actual) && abs(claimed - actual) <= tol * max(abs(claimed), 1e-9)
  cat(sprintf("[%s] 正文=%.4f  数据=%.4f  %s\n", label, claimed, actual, ifelse(ok, "一致", "★不一致★")))
}
g23 <- inc[location == "Global" & year == 2023 & age == "Age-standardized" & metric == "Rate", val]
g13 <- inc[location == "Global" & year == 2013 & age == "Age-standardized" & metric == "Rate", val]
chk_num("全球ASIR 2023", 15.33, g23)
chk_num("全球ASIR 2013", 15.37, g13)
m23 <- dea[location == "Global" & year == 2023 & age == "Age-standardized" & metric == "Rate", val]
chk_num("全球ASMR 2023", 5.29, m23)
d23 <- dal[location == "Global" & year == 2023 & age == "Age-standardized" & metric == "Rate", val]
chk_num("全球ASDR 2023", 97.43, d23)
chk_num("High SDI ASIR 2023", 22.2, w[location == "High SDI", ASR_2023])
chk_num("Low SDI ASIR 2023", 9.69, w[location == "Low SDI", ASR_2023])
chk_num("High-middle SDI ASIR", 9.13, w[location == "High-middle SDI", ASR_2023])
chk_num("Middle SDI ASIR", 7.66, w[location == "Middle SDI", ASR_2023])
chk_num("Low-middle SDI ASIR", 7.29, w[location == "Low-middle SDI", ASR_2023])

cat("\n--- 2023 全球发病/死亡数 (正文 150 万 / 40 万) ---\n")
cat("cases 2023 (all ages):", inc[location == "Global" & year == 2023 & age == "All ages" & metric == "Number", val], "\n")
cat("deaths 2023 (all ages):", dea[location == "Global" & year == 2023 & age == "All ages" & metric == "Number", val], "\n")

cat("\n================ 8. 国家层面 2023 ASIR 排序 ================\n")
cty <- inc[age == "Age-standardized" & metric == "Rate" & year == 2023]
cty <- cty[!grepl("SDI|Global|regions|Asia|Europe|Africa|Caribbean|Oceania|income|Latin", location)]
setorder(cty, -val)
print(head(cty[, .(location, ASIR = round(val, 2))], 12))
cat("\n最低:\n")
print(head(cty[, .(location, ASIR = round(val, 2))], 5)[order(ASIR)])

cat("\n================ 9. 95% UI 宽度合理性 ================\n")
uw <- inc[age == "Age-standardized" & metric == "Rate" & year == 2023]
uw[, width := (upper - lower) / val]
print(summary(uw$width))
cat("UI 宽度 > 100% 的比例:", round(mean(uw$width > 1) * 100, 1), "%\n")

cat("\n完成. 输出目录:", file.path(out, "tables"), "\n")
