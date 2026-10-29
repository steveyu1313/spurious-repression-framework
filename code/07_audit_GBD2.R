suppressPackageStartupMessages(library(data.table))
inc <- fread("./01_GBD_Analysis/data/GBD2023_prostate_incidence.csv")
asr <- inc[age == "Age-standardized" & metric == "Rate"]
w <- dcast(asr, location ~ year, value.var = "val")
mat <- as.matrix(w[, -1]); rownames(mat) <- w$location

cat("=== 位置间时间序列完全相同检测 (可能的重行/错配) ===\n")
loc <- rownames(mat); n <- length(loc); hit <- 0
for (i in 1:(n - 1)) for (j in (i + 1):n) {
  if (max(abs(mat[i, ] - mat[j, ])) < 1e-6) { cat("完全相同:", loc[i], "<->", loc[j], "\n"); hit <- hit + 1 }
}
cat("完全相同的对数:", hit, "\n")

cat("\n=== 各 location 2023 病例数与隐含人口 ===\n")
num <- inc[year == 2023 & age == "All ages" & metric == "Number", .(location, cases = val)]
rt  <- inc[year == 2023 & age == "All ages" & metric == "Rate", .(location, rate = val)]
x <- merge(num, rt, by = "location")
x[, pop_bn := cases / rate * 1e5 / 1e9]
print(as.data.frame(x[order(-pop_bn), .(location, cases = round(cases), rate = round(rate, 2),
                                        pop_bn = round(pop_bn, 3))]))

cat("\n=== SDI 五分位: 2023 病例数与 ASR ===\n")
sdi <- c("High SDI", "High-middle SDI", "Middle SDI", "Low-middle SDI", "Low SDI")
s <- merge(num[location %in% sdi], asr[year == 2023 & location %in% sdi, .(location, ASR = val)], by = "location")
s[, ord := match(location, sdi)]; s <- s[order(ord)]
s[, pct_cases := round(100 * cases / sum(cases), 1)]
print(as.data.frame(s[, .(location, ASR = round(ASR, 2), cases = round(cases), pct_cases)]))

cat("\n=== Low SDI 逐年 ASR vs 其他五分位 ===\n")
print(as.data.frame(dcast(asr[location %in% sdi], location ~ year, value.var = "val")[
  , lapply(.SD, function(z) round(z, 2)), by = location]))

cat("\n=== 正文死亡数核对 ===\n")
dea <- fread("./01_GBD_Analysis/data/GBD2023_prostate_deaths.csv")
print(as.data.frame(dea[location == "Global" & year %in% c(2013, 2023) & age == "All ages" & metric == "Number",
                        .(year, deaths = round(val))]))
print(as.data.frame(dea[location == "Global" & year %in% c(2013, 2023) & age == "Age-standardized" & metric == "Rate",
                        .(year, ASMR = round(val, 2))]))
cat("\n=== 发病率/死亡率比 (2023, 全球) ===\n")
i23 <- inc[location == "Global" & year == 2023 & age == "Age-standardized" & metric == "Rate", val]
m23 <- dea[location == "Global" & year == 2023 & age == "Age-standardized" & metric == "Rate", val]
cat("ASMR/ASIR =", round(m23 / i23, 3), "(前列腺癌通常约 0.2-0.35)\n")
