suppressPackageStartupMessages(library(data.table))
rd <- function(f, m) { d <- fread(file.path("./01_GBD_Analysis/data", f)); d[metric == m] }
inc <- rd("GBD2023_prostate_incidence.csv", "Rate")
dea <- rd("GBD2023_prostate_deaths.csv", "Rate")
dal <- rd("GBD2023_prostate_DALYs.csv", "Rate")
get <- function(d, yr) d[location == "Global" & year == yr & age == "Age-standardized",
                         .(val = round(val, 2), lo = round(lower, 2), hi = round(upper, 2))]
cat("ASIR 2012:", unlist(get(inc, 2012)), "\n")
cat("ASIR 2023:", unlist(get(inc, 2023)), "\n")
cat("ASMR 2012:", unlist(get(dea, 2012)), "\n")
cat("ASMR 2023:", unlist(get(dea, 2023)), "\n")
cat("ASDR 2012:", unlist(get(dal, 2012)), "\n")
cat("ASDR 2023:", unlist(get(dal, 2023)), "\n")
cat("\n2012/2023 相对变化: ASIR", round(100*(15.33/15.54-1),1), "%; ASMR",
    round(100*(5.29/5.6-1),1), "%\n")
# 精确值
cat("\nASIR 2012 精确:", inc[location=="Global" & year==2012 & age=="Age-standardized", val], "\n")
cat("ASMR 2012 精确:", dea[location=="Global" & year==2012 & age=="Age-standardized", val], "\n")
cat("ASDR 2012 精确:", dal[location=="Global" & year==2012 & age=="Age-standardized", val], "\n")
