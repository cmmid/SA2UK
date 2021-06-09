suppressPackageStartupMessages({
  require(data.table)
  require(countrycode)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper")
.args <- if (interactive()) sprintf(c(
  "%s/ins/fertility.csv",
  "%s/gen/fertility.rds"
), .debug[1]) else commandArgs(trailingOnly = TRUE)

#' covidm population not differentiated by sex, so want total value (SexID == 3)
#' exclude aggregate populations (LocID >= 900) and Channel Islands (830)
wpp.dt <- fread(.args[1])[!(LocID >= 900) & LocID != 830]
#' goal is to have best approximation of expected total population fertility per day, for today
fert.dt <- wpp.dt[!is.na(CBR) & MidPeriod <= 2020, .SD[.N], by=LocID][, .(
  per_capita_day = CBR/1000/365.25, # CBR is births per 1k capita per year; need per capita per day
  iso3 = countrycode(LocID, "iso3n", "iso3c")
)]

saveRDS(fert.dt, tail(.args, 1))