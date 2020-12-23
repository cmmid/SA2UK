suppressPackageStartupMessages({
  require(data.table)
  require(lubridate)
})

.debug <- c("~/Dropbox/covidLMIC", "PAK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/interventions.rds",
  "%s/inputs/ecdc_data.rds",
  .debug[2],
  "%s/outputs/intervention_timing/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

outcomes <- readRDS(.args[2])[iso3 == tariso]

min.date <- as.Date("2020-02-01")

ref.time <- readRDS(.args[1])[iso3 == tariso]

if (dim(ref.time)[1]) {
  transition <- ref.time[, .(iso3, start = as.Date(date + intervention), end = as.Date(date + 6 + intervention), era = "transition")]
  post <- copy(transition)[, .(iso3, start = as.Date(end+1), end = as.Date(end + 31), era = "post")]
  censor <- outcomes[which.max(cumsum(cases) > 10)][, .(iso3, start = min.date, end = date + ref.time$censor, era = "censor")]
  pre <- censor[transition[, .(iso3, end = start)], on=.(iso3)][, .(iso3, start = as.Date(end+1), end = as.Date(i.end-1), era = "pre")]
  
  eras <- rbind(
    censor, pre, transition, post
  )
  
  if (ref.time[,!is.na(modification)]) {
    mod <- ref.time[, .(iso3, start = as.Date(modification), end = as.Date(modification) + 30, era = "modification")]
    eras <- rbind(eras, mod)
  }
  
} else {
  eras <- data.table(iso3=character(), start = Date(), end = Date(), era = character())
}

saveRDS(eras, tail(.args, 1))
