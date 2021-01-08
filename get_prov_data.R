suppressPackageStartupMessages({
  require(data.table)
})

.debug <- "~/Dropbox/SA2UK"
.args <- if (interactive()) sprintf(c(
  "%s/inputs/prov_data.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

extract <- function(url, meas) melt(
  fread(url)[, .SD, .SDcols = c(1, 3:11)],
  id.vars = "date", variable.name = "province"
)[, .(
  date = as.Date(date, "%d-%m-%Y"),
  value = c(value[1], diff(value))
  ),
  by = province
][, measure := meas ]
  
cases.dt <- extract(
  "https://raw.githubusercontent.com/dsfsi/covid19za/master/data/covid19za_provincial_cumulative_timeline_confirmed.csv",
  "cases"
)

deaths.dt <- extract(
  "https://raw.githubusercontent.com/dsfsi/covid19za/master/data/covid19za_provincial_cumulative_timeline_deaths.csv",
  "deaths"
)

final.dt <- dcast(rbind(cases.dt, deaths.dt), province + date ~ measure)
final.dt[is.na(deaths), deaths := 0 ]

saveRDS(final.dt, tail(.args, 1))
