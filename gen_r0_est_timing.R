suppressPackageStartupMessages({
  require(data.table)
  require(lubridate)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  .debug[2],
  "%s/outputs/intervention_timing/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

transition_era <- function(
  dates, period, i3=tariso
) data.table(
  iso3=i3,
  start = as.Date(dates[1:3]),
  end = as.Date(c(dates[2:3]-1, dates[4])),
  era = c("pre", "transition", "post"),
  period = period
)

initial_intervention <- transition_era(
  sprintf("2020-%02i-%02i", c(3,3,4,4), c(6,22,5,30)),
  1
)

censor_era <- function(from, first) {
  dt <- copy(first)[1]
  dt[, era := "censor" ]
  dt[, period := period - 1L ]
  dt[, end := start - 1 ]
  dt[, start := as.Date(from) ]
}

variant_intervention <- transition_era(
  c(
    sprintf("2020-%02i-%02i", c(11,12,12), c(22,16,29)),
    sprintf("2021-%02i-%02i", 1, 7)
  ),
  3
)

#' hand specifying eras from other analyses
eras <- data.table(
  iso3 = tariso,
  start = as.Date(sprintf("2020-%02i-%02i",c(2,3,3,4,5,11), c(1,6,22,5,1,22))),
  end = as.Date(sprintf("2020-%02i-%02i",c(3,3,4,4,10,12), c(5,21,4,30,15,16))),
  era = c("censor", "pre", "transition", "post", "relaxation", "variant")
)



saveRDS(eras, tail(.args, 1))
