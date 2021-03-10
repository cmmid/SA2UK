suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper", "GHA")
.args <- if (interactive()) sprintf(c(
  .debug[2],
  "%s/outputs/intervention_timing/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

transition_era <- function(
  dates, period, i3=tariso
) {
  cast.dates <- as.Date(dates)
  data.table(
    iso3=i3,
    start = cast.dates[1:3],
    end = c(cast.dates[2:3]-1, cast.dates[4]),
    era = c("pre", "transition", "post"),
    period = period
  )
}

censor_era <- function(from, first) {
  dt <- copy(first)[1]
  dt[, era := "censor" ]
  dt[, period := period - 1L ]
  dt[, end := start - 1 ]
  dt[, start := as.Date(from) ]
}

relaxation_era <- function(first, second) {
  ret <- copy(second)[1]
  ret[, period := period - 1L ]
  ret[, end := start - 1 ]
  ret[, start := first[.N, end+1] ]
  ret[, era := "relaxation" ]
  ret
}

if (.debug[2] == "ZAF") {
  
  initial_intervention <- transition_era(
    sprintf("2020-%02i-%02i", c(3,3,4,4), c(6,22,5,30)),
    1
  )

  censor <- censor_era(sprintf("2020-%02i-%02i",2, 1), initial_intervention)
  
  variant_intervention <- transition_era(
    c(
      sprintf("2020-%02i-%02i", c(11,12,12), c(22,17,29)),
      sprintf("2021-%02i-%02i", 1, 7)
    ),
    3
  )
  
  relax <- relaxation_era(initial_intervention, variant_intervention)
  
} else if (.debug[2] == "GHA") {
  
  initial_intervention <- transition_era(
    sprintf("2020-%02i-%02i", c(4,5,5,6), c(1,1,15,15)),
    1
  )
  
  censor <- censor_era(sprintf("2020-%02i-%02i",2, 1), initial_intervention)
  
  variant_intervention <- transition_era(
    sprintf("2021-%02i-%02i", c(1,1,2,2), c(1,21,1,15)),
    3
  )
  
  relax <- relaxation_era(initial_intervention, variant_intervention)
}

eras <- rbind(
  censor,
  initial_intervention,
  relax,
  variant_intervention
)

saveRDS(eras, tail(.args, 1))
