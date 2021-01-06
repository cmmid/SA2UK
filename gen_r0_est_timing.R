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

#' hand specifying eras from other analyses
eras <- data.table(
  iso3 = tariso,
  start = as.Date(sprintf("2020-%02i-%02i",c(2,3,3,4,6,11), c(1,6,22,5,1,22))),
  end = as.Date(sprintf("2020-%02i-%02i",c(3,3,4,5,10,12), c(5,21,4,5,15,9))),
  era = c("censor", "pre", "transition", "post", "relaxation", "variant")
)

saveRDS(eras, tail(.args, 1))
