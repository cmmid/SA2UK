suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper", "PAK")
.args <- if (interactive()) sprintf(c(
  "interventions.csv",
  .debug[2],
  "%s/outputs/intervention_timing/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

dt <- fread(.args[1])[iso3 == tariso][order(start)]
dt[era == "pre", period := 1+2*((1:.N)-1)]
dt[era == "post", period := 1+2*((1:.N)-1)]

#' assert:
#'  - multiple pre/post eras correspond to subsequent variants
#'  - pre first period is censoring era
#'  - between multiple pre/post intervals are transitions

#' @param dt data.table, 2 rows
#' @return data.table, the original dt with a transition era inserted
insert_transition_era <- function(dt) {
  ins <- copy(dt[1,])
  ins$start <- ins$end + 1
  ins$end <- dt[2, start-1]
  ins$era <- "transition"
  rbind(dt, ins)
}

prepend_censor_era <- function(dt, from) {
  cn <- copy(dt[order(start)][1,])
  cn$period <- cn$period - 1L
  cn$end <- cn$start - 1L
  cn$start <- as.IDate(from)
  cn$era <- "censor"
  setkeyv(rbind(cn, dt), key(dt))
}

insert_relaxation_eras <- function(dt, pre_variant_offset = 60) {
  if (dt[, max(period)] > 1) {
    ref <- dt[order(start)][which.max(period > 1) + -1:0 ]
    ins <- copy(ref[1,])
    ins$start <- ins$end + 1
    ins$end <- ref[2, start - pre_variant_offset]
    ins$period <- 2
    ins$era <- "relaxation"
    setkeyv(rbind(dt, ins), key(dt))
  } else dt 
}

eras <- prepend_censor_era(
  insert_relaxation_eras(
    setkey(dt[, insert_transition_era(.SD), by=period ], iso3, period, start, era)
  ),
  from = "2020-02-01"
)

saveRDS(eras, tail(.args, 1))
