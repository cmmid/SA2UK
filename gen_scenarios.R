suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/SA2UK/outputs","ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/fits/%s.rds",
  "%s/intervention_timing/%s.rds", #' TODO remove
  "%s/introductions/%s.rds",
  .debug[2],
  "%s/scenarios/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

#' parameter fitting results
fit.dt <- readRDS(.args[1])

day0 <- readRDS(.args[3])[iso3 == tariso, min(date)]
#' target output
outfile <- tail(.args, 1)

translated.dt <- fit.dt[!is.na(work),.(
  q,
  self_iso = sympred,
  school = fifelse(school == "large", largered, smallred),
  home = 0,
  work = fifelse(work == "large", largered, smallred),
  other = fifelse(other == "large", largered, smallred),
  start_day = as.integer(date - day0),
  end_day = Inf,
  travel = 0,
  population = -1, coverage = 1,
  scen_type = "mitigated",
  trigger_type = NA_character_, trigger_value = NA_real_
)]

all.dt <- rbind(
  data.table(scen_type = "unmitigated", scen_id = 1L),
  translated.dt[, scen_id := 2L ],
  fill = TRUE
)

saveRDS(all.dt, outfile)
