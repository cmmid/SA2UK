suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("analysis", "GHA")
.args <- if (interactive()) sprintf(c(
  "%s/ins/interventions.rds",
  "%s/gen/pops/%s.rds",
  .debug[2],
  "%s/scenarios/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]
pop <- readRDS(.args[1])

#' scenarios for population targets

scenarios <- list(
  efficacy = list(),
  rollout = list(),
  rolloutspeed = list()
)


names(pop$pop[[1]]$size) <- pop$pop[[1]]$group_names
scen1 <- round(sum(pop$pop[[1]]$size[13:16])*.9)
scen2 <- scen1 + sum(pop$pop[[1]]$size[4:12])*.2
scen3per <- scen2 / sum(pop$pop[[1]]$size[4:16])

sprintf("%s: %f", tariso, scen3per)

stop()

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
