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
  data.table(scen_type = "unmitigated"), translated.dt, fill = TRUE
)[, scen_id := 1L:.N ]

#' if there's a modification era
if (fit.dt[era %in% c("modification", "variant"), .N]) {
  mod.dt <- fit.dt[era == "modification"][1]
  var.dt <- fit.dt[era == "variant"][1]
  endday <- as.integer(mod.dt[, date] - day0)
  mod <- copy(translated.dt)[, scen_id := (1L:.N)+1L ]
  mod[, c("self_iso", "school", "home", "work", "other") := NA_real_ ]
  mod[, start_day := endday ]
  all.dt[scen_id > 1, end_day := endday ]
  endday <- as.integer(var.dt[, date] - day0)
  mod[, end_day := endday ]
  modv <- copy(translated.dt)[, scen_id := (1L:.N)+1L ]
  modv[, c("self_iso", "school", "home", "work", "other") := NA_real_ ]
  modv[, start_day := endday ]
  all.dt <- rbind(all.dt, mod, modv)
  # set the end date in the translated scenarios
  # for each of the translated scenarios, run to end date
}

saveRDS(all.dt, outfile)
