
require(data.table)

.debug <- file.path("analysis", "ins") 
.args <- if (interactive()) file.path(
  .debug[1],
  c("vaccines.rds", "campaign_function.rds", "reference_scenarios.rds")
) else commandArgs(trailingOnly = TRUE)

vaxf <- readRDS(.args[1])
campf <- readRDS(.args[2])

vax.products <- vaxf()
camp.spread <- campf()
camp.spread$vaccine <- vax.products
camp.spread$vaccine_dur <- c(1, 2.5, 5, Inf)

scen.dt <- data.table(expand.grid(camp.spread))
scen.dt[, target_id := .GRP, by = doses_per_campaign ]
scen.dt[, recurring_id := .GRP, by = number_of_campaigns ]
scen.dt[, camp_dur_id := .GRP, by = campaign_durations ]
scen.dt[, vax_id := .GRP, by = vaccine ]
scen.dt[, epi_scen_id := 1:.N ]

saveRDS(scen.dt, tail(.args, 1))