
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/SA2UK")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/epi_data.rds",
  "%s/outputs/figs/cfr.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

dt <- readRDS(.args[1])[iso3 == "ZAF"]

window <- 7

dt[,c("cases7","deaths7") := .(frollsum(cases, window),frollsum(deaths, window))]

death.delay <- 21

cases <- dt[!is.na(cases7),.(date, cases7)]
deaths <- dt[!is.na(cases7),.(date, deaths7)]
undelaydeaths <- copy(deaths)[, date := date - death.delay ]

cfr_naive <- deaths[cases, on=.(date)][which.max(cases7 > 0):.N]
#' assert: always have more cases than deaths
#' @example 
#' cfr_naive[, all(cases7 >= deaths7)]

cfr_delay <- undelaydeaths[cases, on=.(date)][which.max(cases7 > 0):.N][!is.na(deaths7)]

bino <- function(ci, pos, tot) as.data.table(t(mapply(
  function(x, n, p=x/n) binom.test(x, n, p, conf.level = ci)$conf.int,
  x = pos, n = tot
)))

res.dt <- rbind(
  cfr_naive[, c("cfr","ver") := .(deaths7/cases7, "nCFR")],
  cfr_delay[, c("cfr","ver") := .(deaths7/cases7, "dCFR")]
)

res.dt[, c("lo95","hi95") := bino(.95, deaths7, cases7) ]
res.dt[, c("lo50","hi50") := bino(.50, deaths7, cases7) ]

cfr.p <- force(ggplot(res.dt) + aes(date, cfr) +
  geom_line(aes(color = ver)) +
  geom_ribbon(aes(fill = ver, ymin = lo95, ymax = hi95), alpha = 0.2) +
  geom_ribbon(aes(fill = ver, ymin = lo50, ymax = hi50), alpha = 0.5) +
  coord_cartesian(
    ylim = c(0, .1),
    expand = FALSE
  ) + 
  scale_y_continuous("Case Fatality Rate, CFR", labels = function(v) sprintf("%0.2g%%", v*100)) +
  scale_x_date(name = NULL, date_breaks = "months", date_minor_breaks = "weeks", date_labels = "%b") +
  scale_color_manual(
    "CFR Calculation",
    labels = c(nCFR="naive", dCFR=sprintf("lagged %i days", death.delay)),
    values = c(nCFR="firebrick", dCFR="dodgerblue"),
    aesthetics = c("color", "fill")
  ) + theme_minimal() +
  theme(
    legend.position = c(0.5, 1), legend.justification = c(0.5, 1)
  ))

saveRDS(cfr.p, tail(.args, 1))