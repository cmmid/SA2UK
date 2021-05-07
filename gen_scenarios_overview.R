suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper","PAK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/scenarios/%s.rds",
  .debug[2], # PAK
  "%s/outputs/scenarios/%s.csv"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

est <- readRDS(.args[1])[
  compartment %in% c("cases", "death_o", "icu_p", "nonicu_p"),
  .(rv = sum(rv), value = sum(value)),
  by=.(scenario, sample, date, compartment)
]

est[compartment %in% c("icu_p","nonicu_p"), compartment := "hosp_p" ]

red <- est[between(date, max(date)-364, max(date)), .(rv=sum(rv), value = sum(value)), by=.(scenario, sample, compartment)]
peaks <- est[between(date, max(date)-364, max(date)), .(rv=max(rv), value = max(value)), by=.(scenario, sample, compartment)]


refscen <- red[scenario == "none"]
refscen[, as.list(quantile(rv, probs = c(0.025, 0.5, 0.975))), by=.(compartment)]
refscen[, as.list(quantile(value, probs = c(0.025, 0.5, 0.975))), by=.(compartment)]

refpeaks <- peaks[scenario == "none"]
refpeaks[, as.list(quantile(rv, probs = c(0.025, 0.5, 0.975))), by=.(compartment)]
refpeaks[, as.list(quantile(value, probs = c(0.025, 0.5, 0.975))), by=.(compartment)]

intscen <- red[scenario != "none"][refscen[,.SD,.SDcols=-("scenario")], on=.(sample, compartment)]
intscen[, averted := i.rv - rv ]
intscen[, eff := fifelse(i.rv == 0, 0, averted/i.rv) ]
intscen[, as.list(quantile(averted, probs = c(0.025, 0.5, 0.975))), by=.(scenario, compartment)]
intscen[, as.list(quantile(eff, probs = c(0.025, 0.5, 0.975))), by=.(scenario, compartment)]


intpeaks <- peaks[scenario != "none"][refpeaks[,.SD,.SDcols=-("scenario")], on=.(sample, compartment)]
intpeaks[, averted := i.rv - rv ]
intpeaks[, eff := fifelse(i.rv == 0, 0, averted/i.rv) ]
intpeaks[, as.list(quantile(averted, probs = c(0.025, 0.5, 0.975))), by=.(scenario, compartment)]
intpeaks[, as.list(quantile(eff, probs = c(0.025, 0.5, 0.975))), by=.(scenario, compartment)]

