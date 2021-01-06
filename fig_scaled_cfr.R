
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/SA2UK")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/epi_data.rds",
  "%s/outputs/figs/ccfr.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

dt <- readRDS(.args[1])[iso3 == "ZAF"]

window <- 21

dt[,
   c("cases.win","deaths.win") := 
  .(frollsum(cases, window), frollsum(deaths, window))
]

dd_mean_mid <- 13.0
dd_median_mid <- 9.1

dd_mu_mid <- log(dd_median_mid)
dd_sigma_mid <- sqrt(2 * (log(dd_mean_mid) - dd_mu_mid))

dd_mean_low <- 8.7
dd_median_low <- 6.7

dd_mu_low <- log(dd_median_low)
dd_sigma_low <- sqrt(2 * (log(dd_mean_low) - dd_mu_low))

dd_mean_high <- 20.9
dd_median_high <- 13.7

dd_mu_high <- log(dd_median_high)
dd_sigma_high <- sqrt(2 * (log(dd_mean_high) - dd_mu_low))

delay_distro <- function(mu, sigma) return(
  function(x) plnorm(x+1, mu, sigma) - plnorm(x, mu, sigma)
)

hosp_to_death_mid <- delay_distro(dd_mu_mid, dd_sigma_mid)
hosp_to_death_low <- delay_distro(dd_mu_low, dd_sigma_low)
hosp_to_death_high <- delay_distro(dd_mu_high, dd_sigma_high)

# function to calculate the adjusted number of 'known outcomes'
# methods from Nishiura et al. (2009)
scale_cfr_rolling <- function(
  cases, deaths, delay_fun
) {
  
  point_known_t <- NULL # point estimate of cases with known outcome at time tt
  
  # Sum over cases up to time tt
  for (ii in 1:length(cases)) {
    
    known_i <- 0 # number of cases with known outcome at time ii
    
    for (jj in 0:(ii - 1)) {
      known_jj <- (cases[ii - jj] * delay_fun(jj))
      known_i <- known_i + known_jj
    }
    
    # point estimate of known outcomes
    point_known_t <- c(point_known_t, known_i)
  }
  
  # corrected CFR estimator (rolling)
  p_tt_series <- deaths / point_known_t
  
  list(p_tt_series)
}

resbase <- dt[!is.na(cases.win)][order(date)][, .(date, cases.win, deaths.win) ]

corrected <- dcast(rbind(
  copy(resbase)[, q := "lo" ][, cCFR := scale_cfr_rolling(cases.win, deaths.win, hosp_to_death_low)],
  copy(resbase)[, q := "md" ][, cCFR := scale_cfr_rolling(cases.win, deaths.win, hosp_to_death_mid)],
  copy(resbase)[, q := "hi" ][, cCFR := scale_cfr_rolling(cases.win, deaths.win, hosp_to_death_high)]
), date ~ q, value.var = "cCFR")[!is.na(md)]

corrected[, hi := pmax(hi, lo, md) ][, lo := pmin(hi, lo, md) ][, ver := "cCFR" ]

bino <- function(ci, pos, tot) as.data.table(t(mapply(
  function(x, n, p=x/n) binom.test(x, n, p, conf.level = ci)$conf.int,
  x = pos, n = tot
)))

naive <- copy(resbase)[cases.win > 0][, md := deaths.win / cases.win ][, ver := "nCFR" ]
naive[, c("lo","hi") := bino(0.95, deaths.win, cases.win) ]

plot.dt <- rbind(corrected, naive, fill = TRUE)

cfr.p <- force(ggplot(plot.dt) + aes(date, md) +
  geom_line(aes(color = ver)) +
  geom_ribbon(aes(fill = ver, ymin = lo, ymax = hi), alpha = 0.2) +
#  geom_ribbon(aes(fill = ver, ymin = lo50, ymax = hi50), alpha = 0.5) +
  coord_cartesian(
    ylim = c(0, .075),
    xlim = as.Date(c("2020-04-01", "2021-01-01")),
    expand = FALSE
  ) + 
  scale_y_continuous(
    "Case Fatality Ratio (CFR)",
    breaks = c(0,0.025,0.05,0.075),
    labels = function(v) sprintf("%0.2g%%", v*100)
  ) +
  scale_x_date(name = NULL, date_breaks = "months", date_minor_breaks = "weeks", date_labels = "%b") +
  scale_color_manual(
    NULL,
    labels = c(nCFR="naive", cCFR="corrected"),
    values = c(nCFR="orchid", cCFR="darkorchid4"),
    aesthetics = c("color", "fill")
  ) + theme_minimal())

saveRDS(cfr.p, tail(.args, 1))




