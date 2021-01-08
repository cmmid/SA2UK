suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/SA2UK")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/prov_data.rds",
  "%s/outputs/figs/ccfr_supplement.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

# filtering for Western Cape data for main text figure
#dt <- readRDS(.args[1])[province == "WC"]

# unfiltered data for all provinces, for supplementary facet_wrap figure
dt <- readRDS(.args[1])

window <- 21

dt[,
   c("cases.win","deaths.win") := 
  .(frollsum(cases, window), frollsum(deaths, window))
]

delay_distro <- function(mu, sigma) return(
  function(x) plnorm(x+1, mu, sigma) - plnorm(x, mu, sigma)
)

length_out_arg <- 100000

# hospital-to-death (Lognormal truncated) distribution parameter range
# Linton et al. (2020) - lower limit, closest to SA data
mean_interval <- seq(8.7, 20.9, length.out = length_out_arg)
median_interval <- seq(6.7, 13.7, length.out = length_out_arg)

# onset-to-death (Lognormal truncated) distribution parameter range
# Linton et al. (2020) - upper limit
#mean_interval <- seq(15.1, 29.5, length.out = length_out_arg)
#median_interval <- seq(13.5, 24.1, length.out = length_out_arg)

mu_interval <- seq(
  log(min(median_interval)),
  log(max(median_interval)),
  length.out = length_out_arg
)

sigma_interval  <- seq(
  sqrt(2*(log(min(mean_interval)) - min(mu_interval))),
  sqrt(2*(log(max(mean_interval)) - max(mu_interval))),
  length.out = length_out_arg
)

mu_mean <- mean(mu_interval)
mu_sd <- 2*sd(mu_interval)
sigma_mean <- mean(sigma_interval)
sigma_sd <- 2*sd(sigma_interval)


# NEED TO ADD PROVINCE BIT HERE
mu_sampler <- function(n) rnorm(n, mean = mu_mean, sd = mu_sd)
sigma_sampler <- function(n) rnorm(n, mean = sigma_mean, sd = sigma_sd)

bootstrap.dt <- data.table(sample_id = 1:20000)
bootstrap.dt[, c("mu","sigma") := .(mu_sampler(.N), sigma_sampler(.N))]

# function to calculate the adjusted number of 'known outcomes'
# methods from Nishiura et al. (2009)
scale_cfr_rolling <- function(
  cases, deaths, delay_fun
) {
  
  # point estimate of cases with known outcome at time tt
  point_known_t <- numeric(length(cases))
  dlys <- delay_fun((1:length(cases))-1)
  # Sum over cases up to time tt
  for (ii in 1:length(cases)) {
    inds <- 0:(ii-1)
    point_known_t[ii] <- sum(cases[ii - inds] * dlys[inds+1])
  }
  
  # corrected CFR estimator (rolling)
  p_tt_series <- deaths / point_known_t
  
  list(p_tt_series)
}

resbase <- dt[!is.na(cases.win)][order(date)][, .(date, cases.win, deaths.win) ]

corrected <- bootstrap.dt[,{
  dd <- delay_distro(mu, sigma)
  copy(resbase)[, cCFR := scale_cfr_rolling(cases.win, deaths.win, dd)]
}, by=.(sample_id, province)][, {
  qs <- quantile(cCFR, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  names(qs) <- c("lo","md","hi")
  as.list(qs)
}, by=date][!is.na(md)]
corrected[, ver := "cCFR" ]

bino <- function(ci, pos, tot) as.data.table(t(mapply(
  function(x, n, p=x/n) binom.test(x, n, p, conf.level = ci)$conf.int,
  x = pos, n = tot
)))

#naive <- copy(resbase)[cases.win > 0][, md := deaths.win / cases.win ][, ver := "nCFR" ]
#naive[, c("lo","hi") := bino(0.95, deaths.win, cases.win) ]

#deathdelay <- 21

#delayed <- copy(resbase)[which.max(cases.win > 0):.N][,
  #.(date = head(date, -deathdelay), cases.win = head(cases.win, -deathdelay), deaths.win = tail(deaths.win, -deathdelay))
#][, md := deaths.win / cases.win ][, ver := "dCFR" ]
#delayed[, c("lo","hi") := bino(0.95, deaths.win, cases.win) ]

#plot.dt <- rbind(corrected, naive, delayed, fill = TRUE)
plot.dt <- rbind(corrected, fill = TRUE)
plot.dt[, ver := "dCFR"]


# I turned off the naive and delayed CFRs and have added a 
# facet_wrap at the end of the plot over the difference provinces.
# At least, this is how I imagined it working
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
    breaks = c("dCFR"),
    labels = c(dCFR = "delayed"),
    values = c(dCFR="darkorchid4"),
    aesthetics = c("color", "fill")
  ) + theme_minimal()
    + facet_wrap(~province)

saveRDS(cfr.p, tail(.args, 1))
