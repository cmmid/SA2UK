
suppressPackageStartupMessages({
  require(data.table)
  require(EpiNow2)
  require(qs)
})

#' TODO best way to get these arguments in? dates could come from intervention timing
.debug <- c("~/Dropbox/SA2UK", "ZAF", "2020-06-01", "2020-10-15")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/epi_data.rds",
  "%s/inputs/pops/%s.rds",
  "%s/inputs/covidm_fit_yu.qs",
  .debug[3],
  .debug[4],
  .debug[2],
  "%s/outputs/relaxation/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

est.window <- 14

tariso <- tail(.args, 2)[1]
window.start <- as.Date(.args[2])
window.end <- as.Date(.args[3])

case.dt <- setkey(
  readRDS(.args[1])[iso3==tariso, .(date, confirm=cases)], date
)[between(date, window.start-est.window, window.end+est.window*2)]

params <- readRDS(.args[4])

yu_fits <- qread(.args[5])[order(ll)]
yu_fits[, eqs := (1:.N)/.N ]
#' using the median yu fits
medyu <- yu_fits[which.max(eqs > 0.5)]
yref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("y_",colnames(medyu))]))
uref <- unname(as.matrix(medyu[, .SD, .SDcols = grep("u_",colnames(medyu))]))
ys <- rep(yref[1, ], each = 2)
us <- rep(uref[1, ], each = 2)

params$pop <- lapply(
  params$pop,
  function(x){
    x$y <- ys
    x$u <- us
    return(x)
  }
)

load("NGM.rda")

# Set up example generation time
generation_time <- as.list(
  EpiNow2::generation_times[disease == "SARS-CoV-2",
  .(mean, mean_sd, sd, sd_sd, max=30)
])

tarmcv <- generation_time$mean_sd/generation_time$mean
tarscv <- generation_time$sd_sd/generation_time$sd
tarcv <- generation_time$sd/generation_time$mean

generation_time$mean <- unname(cm_generation_time(params))
generation_time$mean_sd <- generation_time$mean * tarmcv
generation_time$sd <- generation_time$mean * tarcv
generation_time$sd_sd <- generation_time$sd * tarscv

# Set delays between infection and case report
# (any number of delays can be specifed here)
incubation_period <- as.list(
  EpiNow2::incubation_periods[disease == "SARS-CoV-2",
  .(mean, mean_sd, sd, sd_sd, max=30)
])

re.est <- estimate_infections(
  reported_cases = case.dt,
  generation_time = generation_time,
  delays = delay_opts(incubation_period),
  stan = stan_opts(
    samples = 4e3,
    warmup = 200, 
    cores = 4,
    control = list(adapt_delta = 0.9)
  ),
  horizon = 0, CrIs = c(0.5, 0.95),
  verbose = TRUE
)$summarised[variable == "R"][between(date, window.start, window.end)]

saveRDS(re.est, tail(.args, 1))
