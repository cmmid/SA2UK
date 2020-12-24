suppressPackageStartupMessages({
  require(EpiNow2)
  require(data.table)
  require(qs)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/epi_data.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/inputs/pops/%s.rds",
  "%s/inputs/covidm_fit_yu.qs",
  "4", "8e3", # cores, samples
  .debug[2],
  "%s/outputs/r0/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

smps <- as.integer(tail(.args, 3)[1])
crs <- as.integer(tail(.args, 4)[1])

tariso <- tail(.args, 2)[1]

case.dt <- readRDS(.args[1])[iso3 == tariso][, .(date, confirm = cases )]
fill.case <- case.dt[
  case.dt[, .(date = seq(min(date),max(date),by="day"))],
  on=.(date),
  .(date, confirm = fifelse(is.na(confirm), 0, confirm))
  ]

lims.dt <- readRDS(.args[2])

params <- readRDS(.args[3])

yu_fits <- qread(.args[4])[order(ll)]
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
generation_time <- as.list(EpiNow2::generation_times[disease == "SARS-CoV-2",
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
incubation_period <- as.list(EpiNow2::incubation_periods[disease == "SARS-CoV-2",
  .(mean, mean_sd, sd, sd_sd, max=30)
])
# replace mean & sd here with what go into rlnorm meanlog, sdlog
# which is not mean(data), sd(data)

# additional time to include for algorithm
est.window <- 30

early_reported_cases <- fill.case[date <= (lims.dt[era == "post"]$end + est.window)]
early_reported_cases[, era := "tail"]
for (e in c("post", "transition", "pre", "censor")) {
  early_reported_cases[date <= lims.dt[era == e]$end, era := e  ]
}

# Add breakpoints
early_reported_cases[,
  breakpoint := era %in% c("censor", "transition", "tail")
]

re.est <- estimate_infections(
  reported_cases = early_reported_cases,
  generation_time = generation_time,
  delays = delay_opts(incubation_period),
  #rt = NULL, backcalc = backcalc_opts(),
  stan = stan_opts(
    samples = smps,
    warmup = 200, 
    cores = crs,
    control = list(adapt_delta = 0.9)
  ),
  gp = NULL,
  verbose = TRUE
)

results <- re.est$samples[variable == "R", .(value), by=.(sample, date)][
  between(date, lims.dt[era == "pre", end], lims.dt[era == "post", start])
][, {
  qs <- quantile(value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
  as.list(qs)
}, keyby = .(date)][, era:= c("pre",rep("transition",.N-2),"post")]

if (results[era %in% c("pre","post"), sign(diff(med)) != -1]) warning(sprintf("did not observe post-intervention reduction for %s", tariso))

#' if we're considering a modification period as well
if (lims.dt[era == "modification", .N]) {
  mod_reported_cases <- with(lims.dt[era == "modification"], fill.case[between(date, start - 14, end + est.window)])
  mod_reported_cases[, breakpoint := TRUE ]
  with(lims.dt[era == "modification"], mod_reported_cases[between(date, start, end), breakpoint := FALSE ])
  mod.est <- estimate_infections(
    reported_cases = mod_reported_cases,
    generation_time = generation_time,
    delays = delay_opts(incubation_period),
 #   rt = NULL, backcalc = backcalc_opts(),
    stan = stan_opts(
      samples = smps,
      warmup = 200, 
      cores = crs,
      control = list(adapt_delta = 0.9)
    ),
    gp = NULL,
    verbose = TRUE
  )
  results <- rbind(
    results,
    mod.est$samples[variable == "R", .(value), by=.(sample, date)][date == lims.dt[era == "modification", start]][, {
      qs <- quantile(value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
      names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
      as.list(qs)
    }, keyby = .(date)][, era:= "modification" ]
  )
}

if (lims.dt[era == "variant", .N]) {
  mod_reported_cases <- with(lims.dt[era == "variant"], fill.case[date > start - 14])
  mod_reported_cases[, breakpoint := TRUE ]
  with(lims.dt[era == "variant"], mod_reported_cases[between(date, start, end), breakpoint := FALSE ])
  mod.est <- estimate_infections(
    reported_cases = mod_reported_cases,
    generation_time = generation_time,
    delays = delay_opts(incubation_period),
    #   rt = NULL, backcalc = backcalc_opts(),
    stan = stan_opts(
      samples = smps,
      warmup = 200, 
      cores = crs,
      control = list(adapt_delta = 0.9)
    ),
    gp = NULL,
    verbose = TRUE
  )
  results <- rbind(
    results,
    mod.est$samples[variable == "R", .(value), by=.(sample, date)][date == lims.dt[era == "variant", start]][, {
      qs <- quantile(value, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
      names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
      as.list(qs)
    }, keyby = .(date)][, era:= "variant" ]
  )
}

saveRDS(results, tail(.args, 1))
