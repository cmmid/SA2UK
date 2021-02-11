suppressPackageStartupMessages({
  require(EpiNow2)
  require(data.table)
  require(qs)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/epi_data.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/inputs/yuqs/%s.rds",
  "2", "8e3", # cores, samples
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

mean_generation_interval <- readRDS(.args[3])[, mean(si)]

# Set up example generation time
generation_time <- as.list(EpiNow2::generation_times[disease == "SARS-CoV-2",
  .(mean, mean_sd, sd, sd_sd, max=30)
])

tarmcv <- generation_time$mean_sd/generation_time$mean
tarscv <- generation_time$sd_sd/generation_time$sd
tarcv <- generation_time$sd/generation_time$mean

generation_time$mean <- mean_generation_interval
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

est.qs <- unique(c(pnorm(seq(-1,0,by=0.25)), pnorm(seq(0,1,by=0.25))))

Rtcalc <- function(case.dt) estimate_infections(
  reported_cases = case.dt,
  generation_time = generation_time,
  delays = delay_opts(incubation_period),
  #rt = NULL, backcalc = backcalc_opts(),
  stan = stan_opts(
    samples = smps*2,
    warmup = 200, 
    cores = crs,
    control = list(adapt_delta = 0.9, max_treedepth = 20)
  ),
  gp = NULL,
  verbose = TRUE,
  CrIs = est.qs
)

processRt <- function(
  rt, keep.start, keep.end,
  era.labels
) rt$samples[variable == "R", .(value), by=.(sample, date)][
  between(date, keep.start, keep.end)
][, era := eval(era.labels), by= sample ]


results <- processRt(Rtcalc(early_reported_cases),
  lims.dt[era == "pre", end], lims.dt[era == "post", start],
  expression(c("pre",rep("transition",.N-2),"post"))
)

chk <- results[era != "transition",.(med=median(value)), by=era]

if (chk[, sign(diff(med)) != -1]) warning(sprintf("did not observe post-intervention reduction for %s", tariso))

#' if we're considering a modification period as well
#' N.B. a relaxation era is evaluated differently
if (lims.dt[era == "modification", .N]) {
  mod_reported_cases <- with(lims.dt[era == "modification"], fill.case[between(date, start - 14, end + est.window)])
  mod_reported_cases[, breakpoint := TRUE ]
  with(lims.dt[era == "modification"], mod_reported_cases[between(date, start, end), breakpoint := FALSE ])
  results <- rbind(
    results,
    Rtcalc(mod_reported_cases, lims.dt[era == "modification", start], lims.dt[era == "modification", start], "modification")
  )
}

if (lims.dt[era == "variant", .N]) {
  mod_reported_cases <- with(lims.dt[era == "variant"], fill.case[date > start - 14])
  mod_reported_cases[, breakpoint := TRUE ]
  with(lims.dt[era == "variant"], mod_reported_cases[between(date, start, end), breakpoint := FALSE ])
  results <- rbind(
    results,
    processRt(Rtcalc(mod_reported_cases), lims.dt[era == "variant", start], lims.dt[era == "variant", start], "variant")
  )
}

ret <- dcast(results[era != "transition"], sample ~ era, value.var = "value")

#' @examples 
#' require(ggplot2)
#' ggplot(results[era != "variant"]) + aes(date, value, group = sample) +
#'  geom_line(alpha = 0.05) + theme_minimal()
#' ggplot(results[!(era %in% c("variant", "transition"))][sample %in% sample(.N/2, 100)]) + aes(era, value, group = sample) +
#'  geom_line(alpha = 0.05) + theme_minimal()
#' ggplot(ret[, .(preq=order(pre)/.N, postq=order(post)/.N)]) +
#'  aes(preq, postq) + geom_point() + theme_minimal()
#' ggplot(ret[, .(pre, post)]) +
#'  aes(pre, post) + geom_point() + theme_minimal()
#' ggplot(ret[, .(post, variant)]) +
#'  aes(post, variant) + geom_point() + theme_minimal()

saveRDS(ret, tail(.args, 1))

