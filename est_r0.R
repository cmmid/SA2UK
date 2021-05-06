suppressPackageStartupMessages({
  require(data.table)
  require(EpiNow2)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper", "PAK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/adj_data.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/inputs/yuqs/%s.rds",
  "%s/inputs/pops/%s.rds",
  getDTthreads(),
  "4e3", #' cores, samples
  .debug[2],
  "%s/outputs/r0/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

smps <- as.integer(tail(.args, 3)[1])
crs <- as.integer(tail(.args, 4)[1])

tariso <- tail(.args, 2)[1]

case.dt <- readRDS(.args[1])[iso3 == tariso][, .(date, confirm = cases )]
fill.case <- case.dt[
  case.dt[, .(date = seq(min(date),max(date), by="day"))],
  on=.(date),
  .(date, confirm = fifelse(is.na(confirm), 0, confirm))
]

lims.dt <- readRDS(.args[2])

mean_generation_interval <- readRDS(.args[3])[, mean(si)]

simpar <- readRDS(.args[4])
tstep <- simpar$time_step
pop <- simpar$pop[[1]]

#' Set up example generation time
#' TODO: approach this similarly to incubation period?
#' right now allows variation by age distro etc,
#' but the model assumptions mean that's not practically relevant
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

incubation_period <- estimate_delay(
  sample(length(pop$dE), 100000, replace = T, prob = pop$dE)*tstep
)

est.window <- 30

est.qs <- unique(c(pnorm(seq(-1,0,by=0.25)), pnorm(seq(0,1,by=0.25))))

Rtcalc <- function(case.dt, gp = NULL, rt = rt_opts()) estimate_infections(
  reported_cases = case.dt,
  generation_time = generation_time,
  delays = delay_opts(incubation_period),
  rt = rt,
  stan = stan_opts(
    samples = smps*2,
    warmup = 200, 
    cores = crs,
    control = list(adapt_delta = 0.99, max_treedepth = 20)
  ),
  gp = gp,
  verbose = TRUE,
  CrIs = est.qs,
  horizon = 0
)

processRt <- function(
  rt, keep.start, keep.end,
  era.labels, tarvar = "R"
) rt$samples[variable == tarvar, .(value, variable), by=.(sample, date)][
  between(date, keep.start, keep.end)
][, era := eval(era.labels), by= sample ]

results <- list()

for (grpi in lims.dt[era != "censor", sort(unique(period))]) {
  sublims <- lims.dt[period == grpi]
  daterange <- sublims[, range(c(start, end))]
  incslice <- fill.case[between(date, daterange[1]-est.window, daterange[2]+est.window)]
  breakbased <- sublims[era == "transition", .N]
  if (breakbased) {
    incslice[, era := "tail" ]
    for (e in c("post", "transition", "pre")) {
      incslice[date <= sublims[era == e]$end, era := e  ]
    }
    incslice[date < sublims[era == e]$start, era := "censor" ]
    incslice[, breakpoint := era %in% c("censor", "transition", "tail") ]
    incslice[era == "post", breakpoint := c(TRUE, rep(FALSE, .N-1))]
    results[[grpi]] <- processRt(
      Rtcalc(incslice),
      sublims[era == "pre", end], sublims[era == "post", start],
      expression(c("pre",rep("transition",.N-2),"post"))
    )[, period := grpi ]
  } else {
    results[[grpi]] <- processRt(
      Rtcalc(incslice, rt = NULL, gp = gp_opts()),
      sublims[, min(start)], sublims[, max(end)],
      era.labels = "relaxation", tarvar = "infections"
    )[, period := grpi ]
  }
}

ret <- rbindlist(results)

saveRDS(ret, tail(.args, 1))

# ret <- dcast(results[era != "transition"], sample ~ era, value.var = "value")

#' @examples 
#' require(ggplot2)
#' ggplot(results[era != "variant"][sample %in% sample(.N, 2000)]) + aes(date, value, group = sample) +
#'  geom_line(alpha = 0.05) + theme_minimal()
#' ggplot(results[!(era %in% c("variant", "transition"))][sample %in% sample(.N/2, 100)]) + aes(era, value, group = sample) +
#'  geom_line(alpha = 0.05) + theme_minimal()
#' ggplot(ret[, .(preq=order(pre)/.N, postq=order(post)/.N)]) +
#'  aes(preq, postq) + geom_point() + theme_minimal()
#' ggplot(ret[, .(pre, post)]) +
#'  aes(pre, post) + geom_point() + theme_minimal()
#' ggplot(ret[, .(post, variant)]) +
#'  aes(post, variant) + geom_point() + theme_minimal()

#' smp <- sample(16000, 1000)
#' ext <- 16
#' tst <- results[(sample %in% smp) & (era != "variant")][,
#'   .(
#'     date = c(date[1] - c(ext:1), date, date[.N] + c(1:ext)),
#'     v = (c(rep(value[1], ext), value, rep(value[.N], ext)))
#'   ),
#'   by=sample
#' ][, .(date, v = cumsum((v-1)/generation_time$mean)), by=sample]
#' dr <- tst[, range(date)]
#' epi.dt <- readRDS(
#'   "~/Dropbox/Covid_LMIC/All_Africa_paper/inputs/epi_data.rds")[
#'     iso3 == "ZAF" & between(date, dr[1], dr[2]),
#'     .(date, v = log(cases + exp(tst[date == min(date), mean(v)])), sample = 0)
#' ]
#' p <- ggplot(tst) + aes(date, group = sample) +
#'   geom_line(aes(y=v-1.5), alpha = 0.05) +
#'   geom_line(aes(y=v), data = epi.dt, color = "red") +
#'   theme_minimal()
