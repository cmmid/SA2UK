
suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/params/%s_consolidated.rds",
  "%s/inputs/pops/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/inputs/urbanization.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/outputs/variant/%s.rds",
  .debug[2],
  "../covidm",
  "%s/outputs/vax/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 3)[1]

fits <- readRDS(.args[1])
varinfo <- readRDS(.args[6])
intros.dt <- readRDS(.args[3])[iso3 == tariso]
urbfrac <- readRDS(.args[4])[iso3 == tariso, value / 100]
timings <- readRDS(.args[5])

day0 <- as.Date(intros.dt[, min(date)])
intros <- intros.dt[,
   intro.day := as.integer(date - date[1])
][, .(t=Reduce(c, mapply(rep, intro.day, infections, SIMPLIFY = FALSE))), by=sample ]

popsetup <- function(basep, urbanfraction, day0) {
  basep$date0 <- day0
  basep$pop[[1]]$size <- round(basep$pop[[1]]$size*urbanfraction)
  basep$pop[[1]]$dist_seed_ages <- c(rep(0,4), rep(1, 6), rep(0, 6))
  basep
}

base <- popsetup(readRDS(.args[2]), urbfrac, day0)

startrelax <- as.integer(timings[era == "relaxation", start] - day0)
endrelax <- as.integer(timings[era == "relaxation", end] - day0)
startvar <- timings[era == "pre" & period == 3, start]

startpost <- as.integer(timings[era == "transition", start[1]] - day0)

tms <- day0 + startpost
relaxtms <- day0 + startrelax:endrelax

tier2 <- as.Date("2020-08-15")

load("NGM.rda")

# load covidm
cm_path = tail(.args, 2)[1]
cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T
cm_version = 2

suppressPackageStartupMessages({
  source(file.path(cm_path, "R", "covidm.R"))
})

vax_delay <- as.Date("2021-04-01")
t_vax <- as.numeric(vax_delay - day0)
vax_eff <- 0.9
vax_imm_dur_days <- 0
doses_per_day <- 4000
from_age <- 4
to_age <- 16
strategy_str <- 365
t_end <- t_vax + strategy_str
horizon <- 5
base$time1 <- t_vax + horizon*365
strategy <- "campaign"

mk_waning <- function(baseline_dur_days, ages = 16, age_dur_mods = rep(1, ages) ) {
  rep(
    ifelse(baseline_dur_days == 0, 0, 1/baseline_dur_days),
    ages
  ) / age_dur_mods
}

calc_dose_per_by_age <- function(f, t, sizes, total_dpd) {
  doses_per_day <- rep(0, 16)
  tar_ages <- from_age:to_age
  vp <- sizes[tar_ages]; vp <- vp/sum(vp)
  #' TODO potentially make demographic sensitive?
  doses_per_day[tar_ages] <- floor(vp*total_dpd)
  del <- total_dpd - sum(doses_per_day)
  if (del) {
    del_tar <- from_age:(from_age+del-1)
    doses_per_day[del_tar] <- doses_per_day[del_tar] + 1
  }
  doses_per_day
}

vaxschedule <- {if (strategy == "campaign") {
  # set parameters for this set of scenario runs
  base$pop[[1]]$ev = rep(vax_eff, 16) #' TODO mods by age?
  base$pop[[1]]$wv = mk_waning(vax_imm_dur_days)
  
  doses_per_day_base <- calc_dose_per_by_age(from_age, to_age, base$pop[[1]]$size, doses_per_day)
  doses_per_day_later <- calc_dose_per_by_age(4, to_age, base$pop[[1]]$size, doses_per_day)
  
  covaxIncreaseMult <- c(4, 6, 8)
  covaxIncreaseDays <- t_vax + seq(91, by=91, length.out = length(covaxIncreaseMult))
  
  dpd <- lapply(
    1:(length(covaxIncreaseMult)+2),
    function (i) c(1, covaxIncreaseMult, 0)[i]*(if (i==1) doses_per_day_base else doses_per_day_later)
  )
  
  ### NEW BIT IS HERE
  list(list(   # care when adding to schedule
    parameter = 'v',             # impact parameter 'v' (vaccines administered per day for each age group)
    pops = 0,                    # 0th population
    mode = 'assign',             # assign values to v
    times =     c(t_vax, covaxIncreaseDays, t_end) + day0,    # do changes on vax day, vax day + 90
    values = dpd
    # however many doses a day for strategy_str days, then stop
  ))
  
} else {
  list()
}}

vartimes <- seq(as.Date("2020-10-15")-1, startvar+20, by="day")
#' interpolate between 1 and variant_mod over phylo period
#' TODO alternative functional shape?
varmul <- function(vm) lapply(seq(1, vm, along.with = vartimes), function(vm) rep(vm, 16))

stretch <- function(end, along) lapply(along, function(...) end[[1]])

varlater <- seq(tail(vartimes, 1)+1, base$time1 + day0, by=1)
relaxlater <- seq(tail(relaxtms, 1)+1, base$time1 + day0, by=1)

scheduler <- function(large, small, symp, k, shft, variant_mod, vax) {
  cons <- list(1-c(0, small, large, small))
  si <- list(rep(1-symp, 16))
  relaxfact <- 1-(1+exp(-k*as.numeric(relaxtms-tier2-shft)))^-1
  relaxfact <- (1-relaxfact[1]) + relaxfact
  relaxcons <- lapply(relaxfact, function(rf) 1-c(0, small, large, small)*rf)
  relaxsi <- lapply(relaxfact, function(rf) 1-rep(symp, 16)*rf)
  reffact <- rep(1, 16)
  
  vm <- varmul(variant_mod)
  
  ret <- list(
    list(
      parameter = "contact",
      pops = numeric(),
      mode = "multiply",
      values = c(cons, relaxcons, stretch(tail(relaxcons,1), relaxlater)),
      times = c(tms, relaxtms, relaxlater)
    ),
    list(
      parameter = "fIs",
      pops = numeric(),
      mode = "multiply",
      values = c(si, relaxsi, stretch(tail(relaxsi,1), relaxlater)),
      times = c(tms, relaxtms, relaxlater)
    ),
    list(
      parameter = "u",
      pops = numeric(),
      mode = "multiply",
      values = c(vm, stretch(tail(vm, 1), varlater)),
      times = c(vartimes, varlater)
    )
  )
  if (length(vax)) {
    return(c(ret, vax)) 
  } else return(ret)
}

keepoutcomes <- c(
  "cases", "death_o",
  "non_icu_severe_i", "non_icu_critical_i", "icu_critical_i",
  "non_icu_severe_p", "non_icu_critical_p", "icu_critical_p"
)

sims <- fits[, {
  # browser()
  us <- rep(.SD[, as.numeric(.SD), .SDcols = grep("^u_",names(.SD))], each = 2)*umod
  ys <- rep(.SD[, as.numeric(.SD), .SDcols = grep("^y_",names(.SD))], each = 2)
  testpop <- base;
  testpop$pop[[1]]$y <- ys
  testpop$pop[[1]]$u <- testpop$pop[[1]]$u*us
  sid <- sample
  var_mod <- varinfo[sample == sid, withdepl]
  testpop$pop[[1]]$seed_times <- intros[sample == sid, t]
  testpop$schedule <- scheduler(large, small, sympt, k, shft, var_mod, vaxschedule)
  res <- cm_simulate(
    testpop, 1, model_seed = 42L
  )$dynamics[compartment %in% keepoutcomes]
}, by=sample]

res <- sims[, .(
  sample, t,
  group, compartment,
  value
)]

res[order(t), cvalue := cumsum(value), by=.(sample, group, compartment)]

evalts <- t_vax + 0:5*365

int_id <- 2
ret <- res[t %in% evalts][, anni_year := (t - t_vax)/365 ][, date := t + day0 ][, epi_id := 0 ][, intervention_id := int_id ]

saveRDS(ret, sprintf("example_%04i.rds", int_id))

#' @examples 
#' comparison <- res[compartment == "cases" & between(date, tarwindow[1], tarwindow[2]), .(value = sum(value)), by=.(sample, date)][sample == 1]
#' plot.dt <- res[compartment == "cases"][fits[, .(sample, asc)], on=.(sample)][, .(asc.value = sum(value*asc), value = sum(value)), by=.(sample, date)]
#' ggplot(plot.dt) +
#'   aes(date, asc.value, group = sample) +
#'   geom_line(alpha = 0.1) +
#'   geom_line(aes(y=value), alpha = 0.1, color = "red") +
#'   geom_line(
#'     aes(date, cases),
#'     data = readRDS("~/Dropbox/SA2UK/inputs/epi_data.rds")[iso3 == "ZAF"],
#'     color = "black", inherit.aes = FALSE
#'   ) +
#'   annotate("rect", xmin=as.Date("2020-09-01"), xmax =as.Date("2020-10-01"), ymin = 0.01, ymax = Inf, alpha = 0.2, fill = "dodgerblue") +
#'   geom_vline(xintercept = as.Date("2020-11-22"), color = "dodgerblue") +
#'   scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
#'   scale_y_log10() +
#'   theme_minimal() +
#'   coord_cartesian(xlim = as.Date(c(NA, "2021-03-01")))

saveRDS(res, tail(.args, 1))

