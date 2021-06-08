suppressPackageStartupMessages({
  require(data.table)
  require(countrycode)
  require(wpp2019)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper", "GHA")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/mortality.rds",
  "%s/inputs/fertility.rds",
  "%s/inputs/urbanization.rds",
  "%s/inputs/matrices.rds",
  .debug[2],
  "covidm",
  "%s/inputs/pops/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)
#' @examples 
#' .args <- gsub("ZWE","guineabissau",.args)
#' .args <- gsub("ZWE","palestine",.args)

cm_path = tail(.args, 2)[1]
target = tail(.args, 3)[1]
outfile = tail(.args, 1)

mort.dt <- readRDS(.args[1])[iso3 == target]
fert.dt <- readRDS(.args[2])[iso3 == target]
urb.dt <- readRDS(.args[3])[iso3 == target]
mats <- readRDS(.args[4])[[target]]

cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T;
cm_version = 2;

suppressPackageStartupMessages({
  source(file.path(cm_path, "R", "covidm.R"))
})

matref <- country <- cm_populations[
  country_code == countrycode(target, "iso3c", "iso3n"),
  unique(as.character(name))
]

stopifnot(length(country)==1)

#set up node-runs
#TODO: country speci
probs = fread(
  "Age_low,Age_high,Prop_symptomatic,IFR,Prop_inf_hosp,Prop_inf_critical,Prop_critical_fatal,Prop_noncritical_fatal,Prop_symp_hospitalised,Prop_hospitalised_critical
0,9,0.66,8.59E-05,0.002361009,6.44E-05,0.5,0,0,0.3
10,19,0.66,0.000122561,0.003370421,9.19E-05,0.5,9.47E-04,0.007615301,0.3
20,29,0.66,0.000382331,0.010514103,0.000286748,0.5,0.001005803,0.008086654,0.3
30,39,0.66,0.000851765,0.023423527,0.000638823,0.5,0.001231579,0.009901895,0.3
40,49,0.66,0.001489873,0.0394717,0.001117404,0.5,0.002305449,0.018535807,0.3
50,59,0.66,0.006933589,0.098113786,0.005200192,0.5,0.006754596,0.054306954,0.3
60,69,0.66,0.022120421,0.224965092,0.016590316,0.5,0.018720727,0.150514645,0.3
70,79,0.66,0.059223786,0.362002579,0.04441784,0.5,0.041408882,0.332927412,0.3
80,100,0.66,0.087585558,0.437927788,0.065689168,0.5,0.076818182,0.617618182,0.3"
)

#increase CFR
# cfr_RR <- 1.5
# probs[, Prop_critical_fatal := Prop_critical_fatal * cfr_RR]
# probs[, Prop_noncritical_fatal := Prop_noncritical_fatal * cfr_RR]

#min(1, x) does not work for noncritical_fatal, for some reason
# probs[Prop_critical_fatal > 1, Prop_critical_fatal := 1]
# probs[Prop_noncritical_fatal > 1, Prop_noncritical_fatal := 1]

reformat = function(P, lmic_adjust=TRUE) {
  # no info to re-weight these, so assume 70-74 is like 70-79, and 75+ is like 80+
  if(lmic_adjust){
    P <- P[2:length(P)]
    return (rep(P[1:8], each = 2))
  } else {
    return (c(rep(P[1:7], each = 2),P[8:9])) 
  }
}

P.icu_symp     = reformat(probs[, Prop_symp_hospitalised * Prop_hospitalised_critical], lmic_adjust = FALSE);
P.nonicu_symp  = reformat(probs[, Prop_symp_hospitalised * (1 - Prop_hospitalised_critical)], lmic_adjust = FALSE);
P.death    = reformat(probs[, Prop_noncritical_fatal], lmic_adjust = FALSE);

max_time <- 60
tres <- 0.25

ponset2hosp <- cm_delay_gamma(7, 7, max_time, tres)$p
pignore <- cm_delay_skip(max_time, tres)$p
icustay <- cm_delay_gamma(10, 10, max_time, tres)$p
nonicustay <- cm_delay_gamma(8, 8, max_time, tres)$p
ponset2death <- cm_delay_gamma(22, 22, max_time, tres)$p

cm_multinom_process <- function(
  src, outcomes, delays,
  report = ""
) {
  if ("null" %in% names(outcomes)) {
    if (length(report) != length(outcomes)) report <- rep(report, length(outcomes))
    report[which(names(outcomes)=="null")] <- ""
    if (!("null" %in% names(delays))) {
      delays$null <- c(1, rep(0, length(delays[[1]])-1))
    }
  } else if (!all(rowSums(outcomes)==1)) {
    report <- c(rep(report, length(outcomes)), "")
    outcomes$null <- 1-rowSums(outcomes)
    delays$null <- c(1, rep(0, length(delays[[1]])-1))
  }
  nrow <- length(outcomes)
  list(
    source = src, type="multinomial", names=names(outcomes), report = report,
    prob = t(as.matrix(outcomes)), delays = t(as.matrix(delays))
  )
}

cm_track_process <- function(src, name, delays, agecats = 16, report = "p") {
  list(
    source = src, type="multinomial", names = name, report = report,
    prob = matrix(1, nrow = 1, ncol = agecats),
    delays = t(delays)
  )
}

burden_processes = list(
  # process of sending symptomatic cases to hospital icu or ward
  cm_multinom_process(
    "Ip",
    outcomes = data.frame(to_icu = P.icu_symp, to_nonicu = P.nonicu_symp),
    delays = data.frame(to_icu = ponset2hosp, to_nonicu = ponset2hosp)
  ),
  # track icu prevalance
  cm_track_process("to_icu", "icu", icustay),
  # track ward prevalence
  cm_track_process("to_nonicu", "nonicu", nonicustay),
  # track infections - get from delta R prevalence
#  cm_track_process("S", "infection", pignore, report="i"),
  # send some cases to death, tracking outcidence
  cm_multinom_process(
    "Ip",
    outcomes = data.table(death=P.death),
    delays = data.table(death=ponset2death),
    report = "o"
  )
)

popnorm <- function(x, seed_cases = 50, urbfrac) {
  
  #age-specific probability of being symptomatic
  #x$y <- c(rep(0.056, 3), rep(0.49, 8), rep(0.74, 8))
  #new values proposed by nick
  # x$y <- c(
  #   rep(0.2973718, 2), rep(0.2230287, 2), rep(0.4191036, 2),
  #   rep(0.4445867, 2), rep(0.5635720, 2), rep(0.8169443, 6)
  # )
  
  # no cases in empty compartments
  x$dist_seed_ages <- as.numeric(!(x$size == 0))
  
  # seed cases
  x$seed_times <- rep(0, seed_cases)
  
  # incorporate urbanization fraction(s)
  x$size <- round(x$size*urbfrac)
  
  return(x)
}

calllist <- list(
  country, matref,
  dE  = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 0.25)$p,
  dIp = cm_delay_gamma(1.5, 4.0, t_max = 15, t_step = 0.25)$p,
  dIs = cm_delay_gamma(3.5, 4.0, t_max = 15, t_step = 0.25)$p,
  dIa = cm_delay_gamma(5.0, 4.0, t_max = 15, t_step = 0.25)$p,
  A   = c(rep(1/(5*365.25), 15), 0),
  B   = c(fert.dt$per_capita_day, rep(0, 15)),
  D   = mort.dt$per_capita_day
)

if (!is.null(mats)) calllist$matrices <- mats

params <- cm_base_parameters_SEI3R(
  deterministic=FALSE,
  pop=list(do.call(cm_build_pop_SEI3R, calllist))
)

params$schedule = list()

params$processes = burden_processes

params$pop <- lapply(params$pop, popnorm, urbfrac = urb.dt$value)
#params1$time1 <- as.Date(params1$time1)

saveRDS(params, outfile)
