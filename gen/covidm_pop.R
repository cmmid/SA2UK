suppressPackageStartupMessages({
  require(data.table)
  require(countrycode)
  require(wpp2019)
})

.debug <- c("analysis", "GHA")
.args <- if (interactive()) sprintf(c(
  "%s/gen/mortality.rds",
  "%s/gen/fertility.rds",
  "%s/gen/urbanization.rds",
  "%s/gen/matrices.rds",
  .debug[2],
  "covidm",
  "%s/gen/pops/%s.rds"
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

#' IFR from https://doi.org/10.1007/s10654-020-00698-1
#' in that pub
#' log10(IFR in %) = -3.27 + 0.0524*age
#' IFR * 10^2 = 10^(-3.27)*10^(0.0524*age)
#' changing to
#' logit(IFR) basis
ifr_levin = function(age) exp(-7.56 + 0.121 * age) / (100 + exp(-7.56 + 0.121 * age))

#' Infection hospitalisation rate (derived from https://doi.org/10.1126/science.abc3517)
ihr_salje = function(age) exp(-7.37 + 0.068 * age) / (1 + exp(-7.37 + 0.068 * age));

#' Probability of ICU given hospitalisation (derived from CO-CIN)
#' TODO CO-CIN DOI?
hicu_cocin = function(
  age,
  x = c(-0.1309118, 0, 17.2398874, 65.7016492, 100),
  y = c(-2.1825091, -2.1407043, -1.3993552, -1.2344361, -8.8191062)
) {
  inn <- exp(splinefun(x, y)(age))
  inn / (1 + inn)
}

#' using from INFECTION proportions, rather than from ONSET proportions.
#' Means processes source from E (or Ev, Ev2, etc) rather than Ip/Is
#' This is balancing what data are available

relpop <- rbindlist(lapply(c("popF", "popM"), function(pg) as.data.table(get(data(list=pg)))[, iso3 := countrycode(country_code, "iso3n", "iso3c") ][
  !is.na(iso3)
][
  iso3 == target, .(age, `2020`)
]))[, .(pop = sum(`2020`)), by=age ]

groupwidth <- 5
relpop[age != "100+", c("ihr", "ifr", "hicu") := {
  refage <- (.GRP-1)*groupwidth
  .(
    sum(ihr_salje(refage+0:(groupwidth-1)))/groupwidth,
    sum(ifr_levin(refage+0:(groupwidth-1)))/groupwidth,
    sum(hicu_cocin(refage+0:(groupwidth-1)))/groupwidth
  )}, by=age]
#' special case
relpop[age == "100+", c("ihr", "ifr", "hicu") := .(
  ihr_salje(100), ifr_levin(100), hicu_cocin(100)
) ]

maxgroups <- 16
relpop[, group := fifelse(.I < maxgroups, .I, maxgroups) ]

processrates <- relpop[,
  .(
    ihr = sum(ihr*pop)/sum(pop),
    ifr = sum(ifr*pop)/sum(pop),
    hicu = sum(hicu*pop)/sum(pop)
  ),
  by = group
]

P.non_hosp   = 1-processrates$ihr
P.death      = processrates$ifr
P.severe     = processrates$ihr*(1-processrates$hicu) 
P.critical   = processrates$ihr*processrates$hicu

max_time <- 60
tres <- 0.25

pignore <- cm_delay_skip(max_time, tres)$p
#' updated to those used in Sindh analysis
pinf2death <- cm_delay_gamma(26, 5, max_time, tres)$p
pinf2severe <- cm_delay_gamma(8.5, 5, max_time, tres)$p
pinf2critical <- cm_delay_gamma(8.5, 5, max_time, tres)$p

preicustay <- cm_delay_gamma(6.0, 5, max_time, tres)$p
icustay <- cm_delay_gamma(10, 10, max_time, tres)$p
nonicustay <- cm_delay_gamma(9.6, 5, max_time, tres)$p

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

cm_track_process <- function(src, name, delays, agecats = 16, report = "pi") {
  list(
    source = src, type="multinomial", names = name, report = report,
    prob = matrix(1, nrow = 1, ncol = agecats),
    delays = t(delays)
  )
}

burden_processes = list(
  # infections that will ultimate end up in hospital or not
  cm_multinom_process(
    "E",
    outcomes = data.frame(
      to_icu = P.critical,
      to_nonicu = P.severe
    ),
    delays = data.frame(
      to_icu = pinf2critical,
      to_nonicu = pinf2severe
    )
  ),
  # track icu prevalance
  cm_track_process("to_icu", "pre_icu", preicustay),
  cm_track_process("pre_icu", "icu", icustay),
  # track ward prevalence
  cm_track_process("to_nonicu", "nonicu", nonicustay),
  # track infections - get from delta R prevalence
#  cm_track_process("S", "infection", pignore, report="i"),
  # send some cases to death, tracking outcidence
  cm_multinom_process(
    "E",
    outcomes = data.table(death=P.death),
    delays = data.table(death=pinf2death),
    report = "o"
  )
)

popnorm <- function(x, seed_cases = 50, urbfrac, introages = 4:13) {
  
  #' incorporate urbanization fraction(s)
  x$size <- round(x$size*urbfrac)

  #' distribute proportionally to size, though for (15, 65) only [by default]
  x$dist_seed_ages <- rep(0, length(x$size))
  x$dist_seed_ages[introages] <- x$size[introages]/sum(x$size[introages])
  
  #' seed cases
  x$seed_times <- rep(0, seed_cases)
  
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

