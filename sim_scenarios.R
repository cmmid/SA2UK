
suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper", "PAK")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/pops/%s.rds",
  "%s/outputs/params/%s_consolidated.rds",
  "%s/inputs/mobility.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/outputs/variant/%s.rds",
  .debug[2], # PAK
  "covidm",
  "%s/outputs/scenarios/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

fits.dt <- readRDS(.args[2])

tariso <- tail(.args, 3)[1]

mob <- readRDS(.args[3])[iso3 == tariso]
timings <- readRDS(.args[4])
tarwindow <- timings[era == "relaxation", start]
tarwindow[2] <- timings[era == "pre" & period == 3, start]

intros.dt <- readRDS(.args[5])[iso3 == tariso]

day0 <- as.Date(intros.dt[, min(date)])

contact_schedule <- with(mob[date >= day0], mapply(
  function(work, other, school, home = 1) c(home, work, other, school),
  work = workr, other = otherr, school = school_multiplier,
  SIMPLIFY = FALSE
))

intros <- intros.dt[,
                    intro.day := as.integer(date - date[1])
][, .(t=Reduce(c, mapply(rep, intro.day, infections, SIMPLIFY = FALSE))), by=.(sid=sample) ]

popsetup <- function(basep, day0) {
  basep$date0 <- day0
  basep
}

params <- popsetup(readRDS(.args[1]), day0)

#' TODO: fix warning here; providing correct value, however
startpost <- as.integer(timings[era == "transition" & period == 1, start[1]] - day0)
startrelax <- as.integer(timings[period == 2, start[1]] - day0)

tart <- as.numeric(tarwindow - day0)

variants.dt <- readRDS(.args[6])
vocday <- timings[period == 3 & era == "pre", as.numeric(start - day0)]

#' start day of period of interest; in analysis, will work from end
reft <- as.numeric(Sys.Date()-day0)

params$time1 <- reft + 365

alltms <- 0:params$time1
pretms <- (1:startpost)-1
posttms <- startpost:params$time1

params$schedule <- list(
  list(
    parameter = "fIs",
    pops = numeric(),
    mode = "multiply",
    # values = c(lapply(pretms, function(x) rep(1, 16)),lapply(posttms, function(x) rep(0, 16))),
    # times = c(pretms, posttms)
    values = lapply(alltms, function(x) rep(1, 16)),
    times = alltms
  ),
  list(
    parameter = "contact",
    pops = numeric(),
    mode = "multiply",
    values = contact_schedule,
    times = 0:(length(contact_schedule)-1)
  )
)

vocintro <- function(vocmul, tm = vocday) lapply(c("fIs", "fIa", "fIp"), function(par) list(
  parameter = par,
  pops = numeric(),
  mode = "multiply",
  values = list(rep(vocmul, 16)),
  times = vocday
))

redistance <- function(fromt = reft, work = Inf, other = Inf, school = Inf) list(
  parameter = "contact", pops = numeric(),
  mode = "lowerto",
  values = list(c(1, work, other, school)),
  times = fromt
)


noschool <- function(fromt = reft) list(
  parameter = "contact", pops = numeric(),
  mode = "multiply",
  values = list(c(1, 1, 1, 0)),
  times = fromt
)

vaccination <- function(
  dpd = 4000*216.6/47.89*.3691, fromt = reft, ages = 14:16, agedisto = params$pop[[1]]$size
) {
  refpop <- sum(agedisto[ages])
  agedisto[-ages] <- 0
  dpage <- round(agedisto/refpop*dpd)
  list(
    parameter = 'v', pops = numeric(),
    mode = "assign",
    values = list(dpage),
    times = reft
  )
}

#' peak distancing:
peak.dist <- mob[which.min((workr*otherr)^(1/2)), .(workr, otherr)]

cm_path = tail(.args, 2)[1]
cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T
cm_version = 2

underlying <- function(
  p, case0, k, baseline, asc, startt = startpost
) cm_backend_sample_fit_test(
  R_base_parameters = p,
  posterior = data.frame(
    placeholder0=1, placeholder1=1, placeholder2=1, placeholder3=1,
    case0=case0, k=k, baseline=baseline, asc=asc, startt=startt
  ),
  n = 1, seed = 42L
)[[1]]

sim_step <- function(
  p, case0, k, baseline, asc,
  startt = startpost,
  keepoutcomes = "cases",
  idv = c("t", "group")
) melt(underlying(p, case0, k, baseline, asc, startt = startpost), id.vars = c("t", "group"), measure.vars = keepoutcomes, variable.name = "compartment")[
  , .(value = sum(value)), by=c(idv, "compartment")
]

us <- function(sdt.row, umod = sdt.row[1, umod]) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^u_",names(sdt.row))], each = 2)*umod
}

ys <- function(sdt.row) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^y_",names(sdt.row))], each = 2)
}

keepers <- c("cases", "death_o", "R", "nonicu_p", "icu_p")

.cl <- makeCluster(getDTthreads())
clusterExport(.cl, ls(), environment())
clusterEvalQ(.cl, { 
  require(data.table)
  require(optimization)
})

span <- nrow(fits.dt)
#' span <- 2

#' scenarios
#'  - no change (may not appropriately reflect natural shifts in contact patterns)
#'  - no school
#'  - re-introduce peak workplace reductions
#'  - no school, peak workplace / other reductions
#'  - minimal vaccination scenario - slow dosage, but: high eff, no waning

est <- rbindlist(parLapply(.cl, X = 1:span, function(i) {
  suppressPackageStartupMessages({
    source(file.path(cm_path, "R", "covidm.R"))
  })
  sdt <- fits.dt[i]
  testpop <- params;
  testpop$pop[[1]]$y <- ys(sdt)
  testpop$pop[[1]]$u <- us(sdt)
  testpop$pop[[1]]$seed_times <- intros[i == sid, t]
  testpop$schedule <- c(testpop$schedule, vocintro(variants.dt[i, withdepl]))

  scen0 <- sim_step(
    testpop, case0 = sdt$case0, k = sdt$k, baseline = sdt$fIsbaseline, asc = sdt$asc,
    keepoutcomes = keepers
  )[order(t), .(
    sample = i, scenario = "none", date = t + day0, group, compartment,
    rv = value, value = value*sdt$asc, asc = sdt$asc
  )]
  
  testpop$schedule <- c(testpop$schedule, list(noschool()))
  scen1 <- sim_step(
    testpop, case0 = sdt$case0, k = sdt$k, baseline = sdt$fIsbaseline, asc = sdt$asc,
    keepoutcomes = keepers
  )[order(t), .(
    sample = i, scenario = "noschool", date = t + day0, group, compartment,
    rv = value, value = value*sdt$asc, asc = sdt$asc
  )]
  
  testpop$schedule[[length(testpop$schedule)]] <- redistance(work = peak.dist$workr)
  scen2 <- sim_step(
    testpop, case0 = sdt$case0, k = sdt$k, baseline = sdt$fIsbaseline, asc = sdt$asc,
    keepoutcomes = keepers
  )[order(t), .(
    sample = i, scenario = "peakworkred", date = t + day0, group, compartment,
    rv = value, value = value*sdt$asc, asc = sdt$asc
  )]
  
  testpop$schedule[[length(testpop$schedule)]] <- redistance(work = peak.dist$workr, other = peak.dist$otherr, school = 0)
  scen3 <- sim_step(
    testpop, case0 = sdt$case0, k = sdt$k, baseline = sdt$fIsbaseline, asc = sdt$asc,
    keepoutcomes = keepers
  )[order(t), .(
    sample = i, scenario = "peakred", date = t + day0, group, compartment,
    rv = value, value = value*sdt$asc, asc = sdt$asc
  )]
  
  testpop$pop[[1]]$ev = rep(.80, 16) #' TODO mods by age?
  testpop$pop[[1]]$wv = rep(0, 16)
  testpop$schedule[[length(testpop$schedule)]] <- vaccination()
  scen4 <- sim_step(
    testpop, case0 = sdt$case0, k = sdt$k, baseline = sdt$fIsbaseline, asc = sdt$asc,
    keepoutcomes = keepers
  )[order(t), .(
    sample = i, scenario = "vaccination", date = t + day0, group, compartment,
    rv = value, value = value*sdt$asc, asc = sdt$asc
  )]
  
  rbind(scen0, scen1, scen2, scen3, scen4)
  
}))

saveRDS(est, tail(.args, 1))
