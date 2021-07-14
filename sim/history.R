suppressPackageStartupMessages({
  require(data.table)
  require(doParallel)
  require(ggplot2)
})

.debug <- c("analysis", "GHA")
.args <- if (interactive()) sprintf(c(
  "%s/gen/intervention_timing/%s.rds",
  "%s/gen/pops/%s.rds",
  "%s/gen/yuqs/%s.rds", # ignored
  "%s/gen/mobility.rds",
  "%s/est/introductions/%s.rds",
  "%s/est/params/%s.rds",
  .debug[2], # PAK
  "covidm",
  "%s/sim/history/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 3)[1]

timings <- readRDS(.args[1])
tarwindow <- timings[era == "relaxation", start]
tarwindow[2] <- timings[era == "pre" & period == 3, start]

intros.dt <- readRDS(.args[5])[iso3 == tariso]
day0 <- as.Date(intros.dt[, min(date)])

popsetup <- function(basep, day0) {
  basep$date0 <- day0
  basep
}
params <- popsetup(readRDS(.args[2]), day0)

intros <- intros.dt[,
  intro.day := as.integer(date - date[1])
][, .(t=Reduce(c, mapply(rep, intro.day, infections, SIMPLIFY = FALSE))), by=.(sid=sample) ]

mob <- readRDS(.args[4])[iso3 == tariso]

contact_schedule <- with(mob[date >= day0], mapply(
  function(work, other, school, home = 1) c(home, work, other, school),
  work = workr, other = otherr, school = school_multiplier,
  SIMPLIFY = FALSE
))

#' has all the yuqs info as well; TODO: stop duplicating this?
fits.dt <- readRDS(.args[6])

#' TODO: fix warning here; providing correct value, however
startpost <- as.integer(timings[era == "transition" & period == 1, start[1]] - day0)
startrelax <- as.integer(timings[period == 2, start[1]] - day0)

tart <- as.numeric(tarwindow - day0)
params$time1 <- tart[2]

alltms <- 0:params$time1
pretms <- (1:startpost)-1
posttms <- startpost:params$time1

params$schedule <- list(
  list(
    parameter = "fIs",
    pops = numeric(),
    mode = "multiply",
#    values = c(lapply(pretms, function(x) rep(1, 16)),lapply(posttms, function(x) rep(0.5, 16))),
#    times = c(pretms, posttms)
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

cm_path = tail(.args, 2)[1]
cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T
cm_version = 2

underlying <- function(
  p, case0, k, asc, startt = startpost, symp
) cm_backend_sample_fit_test(
  R_base_parameters = cm_check_parameters(p),
  posterior = data.frame(
    placeholder0=1, placeholder1=1, placeholder2=1, placeholder3=1,
    case0=case0, k=k, asc=asc, startt=startt, symp=symp
  ),
  n = 1, seed = 42L
)[[1]]

sim_step <- function(
  p, case0, k, asc,
  startt = startpost, symp,
  keepoutcomes = "cases",
  idv = c("t", "group")
) melt(underlying(p, case0, k, asc, startt, symp), id.vars = c("t", "group"), measure.vars = keepoutcomes, variable.name = "compartment")[
  , .(value = sum(value)), by=c(idv, "compartment")
]

us <- function(sdt.row, umod = sdt.row[1, umod]) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^u_",names(sdt.row))], each = 2)*umod
}

ys <- function(sdt.row) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^y_",names(sdt.row))], each = 2)
}

.cl <- makeCluster(getDTthreads())
clusterExport(.cl, ls(), environment())
clusterEvalQ(.cl, { 
  require(data.table)
})

span <- nrow(fits.dt)
#' span <- 2

est <- rbindlist(parLapply(.cl, X = 1:span, function(i) {
  suppressPackageStartupMessages({
    source(file.path(cm_path, "R", "covidm.R"))
  })
  sdt <- fits.dt[i]
  testpop <- params;
  testpop$pop[[1]]$y <- ys(sdt)
  testpop$pop[[1]]$u <- us(sdt)
  testpop$pop[[1]]$seed_times <- intros[i == sid, t]
  res <- sim_step(
    testpop, case0 = sdt$case0, k = sdt$k, asc = sdt$asc, symp = sdt$asc,
    keepoutcomes = c("cases", "death_o", "R")
  )
  res[order(t), .(
    sample = i, date = t + day0, group, compartment,
    rv = value, value = value*sdt$asc, asc = sdt$asc
  )]
}))

saveRDS(est, tail(.args, 1))
