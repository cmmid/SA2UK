suppressPackageStartupMessages({
  require(data.table)
  require(doParallel)
  require(ggplot2)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper","PAK")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/pops/%s.rds",
  "%s/outputs/params/%s_consolidated.rds",
  "%s/inputs/mobility.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/outputs/introductions/%s.rds",
  .debug[2], # PAK
  "covidm",
  "%s/outputs/projections/%s.rds"
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
params$time1 <- tart[2]

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

.cl <- makeCluster(getDTthreads())
clusterExport(.cl, ls(), environment())
clusterEvalQ(.cl, { 
  require(data.table)
  require(optimization)
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
    testpop, case0 = sdt$case0, k = sdt$k, baseline = sdt$fIsbaseline, asc = sdt$asc,
    keepoutcomes = c("cases", "death_o", "R")
  )
  res[order(t), .(
    sample = i, date = t + day0, group, compartment,
    rv = value, value = value*sdt$asc, asc = sdt$asc
  )]
}))


# sims <- fits[,{
#   us <- rep(.SD[, as.numeric(.SD), .SDcols = grep("^u_",names(.SD))], each = 2)*umod
#   ys <- rep(.SD[, as.numeric(.SD), .SDcols = grep("^y_",names(.SD))], each = 2)
#   testpop <- base;
#   testpop$pop[[1]]$y <- ys
#   testpop$pop[[1]]$u <- testpop$pop[[1]]$u*us
#   sid <- sample
#   testpop$pop[[1]]$seed_times <- intros[sample == sid, t]
#   testpop$schedule <- scheduler(large, small, sympt, k, shft)
#   res <- cm_simulate(
#     testpop, 1, model_seed = 42L
#   )$dynamics[compartment %in% c("cases","death_o","R")]
# }, by=sample]
# 
# res <- sims[,
#             .(
#               sample, date = t + day0,
#               group, compartment,
#               value
#             )
# ]

#' @examples 
#' comparison <- res[compartment == "cases" & between(date, tarwindow[1], tarwindow[2]), .(value = sum(value)), by=.(sample, date)][sample == 1]
#' ggplot(res[compartment == "cases"][fits[, .(sample, asc)], on=.(sample)][, .(asc.value = sum(value)*asc, value = sum(value)), by=.(sample, date)]) +
#'   aes(date, asc.value, group = sample) +
#'   geom_line(alpha = 0.1) +
#'   geom_line(aes(y=value), alpha = 0.1, color = "red") +
#'   geom_line(
#'     aes(date, cases),
#'     data = readRDS("~/Dropbox/Covid_LMIC/All_Africa_paper/inputs/epi_data.rds")[iso3 == .debug[2]],
#'     color = "black", inherit.aes = FALSE
#'   ) +
#'   annotate("rect", xmin=as.Date("2020-09-01"), xmax =as.Date("2020-10-01"), ymin = 0.01, ymax = Inf, alpha = 0.2, fill = "dodgerblue") +
#'   scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
#'   scale_y_log10() +
#'   theme_minimal()

saveRDS(est, tail(.args, 1))

