
suppressPackageStartupMessages({
  require(data.table)
  require(doParallel)
})

.debug <- c(".", "NGA", "1")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/pops/%s.rds",
  "%s/outputs/params/%s_consolidated.rds",
  "%s/inputs/mobility.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/outputs/variant/%s.rds",
  .debug[3],
  .debug[2], # PAK
  "covidm",
  "%s/outputs/scenarios/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

fits.dt <- readRDS(.args[2])

tariso <- tail(.args, 3)[1]
scenid <- as.integer(tail(.args, 4)[1])

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
  basep$pop[[1]]$dist_seed_ages <- c(rep(0,4), rep(1, 6), rep(0, 6))
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
reft <- as.numeric(as.Date("2021-09-01")-day0)

years <- 5
params$time1 <- reft + 365*years

alltms <- 0:params$time1
pretms <- (1:startpost)-1
posttms <- startpost:params$time1

lastmonth_gm_contacts <- list(Reduce(`*`, contact_schedule[length(contact_schedule)-(0:30)])^(1/30))

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
    #values = contact_schedule,
    values = c(contact_schedule, rep(lastmonth_gm_contacts, length(alltms)-length(contact_schedule))),
    times = alltms
  )
)

vocintro <- function(vocmul, tm = vocday) lapply(c("fIs", "fIa", "fIp"), function(par) list(
  parameter = par,
  pops = numeric(),
  mode = "multiply",
  values = list(rep(vocmul, 16)),
  times = vocday
))

vaccination <- function(
  coverage,
  fromt = reft, duration = 365,
  ages = 4:16, agedistro = params$pop[[1]]$size
) {
  agedistro[-ages] <- 0
  agedistro[ages] <- round(agedistro[ages]*coverage/duration)
  list(
    parameter = 'v', pops = numeric(),
    mode = "assign",
    values = list(agedistro, rep(0, length(agedistro))),
    times = c(fromt, fromt+duration)
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
) {
  res <- underlying(p, case0, k, baseline, asc, startt = startpost)
#  browser()
  melt(res, id.vars = c("t", "group"), measure.vars = keepoutcomes, variable.name = "compartment")[
  , .(value = sum(value)), by=c(idv, "compartment")
  ]
}

us <- function(sdt.row, umod = sdt.row[1, umod]) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^u_",names(sdt.row))], each = 2)*umod
}

ys <- function(sdt.row) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^y_",names(sdt.row))], each = 2)
}

keepers <- c("cases", "death_o", "R", "nonicu_p", "nonicu_i", "icu_p", "icu_i", "Sv")

scen.dt <- rbind(
  data.table(vax_mech = "none", vax_eff = 0, coverage = 0),
  data.table(expand.grid(
    vax_mech = c("infection", "disease"),
    vax_eff = c(.5, .9),
    coverage = c(.25, .90)
  ))
)[, epi_id := (1:.N)-1 ][ epi_id == scenid ]

params$processes[[2]]$report <- "ip"
params$processes[[3]]$report <- "ip"

if (scen.dt$epi_id != 0) {
  params$schedule <- c(params$schedule, list(vaccination(scen.dt$coverage)))
  params$pop[[1]]$ev <- rep(scen.dt$vax_eff, params$pop[[1]]$n_groups)
  if (scen.dt$vax_mech == "infection") params$pop[[1]]$uv <- rep(0, params$pop[[1]]$n_groups)
  if (scen.dt$vax_mech == "disease") params$pop[[1]]$yv <- rep(0, params$pop[[1]]$n_groups)
}

.cl <- makeCluster(getDTthreads()-2)
clusterExport(.cl, ls(), environment())
clusterEvalQ(.cl, { 
  require(data.table)
  require(optimization)
})

span <- nrow(fits.dt)
# span <- 2

#' scenarios
#'  - no change (may not appropriately reflect natural shifts in contact patterns)
#'  - no school
#'  - re-introduce peak workplace reductions
#'  - no school, peak workplace / other reductions
#'  - minimal vaccination scenario - slow dosage, but: high eff, no waning

est <- rbindlist(parLapply(.cl, X = 1:span, function(i) with(scen.dt, {
  suppressPackageStartupMessages({
    source(file.path(cm_path, "R", "covidm.R"))
  })
  sdt <- fits.dt[i]
  testpop <- params;
  testpop$pop[[1]]$y <- ys(sdt)
  if (scen.dt$vax_mech != "disease") testpop$pop[[1]]$yv <- testpop$pop[[1]]$y
  testpop$pop[[1]]$u <- us(sdt)
  if (scen.dt$vax_mech != "infection") testpop$pop[[1]]$uv <- testpop$pop[[1]]$u
  testpop$pop[[1]]$seed_times <- intros[i == sid, t]
  testpop$schedule <- c(testpop$schedule, vocintro(variants.dt[i, withdepl]))
  res <- sim_step(
    testpop, case0 = sdt$case0, k = sdt$k, baseline = sdt$fIsbaseline, asc = sdt$asc,
    keepoutcomes = keepers
  )[order(t), .(
    sample = i, date = t + day0, group, compartment,
    rv = value, value = value*sdt$asc, asc = sdt$asc
  )]
  res
})))[, epi_id := scenid ]

#' @examples
#' ggplot(res[,.(value = sum(rv)), by=.(date, compartment)]) + 
#'   aes(date, value) + 
#'   facet_grid(compartment ~ ., scales = "free_y", switch="y") + 
#'   geom_line() + scale_x_date() + 
#'   theme_minimal() + 
#'   scale_y_log10("count of ...", labels = scales::label_number_si(), minor_breaks = NULL) +
#'   geom_vline(xintercept = day0+vocday, color = "red") + 
#'   geom_vline(xintercept = day0+length(contact_schedule)-1, color="green") +
#'   geom_vline(xintercept = day0+reft, color="blue")


saveRDS(est, tail(.args, 1))
