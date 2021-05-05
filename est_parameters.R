suppressPackageStartupMessages({
  require(data.table)
  require(optimization)
  require(doParallel)
})

#' fixed stride of 5, so adjust starting point; TODO: unfix
.stride <- 5
.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper","PAK","0001")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/pops/%s.rds",
  "%s/outputs/r0/%s.rds",
  "%s/inputs/mobility.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/outputs/sample/%s.rds",
  .debug[2], # ZAF
  .debug[3], # the id
  "../covidm",
  "%s/outputs/params/%s_%s.rds"
), .debug[1], .debug[2], .debug[3]) else commandArgs(trailingOnly = TRUE)

fitslc <- seq(as.integer(tail(.args, 3)[1]), by=1, length.out = .stride)
tariso <- tail(.args, 4)[1]

mob <- readRDS(.args[3])[iso3 == tariso]

timings <- readRDS(.args[4])

#' MAGIC NUMBER: average delay from infection -> reporting
rep_delay <- 7

post_contact_reductions <- timings[(period == 1) & (era == "post")][
  mob, on=.(iso3), allow.cartesian = TRUE
][
  between(date, start - rep_delay, end - rep_delay),
  1 - c(
    home = 1, 
    work = prod(work_multiplier)^(1/.N),
    other = prod(other_multiplier)^(1/.N),
    school = mean(school_multiplier)
  )
]

tarwindow <- timings[era == "relaxation", c(start, end)]

case.dt <- readRDS(.args[2])[
  variable == "infections" &
  between(date, tarwindow[1]-6, tarwindow[2]),
  .(croll = median(value)),
  by=.(date)
]

intros.dt <- readRDS(.args[5])[iso3 == tariso][sample %in% fitslc]
bootstrap.dt <- readRDS(.args[6])[sample %in% fitslc][period == 1]

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
params$schedule <- list(
  list(
    parameter = "contact",
    pops = numeric(),
    mode = "multiply",
    values = contact_schedule,
    times = day0 + 0:(length(contact_schedule)-1)
  )
)

#' TODO: fix warning here; providing correct value, however
startpost <- as.integer(timings[era == "transition" & period == 1, start[1]] - day0)
startrelax <- as.integer(timings[period == 2, start[1]] - day0)

tart <- as.numeric(tarwindow - day0)
case.slc <- case.dt[between(date, tarwindow[1], tarwindow[2]), round(croll)]
endrelax <- as.integer(min(timings[era == "relaxation", end], tarwindow[2]) - day0)
params$time1 <- endrelax

#' assert: fIs reduction is related to number of observed cases
#' more cases == higher percent of peak reduction
#' desired constraints:
#'  - TODO: at infinite cases, value is some L <= 1 (<= 100% reduction); CURRENT: L = 1
#'  - at some reference value of cases, has some known reduction value (post-intervention R0)
#'  
#'  use logistic on observed cases to achieve this
#'  think of basic logistic as filling the available intervention space
#'  (1-baseline reduction) = delta
#'  
#'  reduction = 1/(1+exp(-k(cases-cases0)))*delta + baseline reduction
#'  since we know cases & reduction @ the post-intervention R0 estimate 
#'  
#'  ref_red = 1/(1+exp(-k(ref_cases-cases0)))*(1-baseline) + baseline
#'  (ref_red*(1+exp(-k(ref_cases-cases0))) - 1)/(exp(-k(ref_cases-cases0)) = baseline 
#'  

fIs_baseline <- function(
  case0, k, reff, refcase = case.slc[1]
) {
  M <- exp(-k*(refcase-case0))
  if (is.infinite(M)) {
    reff
  } else {
    (reff*(1+M) - 1)/M
  }
}

#' need (reff*(1+exp(-k*(refcase-case0))) - 1)/exp(-k*(refcase-case0)) > epsilon
suggest_k <- function(
  reff, 
  case0 = 10^(mean(log10(range(case.slc)))),
  refcase = case.slc[1]
) -log(1/reff-1)/(refcase-case0) + c(1e-4,-1e-4)[(reff > 0.5)+1]


fIs_amp <- function(
  case0, k, # fit elements
  reff, # sampling element: value of function at css[1]
  css = case.slc, # data element
  baseline = fIs_baseline(case0, k, reff, css[1]) # entailed remaining coefficient
) (1-baseline)/(1+exp(-k*(css-case0))) + baseline

#   case0 = 10^mean(range(log10(case.slc)))

add_fIs <- function(
  case0, k, # fit elements
  fIs_reduction_at_post, # value at post-intervention
  model_t = startpost, d0 = startrelax, df = endrelax
) {
  reds <- fIs_amp(case0, k, fIs_reduction_at_post)
  vals <- c(list(rep(1-fIs_reduction_at_post, 16)), lapply(reds, function(d) rep(1-d, 16)))
  # vals <- c(list(rep(1-fIs_reduction_at_post, 16)), list(rep(0, 16)))
  return(list(
    parameter = "fIs",
    pops = numeric(),
    mode = "multiply",
    values = vals,
    times = day0 + c(model_t, d0:df)
    #times = day0 + c(model_t, model_t+30)
  ))
}

# load covidm
cm_path = tail(.args, 2)[1]
cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T
cm_version = 2

ascll <- function(asc, sim.cases) -sum(dpois(case.slc, sim.cases*asc, log = TRUE))

sim_step <- function(p) cm_simulate(p, 1, model_seed = 42L)$dynamics[
  compartment == "cases",
  .(value = sum(value)), by=t
]

us <- function(sdt.row, umod = sdt.row[1, umod]) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^u_",names(sdt.row))], each = 2)*umod
}

ys <- function(sdt.row) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^y_",names(sdt.row))], each = 2)
}

#' seeds = intros[sample == sid, t]
dtfun <- function(sdt, pars, seeds, post) {
  pop <- pars # copy constructor
  pop$pop[[1]]$y <- ys(sdt)
  pop$pop[[1]]$u <- us(sdt)
  pop$pop[[1]]$seed_times <- seeds
  
  symp <- optimize(function(symp) abs(cm_eigen_ngm(
    pop, contact_reductions = post_contact_reductions,
    # contact_reductions = c(home=0, work=sml, school=lrg, other=sml),
    fIs_reductions = symp
  )$R0/post-1), c(0.1, 0.9))$minimum
  
  ofun <- function(ps) {
    k <- ps[1]; case0 <- ps[2]
    if (fIs_baseline(case0, k, symp) >= 0) {
      #' project and compare to relaxation period
      pop$schedule[[2]] <- add_fIs(case0, k, symp)
      est <- sim_step(pop)[between(t, tart[1], tart[2]), value]
      #' find the best possible ascertainment prob:
      ret <- optimize(ascll, c(1e-6, .2), sim.cases = est)$objective
      if (is.infinite(ret) | is.na(ret)) ret <- .Machine$integer.max
      ret/length(est)
    } else .Machine$integer.max
  }
  
  cmin <- c(10,case.slc[1]+1)[(symp < 0.5)+1]
  cmax <- c(case.slc[1]-1,1e5)[(symp < 0.5)+1]
  cmed <- sqrt(cmin*cmax)
  ks <- sort(suggest_k(symp, c(cmin, cmed, cmax)))
  
  pars_int <- optim(c(
      k = ks[2],
      case0 = cmed
    ),
    ofun,
    lower = c(k = ks[1], case0 = cmin),
    upper = c(k = ks[3], case0 = cmax),
    method = "L-BFGS-B"
    # ,control = list(
    #   nlimit = 300, ac_acc = 1
    # )
  )
  
  #lrg <- pars_int$par[1]; sml <- pars_int$par[2];
  k <- pars_int$par[1]; case0 <- pars_int$par[2];
  #rel_delay <- pars_int$par[4]; relax_frac <- pars_int$par[5]; relax_period <- pars_int$par[6]
  
  pop$schedule[[2]] <- add_fIs(case0, k, symp)# <- scheduler(lrg, sml, symp, rel_delay, relax_frac, relax_period)
  est <- sim_step(pop)[between(t, tart[1], tart[2]), value]
  
  #' find the best possible ascertainment prob:
  asc <- optimize(ascll, c(1e-6, .2), sim.cases = est)$minimum
  
  pars <- c(symp, pars_int$par, asc, fIs_baseline(case0, k, symp, case.slc[1]))
  #names(pars) <- c("large", "small", "sympt", "rel_delay", "rel_frac", "rel_dur", "asc")
  names(pars) <- c("sympt", "k", "case0", "asc", "fIsbaseline")
  as.list(pars)
}

.cl <- makeCluster(getDTthreads())
clusterExport(.cl, ls(), environment())
clusterEvalQ(.cl, { 
  require(data.table)
  require(optimization)
})

span <- nrow(bootstrap.dt)
#' span <- 2

fits.dt <- rbindlist(parLapply(.cl, X = 1:span, function(i) {
  suppressPackageStartupMessages({
    source(file.path(cm_path, "R", "covidm.R"))
  })
  sdt <- bootstrap.dt[i,]
  res <- dtfun(sdt, params, intros[sdt$sample == sid, t], sdt$post)
  res$sample <- sdt$sample
  res
}))

saveRDS(bootstrap.dt[fits.dt, on = .(sample)], tail(.args, 1))

# stop()
# 
# suppressPackageStartupMessages({
#   source(file.path(cm_path, "R", "covidm.R"))
# })
# 
# fits2.dt <- bootstrap.dt[fits.dt, on = .(sample)]
# est <- rbindlist(lapply(1:nrow(fits2.dt), function(i) with(as.list(fits2.dt[i,]), {
#   sdt <- fits2.dt[i]
#   testpop <- params;
#   testpop$pop[[1]]$y <- ys(sdt)
#   testpop$pop[[1]]$u <- us(sdt)
#   testpop$pop[[1]]$seed_times <- intros[i == sid, t]
#   testpop$schedule[[2]] <- add_fIs(sdt$case0, sdt$k, sdt$sympt)
#   cm_simulate(
#     testpop, 1,
#     model_seed = 42L
#   )$dynamics[compartment == "cases", .(rv = sum(value), value = sum(value)*sdt$asc, asc = sdt$asc), by=.(date=t+day0)][, sample := i ]
# })))
# 
# 
# ggplot(est) +
#  aes(date, rv, group = sample) +
#  geom_line(color="red", alpha = 0.2) +
#   geom_line(aes(y=value, alpha = asc)) +
#  # geom_line(
#  #   aes(date, croll),
#  #   data = case.dt,
#  #   color = "grey", inherit.aes = FALSE
#  # ) +
#   geom_line(
#     aes(date, croll),
#     data = case.dt,
#     color = "black", inherit.aes = FALSE
#   ) +
#  scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
#  scale_y_log10() +
#  theme_minimal()
