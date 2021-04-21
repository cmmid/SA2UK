suppressPackageStartupMessages({
  require(data.table)
  require(optimization)
  require(doParallel)
})

#' fixed stride of 20; adjust starting point
.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper","PAK","0001")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/pops/%s.rds",
  "%s/inputs/urbanization.rds",
  "%s/outputs/r0/%s.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/outputs/sample/%s.rds",
  .debug[2], # ZAF
  .debug[3], # the id
  "../covidm",
  "%s/outputs/params/%s_%s.rds"
), .debug[1], .debug[2], .debug[3]) else commandArgs(trailingOnly = TRUE)

fitslc <- seq(as.integer(tail(.args, 3)[1]), by=1, length.out = 20)
tariso <- tail(.args, 4)[1]

urbfrac <- readRDS(.args[2])[iso3 == tariso, value / 100]

timings <- readRDS(.args[4])
tarwindow <- timings[era == "relaxation", c(start, end)]

# case.dt <- readRDS(.args[3])[
#   iso3 == tariso & between(date, tarwindow[1]-6, tarwindow[2]),
#   .(date, croll = frollmean(cases, 7))
# ][!is.na(croll)]

case.dt <- readRDS(.args[3])[
  variable == "infections" &
  between(date, tarwindow[1]-6, tarwindow[2]),
  .(croll = median(value)),
  by=.(date)
]

intros.dt <- readRDS(.args[5])[iso3 == tariso][sample %in% fitslc]
bootstrap.dt <- readRDS(.args[6])[sample %in% fitslc][period == 1]

day0 <- as.Date(intros.dt[, min(date)])
intros <- intros.dt[,
  intro.day := as.integer(date - date[1])
][, .(t=Reduce(c, mapply(rep, intro.day, infections, SIMPLIFY = FALSE))), by=.(sid=sample) ]

popsetup <- function(basep, day0) {
  basep$date0 <- day0
  basep$pop[[1]]$dist_seed_ages <- c(rep(0,4), rep(1, 6), rep(0, 6))
  basep
}

params <- popsetup(readRDS(.args[1]), day0)

tart <- as.numeric(tarwindow - day0)
case.slc <- case.dt[between(date, tarwindow[1], tarwindow[2]), round(croll)]

startrelax <- as.integer(timings[era == "relaxation", start] - day0)
endrelax <- as.integer(min(timings[era == "relaxation", end], tarwindow[2]) - day0)

startpost <- as.integer(timings[era == "transition" & period == 1, start[1]] - day0)

params$time1 <- endrelax

tms <- day0 + startpost
relaxtms <- day0 + startrelax:endrelax

# load covidm
cm_path = tail(.args, 2)[1]
cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T
cm_version = 2

#' reference for all bootstrap evaluation
scheduler <- function(
  large, small, symp,
  relax_delay, # how long until relaxation starts
  relax_fraction, # what will be the total relaxation (% of original value)
  relax_period # how long does relaxation take?
#  k, shft
) {
  cons <- list(1-c(0, small, large, small))
  si <- list(rep(1-symp, 16))
  
  redtime <- relaxtms[-(1:round(relax_delay))]
  
  relaxfact <- 1-(1+exp(-k*as.numeric(relaxtms-tier2-shft)))^-1
  relaxfact <- (1-relaxfact[1]) + relaxfact
  relaxcons <- lapply(relaxfact, function(rf) 1-c(0, small, large, small)*rf)
#  relaxsi <- lapply(relaxfact, function(rf) 1-rep(symp, 16)*rf)
  
  list(
    list(
      parameter = "contact",
      pops = numeric(),
      mode = "multiply",
      values = c(cons, relaxcons),
      times = c(tms, relaxtms)
    ),
    list(
      parameter = "fIs",
      pops = numeric(),
      mode = "multiply",
      values = si, #c(si, relaxsi),
      times = tms #c(tms, relaxtms)
    )
  )
}

ascll <- function(asc, sim.cases) -sum(dbinom(case.slc, sim.cases, asc, log = TRUE))

#' seeds = intros[sample == sid, t]
dtfun <- function(sdt, umod, pars, seeds, post) {
  us <- rep(sdt[, as.numeric(.SD), .SDcols = grep("^u_",names(sdt))], each = 2)*umod
  ys <- rep(sdt[, as.numeric(.SD), .SDcols = grep("^y_",names(sdt))], each = 2)
  pop <- pars # copy constructor
  pop$pop[[1]]$y <- ys
  pop$pop[[1]]$u <- pop$pop[[1]]$u*us
  pop$pop[[1]]$seed_times <- seeds
  pars_int <- optim(
    ofun, pop = pop, post = post,
    par = c(0.8, 0.25, 0.25, 0.01, 0) ,
    lower = c(0.1, 0.01, 0.01, 1e-6, -50),
    upper = c(0.9, 0.9, 0.9, 0.9, 50),
    method = "L-BFGS-B"
  )
  
  lrg <- pars_int$par[1]; sml <- pars_int$par[2]; symp <- pars_int$par[3];
  k <- pars_int$par[4]; shft <- pars_int$par[5]
  
  pop$schedule <- scheduler(lrg, sml, symp, k, shft)
  est <- cm_simulate(
    pop, 1,
    model_seed = 42L
  )$dynamics[
    compartment == "cases",
    .(value = sum(value)), by=t
  ][between(t, tart[1], tart[2]), value]
  
  #' find the best possible ascertainment prob:
  asc <- optimize(ascll, c(1e-6, .99), sim.cases = est)$minimum
  
  pars <- c(pars_int$par, asc)
  names(pars) <- c("large", "small", "sympt", "k", "shft", "asc")
  as.list(pars)
}

ofun <- function(ps, pop, post) {
  lrg <- ps[1]; sml <- ps[2]; symp <- ps[3]; k <- ps[4]; shft <- ps[5]
  if ((lrg < sml) | (lrg < symp)) .Machine$integer.max else {
    #' calculate reduced Rt
    rerr <- (cm_eigen_ngm(
      pop, contact_reductions = c(home=0, work=sml, school=lrg, other=sml),
      fIs_reductions = symp
    )$R0/post-1)^2
    
    #' project and compare to relaxation period
    pop$schedule <- scheduler(lrg, sml, symp, k, shft)
    est <- cm_simulate(
      pop, 1,
      model_seed = 42L
    )$dynamics[
      compartment == "cases",
      .(value = sum(value)), by=t
    ][between(t, tart[1], tart[2]), value]
    #' find the best possible ascertainment prob:
    ret <- optimize(ascll, c(1e-6, .99), sim.cases = est)$objective
    if (is.infinite(ret) | is.na(ret)) ret <- .Machine$integer.max
    rerr + ret/length(est)
  }
}

.cl <- makeCluster(getDTthreads())
clusterExport(.cl, ls(), environment())
clusterEvalQ(.cl, { 
  require(data.table)
  require(optimization)
})

fits.dt <- rbindlist(parLapply(.cl, X = 1:nrow(bootstrap.dt), function(i) {
  suppressPackageStartupMessages({
    source(file.path(cm_path, "R", "covidm.R"))
  })
  sdt <- bootstrap.dt[i,]
  dtfun(sdt, sdt$umod, params, intros[sdt$sample == sid, t], sdt$post)
}))[, sample := 1:.N ]

saveRDS(bootstrap.dt[fits.dt, on = .(sample)], tail(.args, 1))

# suppressPackageStartupMessages({
#   source(file.path(cm_path, "R", "covidm.R"))
# })
# 
# fits2.dt <- bootstrap.dt[fits.dt, on = .(sample)]
# est <- rbindlist(lapply(1:nrow(fits2.dt), function(i) with(as.list(fits2.dt[i,]), {
#   us <- rep(bootstrap.dt[i, as.numeric(.SD)*umod, .SDcols = grep("^u_",names(bootstrap.dt))], each = 2)
#   ys <- rep(bootstrap.dt[i, as.numeric(.SD), .SDcols = grep("^y_",names(bootstrap.dt))], each = 2)
#   testpop <- params; testpop$pop[[1]]$y <- ys
#   testpop$pop[[1]]$u <- testpop$pop[[1]]$u*us
#   testpop$schedule <- scheduler(large, small, sympt, k, shft)
#   cm_simulate(
#     testpop, 1,
#     model_seed = 42L
#   )$dynamics[compartment == "cases", .(value = sum(value)*asc), by=.(date=t+day0)][, sample := i ]
# })))
# 
# 
# ggplot(est) +
#  aes(date, value, group = sample) +
#  geom_line(color="red", alpha = 0.1) +
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
