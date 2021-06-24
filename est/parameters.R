suppressPackageStartupMessages({
  require(data.table)
  require(optimization)
  require(doParallel)
})

.debug <- c("analysis", "NGA")
.args <- if (interactive()) sprintf(c(
  "%s/gen/intervention_timing/%s.rds",
  "%s/gen/pops/%s.rds",
  "%s/gen/yuqs/%s.rds",
  "%s/gen/mobility.rds",
  "%s/est/r0/%s.rds",
  "%s/est/sample/%s.rds",
  "%s/est/introductions/%s.rds",
  "%s/ins/adj_data.rds",
  .debug[2], # ZAF
  "covidm",
  "%s/est/params/%s_%s.rds"
), .debug[1], .debug[2], .debug[3]) else commandArgs(trailingOnly = TRUE)

timings <- readRDS(.args[1])
tariso <- tail(.args, 3)[1]
mob <- readRDS(.args[4])[iso3 == tariso]

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

case.dt <- readRDS(.args[5])[
  variable == "infections",
  .(croll = median(value)),
  by=.(date)
]

adj.case.dt <- readRDS(.args[8])[tariso == iso3 & between(date, tarwindow[1]-rep_delay, tarwindow[2]-rep_delay)]

llsd <- function(sdlog) {
  -sum(dlnorm(adj.case.dt$adj, log(case.dt$croll)-sdlog^2/2, sdlog = sdlog, log = TRUE))
}
sdlog.est <- optimize(llsd, interval = c(0.1, 5))$minimum
ascll <- function(asc, sim.cases) -sum(dlnorm(adj.case.dt$adj, log(sim.cases*asc)-sdlog.est^2/2, sdlog=sdlog.est, log = TRUE))

bootstrap.dt <- readRDS(.args[6])[period == 1]
intros.dt <- readRDS(.args[7])[iso3 == tariso]

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

params <- popsetup(readRDS(.args[2]), day0)

#' TODO: fix warning here; providing correct value, however
startpost <- as.integer(timings[era == "transition" & period == 1, as.Date(start[1])] - day0)
startrelax <- as.integer(timings[period == 2, as.Date(start[1])] - day0)

tart <- as.numeric(tarwindow - day0)
case.slc <- case.dt[between(date, tarwindow[1], tarwindow[2]), round(croll)]
case.ref <- readRDS(.args[8])[tariso == iso3 & date == timings[period == 1 & era == "post", start], adj]
endrelax <- as.integer(min(timings[era == "relaxation", end], tarwindow[2]) - day0)
params$time1 <- endrelax

alltms <- 0:params$time1
pretms <- (1:startpost)-1
posttms <- startpost:params$time1

#' ! n.b.: the fIs multiplier MUST be in the first schedule position
#' TODO: insert this in population construction?
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

#' want:
#'  - target fIs reduction when cases @ the level observed during post-period == sympt estimate
#'  - more reduction when more cases
#'  - less reduction when fewer cases
#' sympt = fIs_reduction => 1-fIs_reduction = multiplier on fIs
#' so fIs = logistic of k > 0 => increasing cases == increasing reduction
#' if post R0 estimate indicates sympt > 0.5, then case0 < reference cases
#' and vice versa

# load covidm
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
  p, case0, k, asc, startt = startpost, sym,
  keepoutcomes = "cases"
) melt(underlying(p, case0, k, asc, startt, symp = sym), id.vars = c("t","group"), measure.vars = keepoutcomes, variable.name = "compartment")[
  , .(value = sum(value)), by=.(t, compartment)
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
    fIs_reductions = symp #+(1-symp)*symp
    #, fIa_reductions = symp, fIp_reductions = symp
  )$R0/post-1), c(0.1, 0.9))$minimum
  
  # cmin <- c(10,case.slc[1]+1)[(symp < 0.5)+1]
  # cmax <- c(case.slc[1]-1,1e5)[(symp < 0.5)+1]
  # cmed <- sqrt(cmin*cmax)
  # ks <- sort(suggest_k(symp, c(cmin, cmed, cmax)))
  
  ofun <- function(ps) {
    k <- exp(ps[1]); case0 <- ps[2]; asc <- exp(ps[3])
    # baseline <- fIs_baseline(case0, k, symp)
    # if (baseline >= 0) {
      #' project and compare to relaxation period
      # pop$schedule[[2]] <- add_fIs(case0, k, symp)
      est <- sim_step(pop, case0, k, asc, sym=symp)[between(t, tart[1], tart[2]), value]
      ret <- ascll(asc, est)
      if (is.infinite(ret) | is.na(ret)) ret <- .Machine$integer.max
      ret
    # } else .Machine$integer.max
  }
  
  kguess <- -6
  case0guess = case.ref + log(1/symp-1)/exp(kguess)
  
  pars_int <- optim(
    par = c(
      k = kguess,
      case0 = case0guess,
      asc = -2
    ),
    fn = ofun,
    #method = "L-BFGS-B",
    lower = c(k = -10, case0 = 0, asc = -6),
    upper = c(k = -2, case0 = case0guess*10, asc = -1)
    # ,control = list(
    #   nlimit = 300, ac_acc = 1
    # )
  )
  
  #lrg <- pars_int$par[1]; sml <- pars_int$par[2];
  k <- exp(pars_int$par[1]); case0 <- pars_int$par[2] # also has asc as par 3

  pars <- c(symp, k, case0, exp(pars_int$par[3]))
  #names(pars) <- c("large", "small", "sympt", "rel_delay", "rel_frac", "rel_dur", "asc")
  names(pars) <- c("sympt", "k", "case0", "asc")
  as.list(pars)
}

.cores <- fcoalesce(
  as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")),
  getDTthreads()
)

.cl <- makeCluster(.cores)
clusterExport(.cl, ls(), environment())
clusterEvalQ(.cl, {
  require(data.table)
  require(optimization)
})

span <- nrow(bootstrap.dt)
#' @example 
#' span <- 3

#' TODO revisit using future?
fits.dt <- bootstrap.dt[rbindlist(parLapply(.cl, X = 1:span, function(i) {
  suppressPackageStartupMessages({
    source(file.path(cm_path, "R", "covidm.R"))
  })
  sdt <- bootstrap.dt[i,]
  res <- dtfun(sdt, pars = params, seeds = intros[sdt$sample == sid, t], post = sdt$post)
  res$sample <- sdt$sample
  res
})), on = .(sample)]

saveRDS(fits.dt, tail(.args, 1))

# stop()
# 
# fits2.dt <- bootstrap.dt[fits.dt, on = .(sample)]
# suppressPackageStartupMessages({
#   source(file.path(cm_path, "R", "covidm.R"))
# })
# 
# fits2.dt <- bootstrap.dt[fits.dt, on = .(sample)]#readRDS("~/Dropbox/Covid_LMIC/All_Africa_paper/outputs/params/PAK_consolidated.rds")[1:10]
# est <- rbindlist(lapply(1:nrow(fits2.dt), function(i) {
#   sdt <- fits2.dt[i]
#   testpop <- params;
#   testpop$pop[[1]]$y <- ys(sdt)
#   testpop$pop[[1]]$u <- us(sdt)
#   testpop$pop[[1]]$seed_times <- intros[i == sid, t]
#   #testpop$schedule[[2]] <- add_fIs(sdt$case0, sdt$k, sdt$sympt)
# #  browser()
# #  thing1 <- underlying(testpop, sdt$case0, sdt$k, sdt$fIsbaseline, sdt$asc)[, .(value = sum(cases), ver = "backend"), by=t]
# #  thing2 <- altunderlying(testpop, sdt$case0, sdt$k, sdt$fIsbaseline, sdt$asc)[, .(value = sum(value), ver = "sim"), by=t]
# #  ggplot(rbind(thing1, thing2)) + aes(t, value, color = ver) + geom_line() + scale_y_log10()
# 
#   ret <- sim_step(testpop, sdt$case0, sdt$k, baseline = sdt$fIsbaseline, sdt$asc)
#   ret[, .(rv = value, value = value*sdt$asc, asc = sdt$asc), by=.(date=t+day0)][, sample := i ]
# }))
# 
# ggplot(est) +
#  aes(date, rv, group = sample) +
#  geom_line(aes(color = factor(sample)), alpha = 0.2) +
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
#   geom_vline(xintercept = day0+startpost, color = "blue") +
#  scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
#  scale_y_log10() +
#  theme_minimal()
1/(1+exp(-0.01*(434-331)))
