suppressPackageStartupMessages({
  require(data.table)
  require(optimization)
})
if (sys.nframe() == 0) {
    #' fixed stride of 20; adjust starting point
  .debug <- c("~/Dropbox/SA2UK", "ZAF", "0001")
  .args <- if (interactive()) sprintf(c(
    "%s/inputs/pops/%s.rds",
    "%s/inputs/urbanization.rds",
    "%s/inputs/epi_data.rds",
    "%s/outputs/intervention_timing/%s.rds",
    "%s/outputs/introductions/%s.rds",
    "%s/outputs/sample/%s.rds",
    .debug[2], # ZAF
    .debug[3], # the id
    "../covidm",
    "%s/outputs/params/%s_%s.rds"
  ), .debug[1], .debug[2], .debug[3]) else commandArgs(trailingOnly = TRUE)
  starting_step <- as.integer(tail(.args, 3)[1])
  tariso <- tail(.args, 4)[1]
  params <- readRDS(.args[1])
  urbfrac <- readRDS(.args[2])[iso3 == tariso, value / 100]
  case.dt <- readRDS(.args[3])[iso3 == tariso, .(date, cases)]
  timings <- readRDS(.args[4])
  intros.dt <- readRDS(.args[5])[iso3 == tariso]
  sample <- readRDS(.args[6])
  cm_path = tail(.args, 2)[1]
  outfile <- tail(.args, 1)

  cm_force_rebuild = F;
  cm_build_verbose = F;
  cm_force_shared = T
  cm_version = 2

  load("NGM.rda")
}
fitslc <- seq(starting_step, by = 1, length.out = 20)
bootstrap.dt <- sample[fitslc]
case.dt[, croll := frollmean(cases, align = "center", 7)]

day0 <- as.Date(intros.dt[, min(date)])
intros <- intros.dt[,
  intro.day := as.integer(date - date[1])
][, Reduce(c, mapply(rep, intro.day, infections, SIMPLIFY = FALSE))]

params$date0 <- day0
params$pop[[1]]$seed_times <- intros
params$pop[[1]]$size <- round(params$pop[[1]]$size * urbfrac)
params$pop[[1]]$dist_seed_ages <- c(rep(0, 4), rep(1, 6), rep(0, 6))

tarwindow <- as.Date(c("2020-09-01", "2020-10-01"))
tart <- as.numeric(tarwindow - day0)
case.slc <- case.dt[between(date, tarwindow[1], tarwindow[2]), cases]

startrelax <- as.integer(timings[era == "relaxation", start] - day0)
endrelax <- as.integer(min(timings[era == "relaxation", end], tarwindow[2]) - day0)

startpost <- as.integer(timings[era == "transition", start[1]] - day0)

params$time1 <- endrelax

tms <- day0 + startpost
relaxtms <- day0 + startrelax:endrelax
tier2 <- as.Date("2020-08-15")

peakday <- case.dt[date <= "2020-10-01"][which.max(croll), date]
peakt <- as.numeric(peakday - day0)

# load covidm
# suppressPackageStartupMessages({
source(file.path(cm_path, "R", "covidm.R"))
# })


#' reference for all bootstrap evaluation
scheduler <- function(large, small, symp, k, shft) {
  cons <- list(1 - c(0, small, large, small))
  si <- list(rep(1 - symp, 16))

  relaxfact <- 1 - (1 + exp(-k * as.numeric(relaxtms - tier2 - shft)))^-1
  relaxfact <- (1 - relaxfact[1]) + relaxfact
  relaxcons <- lapply(relaxfact, function(rf) 1 - c(0, small, large, small) * rf)
  relaxsi <- lapply(relaxfact, function(rf) 1 - rep(symp, 16) * rf)

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
      values = c(si, relaxsi),
      times = c(tms, relaxtms)
    )
  )
}

pb = txtProgressBar(min = 1, max = length(fitslc), initial = 1)
#' TODO expand sampling
# brower()
print("outer_scope")
fits.dt <- bootstrap.dt[1,
{
  us <- rep(.SD[, as.numeric(.SD), .SDcols = grep("^u_", names(.SD))], each = 2) * umod
  ys <- rep(.SD[, as.numeric(.SD), .SDcols = grep("^y_", names(.SD))], each = 2)
  pop <- params # copy constructor
  pop$pop[[1]]$y <- ys
  pop$pop[[1]]$u <- pop$pop[[1]]$u * us

  pars_int <- optim_sa(function(ps) {
    lrg <- ps[1]
    sml <- ps[2]
    symp <- ps[3]
    if ((lrg < sml) | (lrg < symp)) NA_real_ else {
                                                                                    #' calculate reduced Rt
      (cm_ngm(
        pop, contact_reductions = c(home = 0, work = sml, school = lrg, other = sml),
        fIs_reductions = symp
      )$R0 / post - 1)^2
    } },
                       start = c(0.8, 0.25, 0.25),
                       lower = c(0.1, 0.01, 0.01),
                       upper = c(0.9, 0.9, 0.9)
  )$par
  setTxtProgressBar(pb, .GRP - 0.5)
  lrg <- pars_int[1]; sml <- pars_int[2]; symp <- pars_int[3]
  pars_relax <- optim_sa(function(ps) {
    k <- ps[1]
    shft <- as.integer(ps[2])
    asc <- ps[3]
    pop$schedule <- scheduler(lrg, sml, symp, k, shft)
    print("here")
    browser()
    save(pop,file="/data/lshtm-batch/debug_est_parameters.Rdata")
    sim <- cm_simulate(
      pop, 1,
      model_seed = 42L
    )$dynamics[
      compartment == "cases",
      .(value = sum(value) * asc), by = t
    ]
    print("done simulatin")
    est <- sim[between(t, tart[1], tart[2]), value]
    casefact <- sum((1 - est / case.slc)^2) / length(est)
    pfact <- (1 - sim[which.max(value), t] / peakt)^2
    casefact + pfact
  },
                         start = c(0.1, 0, 0.10),
                         lower = c(0.001, -14, 0.01),
                         upper = c(0.2, 28, 0.5)
  )$par
  pars <- c(pars_int, pars_relax)
  names(pars) <- c("large", "small", "sympt", "k", "shft", "asc")
  setTxtProgressBar(pb, .GRP)
  as.list(pars)
}, by = sample]

saveRDS(bootstrap.dt[fits.dt, on = .(sample)], outfile)

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
#  geom_line(
#    aes(date, cases),
#    data = case.dt,
#    color = "grey", inherit.aes = FALSE
#  ) +
#   geom_line(
#     aes(date, croll),
#     data = case.dt,
#     color = "black", inherit.aes = FALSE
#   ) +
#  scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
#  scale_y_log10() +
#  theme_minimal()
