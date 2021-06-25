suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper", "PAK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/params/%s_consolidated.rds",
  "%s/inputs/pops/%s.rds",
  "%s/outputs/projections/%s.rds",
  "%s/outputs/sample/%s.rds",
  "%s/inputs/mobility.rds",
  "%s/outputs/adj_data.rds",
  .debug[2],
  "covidm",
  "%s/outputs/variant/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 3)[1]

#' targets included here in `variant` column
ref <- readRDS(.args[4])[period == 3][, .(variant = pre), by=sample]
fits <- readRDS(.args[1])[ref, on=.(sample), nomatch = 0]

#' get susceptible depletion
proj.dt <- readRDS(.args[3])[compartment == "R" & date == max(date) ]
mob <- readRDS(.args[5])[iso3 == tariso & between(date, proj.dt$date[1], proj.dt$date[1]+7)]

casef <- readRDS(.args[6])[
  (iso3 == tariso),
  .(date, croll = frollmean(cases, 7, align = "center"))
][date == c(proj.dt$date[1]), croll]

contact_schedule <- mob[,c(
  home = 1, work = prod(work_multiplier)^(1/.N),
  other = prod(other_multiplier)^(1/.N), school = prod(school_multiplier)^(1/.N)
)]

cmpdate <- proj.dt[1, date]

params <- readRDS(.args[2])
proj.dt[, AR := rv / params$pop[[1]]$size[as.numeric(group)] ]

us <- function(sdt.row, umod = sdt.row[1, umod]) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^u_",names(sdt.row))], each = 2)*umod
}

ys <- function(sdt.row) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^y_",names(sdt.row))], each = 2)
}

fIs_baseline <- function(
  case0, k, reff, refcase = case.slc[1]
) (reff*(1+exp(-k*(refcase-case0))) - 1)/exp(-k*(refcase-case0))

fIs_amp <- function(
  case0, k, # fit elements
  reff, # sampling element: value of function at css[1]
  css = case.slc, # data element
  baseline = fIs_baseline(case0, k, reff, css[1]) # entailed remaining coefficient
) (1-baseline)/(1+exp(-k*(css-case0))) + baseline

cm_path = tail(.args, 2)[1]
cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T
cm_version = 2

suppressPackageStartupMessages({
  source(file.path(cm_path, "R", "covidm.R"))
})

sims <- rbindlist(lapply(1:nrow(fits), function(i) {
  AR <- proj.dt[sample == i, AR]
  sdt <- fits[i]
  testpop <- params;
  testpop$pop[[1]]$y <- ys(sdt)
  testpop$pop[[1]]$u <- us(sdt)
  fIsred <- fIs_amp(sdt$case0, sdt$k, sdt$sympt, casef, baseline = sdt$fIsbaseline)
  withdepl <- sdt$variant/cm_eigen_ngm(
    testpop, contact_reductions = 1-contact_schedule,
    u_multiplier = 1-AR,
    fIs_reductions = fIsred
  )$R0
  fixeddepl <- sdt$variant/cm_eigen_ngm(
    testpop, contact_reductions = 1-contact_schedule,
    u_multiplier = 1-AR*.75,
    fIs_reductions = fIsred
  )$R0
  eff <- optimize(function(eff) (sdt$variant/cm_eigen_ngm(
    testpop,
    u_multiplier = 1-AR*eff,
    contact_reductions = 1-contact_schedule,
    fIs_reductions = fIsred
  )$R0 - 1)^2, interval = c(0.01,0.99))$minimum
  
  list(sample=i, withdepl = withdepl, fixeddepl = fixeddepl, tarescape = 1-eff)
}))

saveRDS(sims, tail(.args, 1))
