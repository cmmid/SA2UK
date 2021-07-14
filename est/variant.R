suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("analysis", "GHA")
.args <- if (interactive()) sprintf(c(
  "%s/gen/intervention_timing/%s.rds",
  "%s/gen/pops/%s.rds",
  "%s/gen/yuqs/%s.rds",
  "%s/gen/mobility.rds",
  "%s/est/r0/%s.rds",
  "%s/est/params/%s.rds",
  "%s/sim/history/%s.rds",
  "%s/est/sample/%s.rds",
  .debug[2],
  "covidm",
  "%s/est/variant/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

timing <- readRDS(.args[1])[(period == 1 & era == "transition") | (period == 3 & era == "pre")]
pop <- readRDS(.args[2])

tariso <- tail(.args, 3)[1]

#' targets included here in `variant` column
mob <- readRDS(.args[4])[iso3 == tariso & between(date, timing[period == 3]$start[1], timing[period == 3]$end[1])]
variantRt <- readRDS(.args[5])[period == 3 & era == "pre", .(sample, variant = value)]
fits <- readRDS(.args[6])[variantRt, on=.(sample), nomatch = 0]

#' get susceptible depletion
AR.dt <- readRDS(.args[7])[compartment == "R" & date == max(date), .(sample, group, value = rv) ]
#' get case references for microdistancing
cases.dt <- readRDS(.args[7])[
  compartment == "cases" & (date == timing[period == 3, start] | date == timing[period == 1, start]),
  .(value = sum(value)),
  keyby = .(sample, date)
]

contact_schedule <- mob[,c(
  home = 1, work = prod(work_multiplier)^(1/.N),
  other = prod(other_multiplier)^(1/.N), school = prod(school_multiplier)^(1/.N)
)]

AR.dt[, AR := value / pop$pop[[1]]$size[as.numeric(group)] ]

us <- function(sdt.row, umod = sdt.row[1, umod]) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^u_",names(sdt.row))], each = 2)*umod
}

ys <- function(sdt.row) {
  rep(sdt.row[, as.numeric(.SD), .SDcols = grep("^y_",names(sdt.row))], each = 2)
}

cm_path = tail(.args, 2)[1]
cm_force_rebuild = F;
cm_build_verbose = F;
cm_force_shared = T
cm_version = 2

suppressPackageStartupMessages({
  source(file.path(cm_path, "R", "covidm.R"))
})

fIsfun <- function(case0, k, case) { 1/(1+exp(-k*(case-case0))) }

fIsbaseline <- function(fIs, sympt) { (sympt - fIs)/(1 - fIs) }

cases.dt[fits[, .(sample, case0, k)], fIs := fIsfun(case0, k, value), on=.(sample)]
cases.dt[fits[, .(sample, sympt)], baseline := fIsbaseline(fIs, sympt) ]

sims <- rbindlist(lapply(1:nrow(fits), function(i) {
  AR <- AR.dt[sample == i, AR]
  sdt <- fits[i]
  testpop <- params;
  testpop$pop[[1]]$y <- ys(sdt)
  testpop$pop[[1]]$u <- us(sdt)
  fIsbase <- fIsbaseline(case0, k, sympt, caseref.dt[sample == i & date == min(date), value])
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
