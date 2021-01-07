suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/scenarios/%s.rds",
  "%s/inputs/pops/%s.rds",
  "%s/inputs/yuqs/%s.rds",
  "%s/outputs/r0/%s.rds",
  "%s/outputs/introductions/%s.rds",
  "%s/inputs/urbanization.rds",
  "%s/outputs/projections/%s.rds",
  .debug[2],
  "../covidm",
  "%s/variant/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

#' targets
Rts <- readRDS(.args[4])

#' get susceptible depletion
proj.dt <- readRDS(.args[7])[compartment == "R" & date == "2020-10-15"]
params <- readRDS(.args[2])
urbfrac <- readRDS(.args[6])[iso3 == tariso, value / 100]
params$pop[[1]]$size <- round(params$pop[[1]]$size*urbfrac)
proj.dt[, AR := value / params$pop[[1]]$size[as.numeric(group)] ]
proj.dt[q == "md", q := "med" ]

#' get inventions at post, then relax
day0 <- readRDS(.args[5])[, min(date)]
scens <- readRDS(.args[1])[scen_id == 2]
relaxend <- as.Date("2020-10-15")
relaxstart <- as.Date("2020-05-01")
relaxrate <- c(0.03, 0.035, 0.075)
relaxinit <- 1-(1+exp(-relaxrate*as.numeric(relaxstart-tier2)))^-1
relaxfact <- 1-(1+exp(-relaxrate*as.numeric(relaxend-tier2)))^-1
relaxfact <- 0.5*(relaxfact-0.5)/(relaxinit-0.5)+0.5

scens[, schoolr := school*relaxfact ][, workr := work*relaxfact ][, otherr := other*relaxfact ][, self_isor := self_iso*relaxfact ]

load("NGM.rda")

#' establish baseline pop
yuref <- readRDS(.args[3])[order(eqs)]
qs.inds <- with(yuref, c(which.max(eqs >= 0.25),which.max(eqs >= 0.5),which.max(eqs >= 0.75)))
yuuse <- yuref[qs.inds][, variable := c("lo","med","hi") ][,-c("trial","chain","lp","ll","mult","size")]

run_options <- melt(
  Rts[era == "pre"],
  id.vars = "era",
  measure.vars = c("lo", "med", "hi"),
  value.name = "r0"
)[yuuse, on=.(variable) ][,
  umul := r0 / baseR
][
  melt(
    Rts[era == "variant"],
    id.vars = "era",
    measure.vars = c("lo", "med", "hi"),
    value.name = "varr0"
  ), on=.(variable)
]

withdepletion <- list()
nodepletion <- list()

for (qtar in run_options$variable) {
  opts <- run_options[variable == qtar]
  ar <- proj.dt[q == qtar]$AR
  uf <- opts[, umul*rep(as.numeric(.SD), each = 2), .SDcols = grep("^u_", names(run_options))]
  ys <- opts[, rep(as.numeric(.SD), each = 2), .SDcols = grep("^y_", names(run_options))]
  scn <- scens[q==qtar, .(home = home, school = schoolr, work = workr, other = otherr, self_iso = self_isor)]
  
  withdepletion[[qtar]] <- with(scn,
    opts$varr0/cm_ngm(
      params,
      uf*(1-ar),
      contact_reductions = c(home,work,school,other),
      fIs_reductions = self_iso, ymod = ys
    )$R0
  )
  
  nodepletion[[qtar]] <- with(
    scn,
    opts$varr0/cm_ngm(
      params,
      uf*(1-ar*.75),
      contact_reductions = c(home,work,school,other),
      fIs_reductions = self_iso,
      ymod = ys
    )$R0
  )

}
  
#' now we have the fully depleted population
#' + observed Reff
#' 
#' option 1: full cross protective - what multiple of uf would be required to observe
#' this Reff, given depletion & modification era reductions?
#' that multiple => increase of R0

ufs <- mapply(
  function(br, sfrac) br/refR0 * sfrac,
  br = baseRs,
  sfrac = lapply(1:5, function(ri) 1 - allbind[r_id == ri, AR]),
  SIMPLIFY = FALSE
)

iv_data <- scenario[scen_id != 1 & !is.na(self_iso)][order(start_day)][.N]

varR <- melt(Rts[era == "variant"], id.vars = "era", measure.vars = 2:6)$value

Rmultipliers_with_depletion <- mapply(function(uf, R) with(
  iv_data,
  R/cm_ngm(params_back, R0_multiplier = uf, contact_reductions = c(home,work,school,other), fIs_reductions = self_iso)$R0
), uf = ufs, R = varR)

ufs_non_depl <- mapply(
  function(br, sfrac) br/refR0 * sfrac,
  br = baseRs,
  sfrac = 1,
  SIMPLIFY = FALSE
)

Rmultipliers_non_depletion <- mapply(function(uf, R) with(
  iv_data,
  R/cm_ngm(params_back, R0_multiplier = uf, contact_reductions = c(home,work,school,other), fIs_reductions = self_iso)$R0
), uf = ufs_non_depl, R = varR)

res <- data.table(
  model = c(rep("cross-protected", 5), rep("susceptible", 5)),
  Rfactor = c(Rmultipliers_with_depletion, Rmultipliers_non_depletion)
)

saveRDS(res, tail(.args, 1))