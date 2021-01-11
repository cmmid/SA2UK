suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/params/%s.rds",
  "%s/inputs/pops/%s.rds",
  "%s/inputs/urbanization.rds",
  "%s/outputs/projections/%s.rds",
  .debug[2],
  "%s/variant/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

#' targets included here in `variant` column
fits <- readRDS(.args[1])

#' get susceptible depletion
proj.dt <- readRDS(.args[4])[compartment == "R" & date == max(date) ]
cmpdate <- proj.dt[1, date]

params <- readRDS(.args[2])
urbfrac <- readRDS(.args[3])[iso3 == tariso, value / 100]
params$pop[[1]]$size <- round(params$pop[[1]]$size*urbfrac)
proj.dt[, AR := value / params$pop[[1]]$size[as.numeric(group)] ]

load("NGM.rda")

tier2 <- as.Date("2020-08-15")

rfs <- function(k, shft) {
  relaxref <- (1+exp(-k*as.numeric(as.Date("2020-05-01")-tier2-shft)))^-1
  relaxref + 1 - (1+exp(-k*as.numeric(cmpdate-tier2-shft)))^-1
}

fits[, relaxfactor := rfs(k, shft) ]

sims <- rbindlist(lapply(1:nrow(fits), function(i) with(as.list(fits[i,.(variantR=variant, large, small, sympt, relaxfactor)]), {
  AR <- proj.dt[sample == i, AR]
  us <- rep(fits[i, as.numeric(.SD)*umod, .SDcols = grep("^u_",names(fits))], each = 2)
  ys <- rep(fits[i, as.numeric(.SD), .SDcols = grep("^y_",names(fits))], each = 2)
  testpop <- params; testpop$pop[[1]]$y <- ys
  testpop$pop[[1]]$u <- testpop$pop[[1]]$u*us
  withdepl <- variantR/cm_ngm(
    testpop,
    1-AR,
    contact_reductions = c(0, small, large, small)*relaxfactor,
    fIs_reductions = sympt*relaxfactor
  )$R0
  fixeddepl <- variantR/cm_ngm(
    testpop,
    1-AR*.75,
    contact_reductions = c(0, small, large, small)*relaxfactor,
    fIs_reductions = sympt*relaxfactor
  )$R0
  eff <- optimize(function(eff) (variantR/cm_ngm(
    testpop,
    1-AR*eff,
    contact_reductions = c(0, small, large, small)*relaxfactor,
    fIs_reductions = sympt*relaxfactor
  )$R0 - 1)^2, interval = c(0.01,0.99))$minimum
  
  list(sample=i, withdepl = withdepl, fixeddepl = fixeddepl, tarescape = 1-eff)
})))

saveRDS(sims, tail(.args, 1))