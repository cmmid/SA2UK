suppressPackageStartupMessages({
  require(data.table)
  require(countrycode)
})

.debug <- c("analysis")
.args <- if (interactive()) sprintf(c(
  "%s/ins/lifetables.csv",
  "%s/gen/mortality.rds"
), .debug[1]) else commandArgs(trailingOnly = TRUE)

#' covidm population not differentiated by sex, so want total value (SexID == 3)
#' exclude aggregate populations (LocID >= 900) and Channel Islands (830)
wpp.dt <- fread(.args[1])[SexID == 3][!(LocID >= 900) & LocID != 830]

#' TODO: is this the correct interpretation of WPP life tables?
#' an alternative reading of data description suggests that different ages values
#' should come from different time periods (0-4 from 2015-2020, 5-9 from 2010-2015, etc)
#' goal is to have best approximation of expected non-COVID mortality, today, by age
#' since model doesn't have a 0-1 age group, combine that with 1-4 category
#' model stops at 75+, so use expected additional life years at 75 for death rate
mort.dt <- wpp.dt[Time == "2015-2020" & AgeGrpStart <= 75, .(
  AgeGrp = AgeGrp[-1], AgeGrpStart = AgeGrpStart[-2], AgeGrpSpan = 5,
  px = c(px[1]*px[2], px[-c(1:2)]),
  ex = ex[-1]
), by=.(iso3 = countrycode(LocID, "iso3n", "iso3c"))]

#' now need to convert this into ODE compartment exit rate for this age strata
#' px = probability of not dying in this age strata
#' the aging rate is 1/(5*365.25) (for daily time step)
#' px = aging rate / (aging rate + death rate)
#' death rate = aging*(1 / px - 1)  = 

mort.dt[, per_capita_day := (1/(5*365.25))*(1/px - 1) ]
mort.dt[AgeGrpStart == 75, per_capita_day := 1/(ex*365.25) ]
mort.dt[, per_capita_day_alt := -log(px)/AgeGrpSpan/365.25 ]

saveRDS(mort.dt, tail(.args, 1))