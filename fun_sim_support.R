
suppressPackageStartupMessages({
  require(data.table)
})

.args <- if (interactive()) c(
  "sim_support.rda"
) else commandArgs(trailingOnly = TRUE)

#' @param basep a covidm population object
#' @param urbanfaction a numeric vector; either length 1 or length == size of age groups in basep
#' @param day0 a Date scalar
#' 
#' @return a covidm population object, updated with parameters
popsetup <- function(basep, urbanfraction, day0, intro_ages = c(rep(0,4), rep(1, 6), rep(0, 6))) {
  basep$date0 <- day0
  basep$pop[[1]]$size <- round(basep$pop[[1]]$size*urbanfraction)
  basep$pop[[1]]$dist_seed_ages <- intro_ages
  basep
}

scheduler <- function(large, small, symp, k, shft) {
  cons <- list(1-c(0, small, large, small))
  si <- list(rep(1-symp, 16))
  
  relaxfact <- 1-(1+exp(-k*as.numeric(relaxtms-tier2-shft)))^-1
  relaxfact <- (1-relaxfact[1]) + relaxfact
  relaxcons <- lapply(relaxfact, function(rf) 1-c(0, small, large, small)*rf)
  relaxsi <- lapply(relaxfact, function(rf) 1-rep(symp, 16)*rf)
  
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

cm_modifier <- function(parameter, pops, mode, values, times) as.list(environment())

cm_multiplier <- function(parameter, values, times, pops = numeric()) cm_modifier(parameter, pops, mode = "multiply", values, times)

