suppressPackageStartupMessages({
  require(data.table)
  require(wpp2019)
  require(countrycode)
})

if (sys.nframe() == 0) {
  .debug <- "~/Dropbox/covidLMIC"
  .args <- if (interactive()) sprintf(c(
    "%s/inputs/populations.rds"
  ), .debug) else commandArgs(trailingOnly = TRUE)
  outfile <- tail(.args, 1)
}
data(popF)
data(popM)

#' MAGIC NUMBER WARNING
#' other elements of analysis use 16 age categories, w/ 5 year intervals
#' however, the IFR data is 10 age categories
#' here, we just get into manageable format relative to wpp2019 raw
#' wpp is denominated in 1000s of capita

reformat <- function(dt, gen) as.data.table(dt)[, .(
  age = factor(age, levels = unique(age), ordered = TRUE),
  pop = `2020`*1000
), keyby=.(iso3 = countrycode(country_code, "iso3n", "iso3c"))][
  !is.na(iso3)
][, gender := gen ]

pop.dt <- rbind(
  reformat(popF, "F"),
  reformat(popM, "M")
)[, .(pop = sum(pop)), keyby=.(iso3, age)]

saveRDS(pop.dt, outfile)