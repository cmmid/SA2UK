suppressPackageStartupMessages({
  require(data.table)
  require(countrycode)
  require(wpp2019)
  require(socialmixr)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/mortality.rds"
), .debug[1]) else commandArgs(trailingOnly = TRUE)

get_sex_data <- function(sx, covidm_age_ulim = 75, covid_m_age_w = 5) {
  age_low_lims <- seq(0, covidm_age_ulim, by=covid_m_age_w)
  
  # need pop to weight 75+ mortality
  pop <- as.data.table(get(data(list=list(sprintf("pop%s", sx)))[[1]]))[,.(
    iso3c = countrycode(country_code, "iso3n", "iso3c"),
    pop = `2020`*1000,
    age = as.integer(gsub("^(\\d+)[^\\d]+$","\\1",age))
  )][!is.na(iso3c)]
  
  # get unweighted mortality
  # going to use [1,5) mortality for infants as well,
  # and deduct "excess" infant mortality from birthrate
  mort <- as.data.table(get(data(list=list(sprintf("mx%s", sx)))[[1]]))[,.(
    iso3c = countrycode(country_code, "iso3n", "iso3c"),
    mortality = `2015-2020`, age
  )][!is.na(iso3c)][order(age),
    .(
      age = age[-2],
      mortality = mortality[-1]
    ),
    by=iso3c
  ]
  
  mort[pop, on=.(iso3c, age)][, sex := sx ]
}

f <- get_sex_data("F")
m <- get_sex_data("M")
