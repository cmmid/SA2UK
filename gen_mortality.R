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

data(mxM)
m_mortality <- as.data.table(mxM)[,.(
  iso3c = countrycode(country_code, "iso3n", "iso3c"),
  male = `2015-2020`, age
)][!is.na(iso3c)]

data(mxF)
f_mortality <- as.data.table(mxF)[,.(
  iso3c = countrycode(country_code, "iso3n", "iso3c"),
  female = `2015-2020`, age
)][!is.na(iso3c)]

get_sex_data <- function(sx, covidm_age_ulim = 75, covid_m_age_w = 5) {
  age_low_lims <- seq(0, covidm_age_ulim, by=covid_m_age_w)
  pop <- as.data.table(get(data(list=list(sprintf("pop%s", sx)))[[1]]))[,.(
    iso3c = countrycode(country_code, "iso3n", "iso3c"),
    pop = `2020`*1000,
    age = as.integer(gsub("^(\\d+)[^\\d]+$","\\1",age))
  )][
    !is.na(iso3c),
    pop_age(.SD, age_low_lims, "age", "pop"),
    by=iso3c
  ]
  
  mort <- as.data.table(get(data(list=list(sprintf("mx%s", sx)))[[1]]))[,.(
    iso3c = countrycode(country_code, "iso3n", "iso3c"),
    mortality = `2015-2020`, age
  )][!is.na(iso3c)][order(age),{
      # merging the first two age categories
      # don't have explicit <1 and [1,5) populations
      # so assume the mortality in under 0 determines relative population
      # in [1,5)
    .(
      sex = sx, age = age[-2]
    )},
    by=iso3c
  ]
  
  
  
  mort[pop, on=.(iso3c, age, sex)]
}

as.data.table(get(data(popF)))[,.(sex=)]
