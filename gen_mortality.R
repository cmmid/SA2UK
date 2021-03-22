suppressPackageStartupMessages({
  require(data.table)
  require(countrycode)
  require(wpp2019)
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

get_sex_data <- function(sx) {
  mort <- as.data.table(get(data(list=list(sprintf("mx%s", sx)))[[1]]))[,.(
    iso3c = countrycode(country_code, "iso3n", "iso3c"),
    sex = sx,
    mortality = `2015-2020`, age
  )][!is.na(iso3c)]
  
  pop <- as.data.table(get(data(list=list(sprintf("pop%s", sx)))[[1]]))[,.(
    iso3c = countrycode(country_code, "iso3n", "iso3c"),
    sex = sx,
    pop = `2020`*1000, age
  )][!is.na(iso3c)]
  
  mort[pop, on=.(iso3c, age, sex)]
}

as.data.table(get(data(popF)))[,.(sex=)]
