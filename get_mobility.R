
suppressPackageStartupMessages({
  require(data.table)
  require(countrycode)
})

.oxpat <- "https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/timeseries/c1_%s.csv"
.debug <- "~/Dropbox/Covid_LMIC/All_Africa_paper"
.args <- if (interactive()) sprintf(c(
  "https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv",
  sprintf(.oxpat, c("school_closing", "flag")),
  "%s/inputs/mobility.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

dt <- fread(.args[1])

readox <- function(url) setkey(melt(
  fread(url, drop = c(1,3)),
  id.vars = "country_code"
)[, .(
  iso3 = country_code,
  date = as.Date(variable, "%d%b%Y"),
  value
)], iso3, date)

sch.dt <- readox(.args[2])
#' the warning here is because some of the early all-NA columns are
#' introspected as logical, and then subsequently changed to integers
#' this is fine; the alternative explicitly encoding that is substantial work
schflag.dt <- readox(.args[3])

sch.dt[schflag.dt, flag := i.value ]

sch.dt[, school_status := fifelse(is.na(flag), value / 3, value / 3 * flag) ]

national <- setkey(dt[
  sub_region_1 == "" & metro_area == "",
  .(
    iso3 = fifelse(is.na(country_region_code), "NAM", countrycode(country_region_code, "iso2c", "iso3c")),
    date,
    work_multiplier = workplaces_percent_change_from_baseline/100 + 1,
    other_multiplier = (
      grocery_and_pharmacy_percent_change_from_baseline +
      retail_and_recreation_percent_change_from_baseline +
      transit_stations_percent_change_from_baseline  
    )/300 + 1
  )
], iso3, date)

national[sch.dt, school_multiplier := 1 - school_status ]

window.width <- 7

national[order(date),
  c("workr", "otherr") := .(
    frollapply(work_multiplier, window.width, function(x) prod(x, na.rm = TRUE)^(1/window.width)),
    frollapply(other_multiplier, window.width, function(x) prod(x, na.rm = TRUE)^(1/window.width))
  ),
  by = iso3
]

#' going to use these to adjust contact matrices

saveRDS(national, tail(.args, 1))

#' @example 
#' require(ggplot2)
#' ggplot(national[iso3=="PAK"]) + aes(date) +
#' geom_line(aes(y=otherr, color = "other")) +
#' geom_line(aes(y=workr, color = "work")) +
#' geom_line(aes(y=school_multiplier, color = "school")) +
#' scale_x_date(NULL, date_breaks = "month", date_label = "%b %Y") +
#' theme_minimal()