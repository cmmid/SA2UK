
suppressPackageStartupMessages({
  require(data.table)
  require(countrycode)
})

.debug <- "~/Dropbox/Covid_LMIC/All_Africa_paper/inputs"
.args <- if (interactive()) sprintf(c(
  "https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv",
  "%s/mobility.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

dt <- fread(.args[1])

national <- dt[
  sub_region_1 == "" & metro_area == "",
  .(
    iso3 = countrycode(country_region_code, "iso2c", "iso3c"),
    date,
    work_multiplier = workplaces_percent_change_from_baseline/100 + 1,
    other_multiplier = (
      grocery_and_pharmacy_percent_change_from_baseline +
      retail_and_recreation_percent_change_from_baseline +
      transit_stations_percent_change_from_baseline  
    )/300 + 1
  )
]

window.width <- 7

national[order(date),
  c("workr", "otherr") := .(
    frollapply(work_multiplier, window.width, function(x) prod(x, na.rm = TRUE)^(1/window.width)),
    frollapply(other_multiplier, window.width, function(x) prod(x, na.rm = TRUE)^(1/window.width))
  ),
  by = iso3
]

#' @example 
#' require(ggplot2)
#' ggplot(national[iso3=="PAK"]) + aes(date) +
#' geom_line(aes(y=otherr, color = "other")) +
#' geom_line(aes(y=workr, color = "work")) +
#' scale_x_date(NULL, date_breaks = "month", date_label = "%b %Y") +
#' theme_minimal()