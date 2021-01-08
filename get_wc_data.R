suppressPackageStartupMessages({
  require(data.table)
})

.debug <- "~/Dropbox/SA2UK"
.args <- if (interactive()) sprintf(c(
  "%s/inputs/epi_data.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

cases.dt <- readr::read_csv("https://raw.githubusercontent.com/dsfsi/covid19za/master/data/covid19za_provincial_cumulative_timeline_confirmed.csv") %>% 
    dplyr::select(date, total_cases = WC) %>%
    data.table::data.table()

deaths.dt <- readr::read_csv("https://raw.githubusercontent.com/dsfsi/covid19za/master/data/covid19za_provincial_cumulative_timeline_deaths.csv") %>% 
    dplyr::select(date, total_deaths = WC) %>%
    data.table::data.table()

cases.dt <- cases.dt[, "cases" := c(NA, diff(total_cases))][, c("date", "cases")]
deaths.dt <- deaths.dt[, "deaths" := c(NA, diff(total_deaths))][,c("date", "deaths")]

final.dt <- merge(cases.dt, deaths.dt)

saveRDS(final.dt, here::here("epi_wc_data.RDS"))
