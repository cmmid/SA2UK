
suppressPackageStartupMessages({
  require(data.table)
})

case.dt <- readRDS("~/Dropbox/SA2UK/inputs/prov_data.rds")[,
  .(cases = sum(cases)), by=date
]

testing.dt <- fread(
  "https://raw.githubusercontent.com/dsfsi/covid19za/master/data/covid19za_timeline_testing.csv"
)[, .(date = as.Date(date, "%d-%m-%Y"), tests = c(cumulative_tests[1], diff(cumulative_tests)))]

comb.dt <- case.dt[testing.dt, on=.(date)]
comb.dt[is.na(cases), cases := 0 ]

comb.dt[, positivity := cases / tests ]

ggplot(comb.dt) + aes(date, positivity) +
  geom_line() + theme_minimal() +
  scale_x_date(NULL, date_breaks = "months", date_labels = "%b") + scale_y_log10()
