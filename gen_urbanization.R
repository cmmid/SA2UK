suppressPackageStartupMessages({
  require(data.table)
})

.args <- if (interactive()) c(
  "tmp.csv", "urbanization.rds"
) else commandArgs(trailingOnly = TRUE)

in.dt <- fread(.args[1], header = T, skip = 4)[,
  .SD, .SDcols = -c("Country Name", "Indicator Name", "Indicator Code", "V66")
]

in.mlt <- melt(in.dt, id.vars = "Country Code", variable.name = "Year")

res <- in.mlt[
  order(Year, decreasing = TRUE),
  .SD[which.max(!is.na(value))],
  by=.(iso3=`Country Code`)
][!is.na(value)]

saveRDS(res, tail(.args, 1))