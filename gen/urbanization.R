suppressPackageStartupMessages({
  require(data.table)
})

.debug <- "~/Dropbox/Covid_LMIC/All_Africa_paper"
.args <- if (interactive()) sprintf(c(
  "%s/inputs/urbanization.rds"
), .debug[1]) else commandArgs(trailingOnly = TRUE)

wburl <- "http://api.worldbank.org/v2/en/indicator/SP.URB.TOTL.IN.ZS?downloadformat=csv"
wbfile <- "API_SP.URB.TOTL.IN.ZS_DS2_*.csv"

tmp <- tempfile()
download.file(wburl, tmp)
in.dt <- fread(cmd=sprintf("unzip -p %s %s", tmp, wbfile), header = T, skip = 4)[,
  .SD, .SDcols = -c("Country Name", "Indicator Name", "Indicator Code", "V66")
]
unlink(tmp)

in.mlt <- melt(in.dt, id.vars = "Country Code", variable.name = "Year")

res <- setkey(
  in.mlt[order(Year, decreasing = TRUE),
    .SD[which.max(!is.na(value))],
    by=.(iso3=`Country Code`)
  ][
    !is.na(value)
  ][, .(
    iso3, value = value/100, year = Year
  )],
  iso3
)

saveRDS(res, tail(.args, 1))