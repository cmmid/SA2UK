suppressPackageStartupMessages({
  require(data.table)
  require(countrycode)
})

.args <- if (interactive()) sprintf(c(
  "africaisos.txt"
)) else commandArgs(trailingOnly = TRUE)

#' is.na filter excludes Somaliland
res <- as.data.table(countrycode::codelist)[continent == "Africa", .(iso3 = iso3c)][!is.na(iso3)]

fwrite(
  res,
  tail(.args, 1),
  col.names = FALSE
)