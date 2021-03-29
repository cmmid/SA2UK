suppressPackageStartupMessages({
  require(data.table)
  require(countrycode)
})

.debug <- "~/Dropbox/Covid_LMIC/All_Africa_paper"
.args <- if (interactive()) sprintf(c(
  "%s/inputs/epi_data.rds",
  "africaisos.txt"
), .debug[1]) else commandArgs(trailingOnly = TRUE)

epiisos <- readRDS(.args[1])[continent == "Africa", unique(iso3)]

#' is.na filter excludes Somaliland
res <- as.data.table(countrycode::codelist)[
  continent == "Africa", .(iso3 = iso3c)
][
  !is.na(iso3)
]

warning(sprintf(
  "Excluding ISOs based on availability of epi data; none for: %s.",
  res[!(iso3 %in% epiisos), paste(iso3, collapse=", ")] 
))

res <- res[iso3 %in% epiisos]

fwrite(
  res,
  tail(.args, 1),
  col.names = FALSE
)