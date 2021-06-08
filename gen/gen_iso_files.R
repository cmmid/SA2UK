suppressPackageStartupMessages({
  require(countrycode)
  require(data.table)
})

.debug <- "africa"
.args <- if (interactive()) sprintf(c(
  "~/Dropbox/covidLMIC/inputs/isos/all.iso",
  "%s",
  "~/Dropbox/covidLMIC/inputs/isos/%s.iso"
), .debug) else commandArgs(trailingOnly = TRUE)

base.iso <- fread(.args[1], header = FALSE, col.names = "iso3")

fwrite(
  base.iso[,
  .(iso3, continent=countrycode(iso3, "iso3c", "continent"))
  ][tolower(continent) == .args[2], .(iso3)],
  tail(.args, 1), col.names = FALSE
)