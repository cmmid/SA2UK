
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c(".", "NGA_0")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/scenario/%s.rds",
  .debug[2], # PAK
  "%s/outputs/acc_scen/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

refdate <- as.Date("2021-09-01")

raw <- readRDS(.args[1])

raw[, iyear := ceiling(as.numeric(date - refdate)/365) ]

res <- raw[,.(value = sum(rv)), keyby=.(epi_id, sample, iyear, compartment, group)]

saveRDS(res, tail(.args, 1))