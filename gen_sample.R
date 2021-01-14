suppressPackageStartupMessages({
  require(data.table)
  require(optimization)
})

#' fixed stride of 20; adjust starting point
.debug <- c("~/Dropbox/SA2UK","ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/yuqs/%s.rds",
  "%s/outputs/r0/%s.rds",
  "%s/outputs/sample/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

yuqs <- readRDS(.args[1])
Rts <- readRDS(.args[2])

yusamp <- yuqs[
  sample(.N, max(Rts$sample))
][,
  .SD, .SDcols = -c("trial","chain","eqs","lp", "ll", "mult", "size")
][, sample := 1L:.N ]

bootstrap.dt <- Rts[yusamp, on=.(sample)]
bootstrap.dt[, umod := pre/baseR ]

saveRDS(bootstrap.dt[fits.dt, on = .(sample)], tail(.args, 1))
