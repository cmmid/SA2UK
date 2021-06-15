suppressPackageStartupMessages({
  require(data.table)
})

#' fixed stride of 20; adjust starting point
.debug <- c("analysis/est", "PAK")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/yuqs/%s.rds",
  "%s/outputs/r0/%s.rds",
  "%s/outputs/sample/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

yuqs <- readRDS(.args[1])
Rts <- dcast(
  readRDS(.args[2])[era %in% c("pre", "post")],
  sample + period ~ era, value.var = "value"
)

yusamp <- yuqs[sample_id <= max(Rts$sample)][,
  .SD, .SDcols = -c("trial","chain", "lp", "ll", "mult", "size")
]

bootstrap.dt <- Rts[yusamp, on=.(sample = sample_id), nomatch = 0]
bootstrap.dt[, umod := pre/baseR ]

saveRDS(bootstrap.dt, tail(.args, 1))
