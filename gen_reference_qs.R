
suppressPackageStartupMessages({
  require(data.table)
  require(qs)
})
if (sys.nframe() == 0) {
  .debug <- c("~/Dropbox/SA2UK", "ZAF")
  .args <- if (interactive()) sprintf(c(
    "%s/inputs/covidm_fit_yu.qs",
    "%s/inputs/pops/%s.rds",
    "%s/inputs/yuqs/%s.rds"
  ), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)
  load("NGM.rda")
}

yu <- qread(.args[1])
pop <- readRDS(.args[2])


uids <- rep(grep("^u_", colnames(yu)), each = 2)
yids <- rep(grep("^y_", colnames(yu)), each = 2)

qs.dt <- yu[, {
  umod <- as.numeric(.SD[1,])[1:16]
  ymod <- as.numeric(.SD[1,])[17:32]
  ngm <- cm_ngm(
    pop,
    R0_multiplier = umod, ymod = ymod
  )
  .(baseR = ngm$R0, si = cm_generation_time(pop, ymod = ymod, ngm = ngm))
}, by=.(trial, chain), .SDcols = c(uids, yids)]

ret <- qs.dt[order(baseR), eqs := (1:.N)/.N ][yu, on=.(trial, chain), nomatch = 0]

saveRDS(ret, tail(.args, 1))
