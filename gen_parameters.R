suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/SA2UK/outputs/params", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s",
  "%s/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

fns <- list.files(.args[1], "\\d+\\.rds", full.names = TRUE)

res <- rbindlist(lapply(fns, readRDS))

saveRDS(res, tail(.args, 1))
