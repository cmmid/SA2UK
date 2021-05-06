suppressPackageStartupMessages({
  require(data.table)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper/outputs", "PAK")
.args <- if (interactive()) sprintf(c(
  "%s/params", .debug[2], "%s/params/%s_consolidated.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

fls <- list.files(.args[1], sprintf("%s_\\d+\\.rds", .args[2]), full.names = TRUE)

saveRDS(rbindlist(lapply(fls, readRDS)), tail(.args, 1))
