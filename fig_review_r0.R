
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/SA2UK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/r0", "summary_r0.png"
), .debug) else commandArgs(trailingOnly = TRUE)

r0s <- list.files(.args[1])

extractR0 <- function(fn, pth = .args[1]) {
  readRDS(file.path(pth, fn))[, iso3 := gsub("^.*(.{3})\\.rds$", "\\1", fn)]
}

allR0s <- rbindlist(lapply(r0s, extractR0))
R0.dt <- melt.data.table(allR0s, id.vars = c("iso3","sample"), variable.name = "era")

p <- ggplot(R0.dt) + aes(era, value) +
  geom_point() +
  theme_minimal()
