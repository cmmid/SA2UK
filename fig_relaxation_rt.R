suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/relaxation/%s.rds",
  "covidm",
  .debug[2],
  "thing.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

rt.dt <- readRDS(.args[1])

relax.p <- force(ggplot(rt.dt) + aes(date) +
  geom_ribbon(aes(ymin=lower_95, ymax=upper_95), alpha = 0.1) +
  geom_ribbon(aes(ymin=lower_50, ymax=upper_50), alpha = 0.2) +
  geom_line(aes(y=median)) +
  coord_cartesian(ylim = c(0.5, 1.5), expand = FALSE) +
  scale_y_continuous(expression(R[t])) +
  scale_x_date(name=NULL, date_breaks = "months", date_minor_breaks = "weeks", date_labels = "%b") +
  theme_minimal())

saveRDS(relax.p, tail(.args, 1))