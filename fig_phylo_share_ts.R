suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- "~/Dropbox/SA2UK"
.args <- if (interactive()) sprintf(c(
  "%s/outputs/phylo.rds",
  "%s/outputs/figs/phylo.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

plot.dt <- readRDS(.args[1])

var.label <- "501Y.V2"

p.phylo <- force(ggplot(plot.dt[!is.na(rolling.var)]) + aes(date) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.1, fill = "red") +
  geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.2, fill = "red") +
  geom_line(aes(y=binop), color = "red") +
  scale_x_date(
    name = NULL,
    date_breaks = "months", date_minor_breaks = "weeks",
    date_labels = "%b"
  ) +
  scale_y_continuous(sprintf("%s Fraction", var.label)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(as.Date("2020-10-01"), as.Date("2021-01-01")), expand = FALSE) +
  theme_minimal())

saveRDS(p.phylo, tail(.args, 1))