
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/SA2UK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/cfrs.rds",
  "%s/outputs/figs/ccfr.rds"
), .debug) else commandArgs(trailingOnly = TRUE)

plot.dt <- readRDS(.args[1])[province == "WC" & date > "2020-05-15"]

cfr.p <- force(ggplot(plot.dt) + aes(date, md) +
  geom_line(aes(color = ver)) +
  geom_ribbon(aes(fill = ver, ymin = lo.lo, ymax = hi.hi), alpha = 0.1) +
  geom_ribbon(aes(fill = ver, ymin = lo, ymax = hi), alpha = 0.25) +
  
  coord_cartesian(
    ylim = c(0, .075),
    xlim = as.Date(c("2020-04-01", NA)),
    expand = FALSE
  ) +
  scale_y_continuous(
    "Case Fatality Ratio (CFR)",
    breaks = c(0,0.025,0.05,0.075),
    labels = function(v) sprintf("%0.2g%%", v*100)
  ) +
  scale_x_date(name = NULL, date_breaks = "months", date_minor_breaks = "weeks", date_labels = "%b") +
  scale_color_manual(
    NULL,
    breaks = c("dCFR","cCFR","nCFR"),
    labels = c(nCFR="naive", dCFR = "delayed", cCFR = "corrected"),
    values = c(nCFR="orchid", dCFR = "mediumorchid", cCFR="darkorchid4"),
    aesthetics = c("color", "fill")
  ) + theme_minimal())

saveRDS(cfr.p, tail(.args, 1))
