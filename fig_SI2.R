suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(patchwork)
})

.args <- if (interactive()) sprintf(c(
  "%s/eras.rds",
  "~/Dropbox/SA2UK/outputs/cfrs.rds",
  "%s/SI2.png"
), "~/Dropbox/SA2UK/outputs/figs") else commandArgs(trailingOnly = TRUE)

interventions <- readRDS(.args[1]) + 
  coord_cartesian(
    ylim = c(1, NA),
    xlim = as.Date(c("2020-04-01", "2021-01-07")), expand = FALSE
  ) + scale_fill_manual(guide = "none", values = rep(NA, 4))

plot.dt <- readRDS(.args[2])[province == "all" & date > "2020-05-15"]

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
                   labels = c(nCFR="naive", dCFR = "lagged", cCFR = "delay-adjusted"),
                   values = c(nCFR="firebrick", dCFR = "green", cCFR="darkorchid4"),
                   aesthetics = c("color", "fill")
                 ) + theme_minimal())


res <- interventions / cfr.p + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(text = element_text(size = 7), panel.grid.minor = element_blank())

ggsave(tail(.args, 1), res, height = 6, width = 7, units = "in", dpi = 300)
