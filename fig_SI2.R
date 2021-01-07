suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(patchwork)
})

.args <- if (interactive()) sprintf(c(
  "%s/eras.rds",
  "%s/cfr.rds",
  "%s/SI2.png"
), "~/Dropbox/SA2UK/outputs/figs") else commandArgs(trailingOnly = TRUE)

interventions <- readRDS(.args[1]) + 
  coord_cartesian(
    ylim = c(1, NA),
    xlim = as.Date(c("2020-04-01", "2021-01-01")), expand = FALSE
  ) + scale_fill_manual(guide = "none", values = rep(NA, 4))
cfr <- readRDS(.args[2]) + scale_color_manual(
  NULL,
  breaks = c("dCFR","cCFR","nCFR"),
  labels = c(nCFR="naive", dCFR = "delayed", cCFR = "corrected"),
  values = c(nCFR="firebrick", dCFR = "green", cCFR="darkorchid4"),
  aesthetics = c("color", "fill")
)

res <- interventions / cfr + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(text = element_text(size = 7), panel.grid.minor = element_blank())

ggsave(tail(.args, 1), res, height = 6, width = 7, units = "in", dpi = 300)
