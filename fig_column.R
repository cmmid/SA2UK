
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(patchwork)
})

.args <- if (interactive()) sprintf(c(
  "%s/ccfr.rds", "%s/timeseries.rds", "%s/AR.rds",
  "%s/combined.png"
), "~/Dropbox/SA2UK/outputs/figs") else commandArgs(trailingOnly = TRUE)

zlow <- as.Date("2020-10-15")
zend <- as.Date("2021-01-01")

altcfrscale <-   list(scale_color_manual(
  NULL,
  breaks = "cCFR",
  values = c(cCFR = "darkorchid4", nCFR=NA, dCFR=NA),
  aesthetics = c("color", "fill"), guide = "none"
), scale_y_continuous(
  "Corrected Case Fatality Ratio",
  breaks = c(0,0.025,0.05,0.075),
  labels = function(v) sprintf("%0.2g%%", v*100)
))

phylo <- readRDS(.args[1]) + coord_cartesian(
  xlim = c(zlow, zend),
  ylim = c(0, 1),
  expand = FALSE
)
cfr <- readRDS(.args[1]) + theme(
  legend.position = c(0.4, 1), legend.justification = c(0.4, 1),
  legend.key.height = unit(.75, "line")
) + altcfrscale

ts <- readRDS(.args[2]) + theme(
  legend.direction = "horizontal",
  legend.position = c(0.5, 0), legend.justification = c(0.5, 0),
  legend.margin = margin(t=-1, b=-.25, unit = "line")
)
enddate <- "2021-01-07"

ar <- readRDS(.args[3]) + coord_cartesian(
  xlim = as.Date(c("2020-04-01",enddate)),
  ylim = c(0, 1),
  expand = FALSE
) + theme(
  legend.position = c(0.05, 0.95), legend.justification = c(0, 1),
  legend.key.height = unit(.75, "line"),
  legend.margin = margin(b=-.25, unit = "line"),
  legend.title = element_text(size=rel(0.9))
)

res <- ts / ar / cfr + plot_annotation(tag_levels = "A") & 
  theme(text = element_text(size = 7), panel.grid.minor = element_blank())

ggsave(tail(.args, 1), res, height = 6, width = 5, units = "in", dpi = 300)
