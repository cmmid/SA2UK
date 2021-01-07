
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(patchwork)
})

.args <- if (interactive()) sprintf(c(
  "%s/phylo.rds", "%s/cfr.rds", "%s/timeseries.rds", "%s/AR.rds"
), "~/Dropbox/SA2UK/outputs/figs") else commandArgs(trailingOnly = TRUE)

zlow <- as.Date("2020-10-15")
zend <- as.Date("2021-01-01")

zoom <- annotate(
  "rect",
  xmin = zlow, xmax = zend,
  ymin = 0, ymax = Inf,
  fill = "grey", alpha = 0.1
)

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
cfr <- readRDS(.args[2]) + theme(
  legend.position = c(0.4, 1), legend.justification = c(0.4, 1),
  legend.key.height = unit(.75, "line")
) + altcfrscale

ts <- readRDS(.args[3]) + theme(
  legend.direction = "horizontal",
  legend.position = c(0.5, 0), legend.justification = c(0.5, 0),
  legend.margin = margin(t=-1, b=-.25, unit = "line")
)
ar <- readRDS(.args[4]) +
coord_cartesian(
  xlim = as.Date(c("2020-04-01","2021-01-01")),
  ylim = c(0, 1),
  expand = FALSE
) + theme(
  legend.position = c(0, 1), legend.justification = c(0, 1),
  legend.key.height = unit(.75, "line")
)

playout <- "
AAAA
BBBB
CCCC
DDEE
"

res <- ts / ar / cfr + plot_annotation(tag_levels = "A") & 
  theme(text = element_text(size = 7), panel.grid.minor = element_blank())

ggsave("combined.png", res, height = 6, width = 5, units = "in", dpi = 300)
