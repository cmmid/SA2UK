suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(patchwork)
})

.args <- if (interactive()) sprintf(c(
  "%s/phylo.rds", "%s/eras.rds",
  "%s/SI1.png"
), "~/Dropbox/SA2UK/outputs/figs") else commandArgs(trailingOnly = TRUE)

zlow <- as.Date("2020-10-15")
zend <- as.Date("2020-12-01")

phylo <- readRDS(.args[1]) + coord_cartesian(
  xlim = c(zlow, zend),
  ylim = c(0, 1),
  expand = FALSE
)

interventions <- readRDS(.args[2])

playout <- c("ABBBB")

res <- phylo + interventions + plot_layout(design = playout) +
  plot_annotation(tag_levels = "A") & 
  theme(text = element_text(size = 7), panel.grid.minor = element_blank())

ggsave(tail(.args, 1), res, height = 4, width = 7, units = "in", dpi = 300)
