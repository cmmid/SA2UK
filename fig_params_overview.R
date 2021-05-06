suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper","PAK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/projections/%s.rds",
  "%s/outputs/adj_data.rds",
  .debug[2], # PAK
  "%s/outputs/params/%s.png"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

est <- readRDS(.args[1])[
  compartment == "cases",
  .(rv = sum(rv), value = sum(value), asc = unique(asc)),
  by=.(sample, date)
]

tariso <- tail(.args, 2)[1]

case.dt <- readRDS(.args[2])[
  (iso3 == tariso),
  .(date, croll = frollmean(cases, 7, align = "center"))
]

parplot <- ggplot(est2) +
 aes(date, rv, group = sample) +
 geom_line(aes(color="sim. total"), alpha = 0.2) +
 geom_line(aes(y=value, color="sim. ascertained"), alpha = 0.2) +
 geom_line(
    aes(date, croll, color = "observed"),
    data = case.dt,
    inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = as.Date(c("2020-03-01",NA))) +
 scale_color_manual(
   name=NULL,
   breaks = c("sim. total", "sim. ascertained", "observed"),
   values = c("firebrick", "goldenrod", "black")
 ) +
 scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
 scale_y_log10(
   "Cases",
   breaks = 10^c(0,3,6),
   minor_breaks = 10^(1:6),
   labels = scales::label_number_si()
  ) +
 theme_minimal()

ggsave(tail(.args, 1), parplot, height = 6.5, width = 15, units = "in", dpi = 900)
