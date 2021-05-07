suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper","PAK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/scenarios/%s.rds",
  "%s/outputs/adj_data.rds",
  .debug[2], # PAK
  "%s/outputs/scenarios/%s.png"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

est <- readRDS(.args[1])[
  compartment == "cases",
  .(rv = sum(rv), value = sum(value), asc = unique(asc)),
  by=.(scenario, sample, date)
]

tariso <- tail(.args, 2)[1]

case.dt <- readRDS(.args[2])[
  (iso3 == tariso),
  .(date, croll = frollapply(adj, 7, function(x) prod(x)^(1/7), align = "center"))
]

parplot <- ggplot(est) +
 facet_grid(scenario ~ ., labeller = labeller(
   scenario = c(
     none = "No Change", noschool = "No School", peakworkred = "Highest Work-Based",
     peakred = "Peak NPIs", vaccination = "Vaccination Only"
   )
 )) +
 aes(date, rv, group = sample) +
 geom_line(aes(color="sim. total"), alpha = 0.05) +
 geom_line(
   aes(color="sim. total", group = NULL),
   data = function(dt) dt[, .(rv=median(rv)), by=.(scenario, date)]
 ) +
 geom_line(aes(y=value, color="sim. ascertained"), alpha = 0.05) +
  geom_line(
  aes(y=value, color="sim. ascertained", group = NULL),
  data = function(dt) dt[, .(value=median(value)), by=.(scenario, date)]
  ) +
 geom_line(
    aes(date, croll, color = "observed"),
    data = case.dt,
    inherit.aes = FALSE
  ) +
 coord_cartesian(xlim = as.Date(c("2020-03-01",NA)), ylim=c(1,1e6), expand = F) +
 scale_color_manual(
   name=NULL,
   breaks = c("sim. total", "sim. ascertained", "observed"),
   values = c("firebrick", "goldenrod", "black")
 ) +
 scale_x_date(
   NULL,
   date_breaks = "month", date_minor_breaks = "week",
   labels = function(d) c("J","F","M","A","M","J","J","A","S","O","N","D")[data.table::month(d)]
 ) +
 scale_y_log10(
   "Cases (per day)",
   breaks = 10^c(0,3,6),
   minor_breaks = 10^(1:6),
   labels = scales::label_number_si()
  ) +
 theme_minimal()

ggsave(tail(.args, 1), parplot, height = 6.5, width = 7.5, units = "in", dpi = 900)
