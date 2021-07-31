suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c(".","NGA")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/scenario/%s",
  "%s/outputs/adj_data.rds",
  .debug[2], # PAK
  "%s/outputs/scenarios/%s.png"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

est <- rbindlist(lapply(
  list.files(dirname(.args[1]), sprintf("%s_\\d\\.rds", basename(.args[1])), full.names = TRUE),
  function(fn) readRDS(fn)[compartment == "cases", .(epi_id, sample, date, group, cases = rv)]
))

allages <- est[, .(cases = sum(cases)), by=.(epi_id, sample, date)]
allages[
  date > "2021-09-01", .(
    date, ccases=cumsum(cases)
  ), by=.(epi_id, sample)
]


[date == max(date), as.list(
  quantile(ccases, probs = c(0.025,0.25,0.5,0.75,0.975))
), by=epi_id]
ggplot(allages) + aes(date, cases, group=sample) +
  facet_grid(epi_id ~ .) +
  geom_line(alpha = 0.1) +
  scale_x_date() +
  scale_y_log10(labels = scales::label_number_si()) +
  coord_cartesian(xlim = as.Date(c("2021-09-01","2022-01-01")), ylim = c(1e3, NA)) +
  theme_minimal()

readRDS(.args[1])[
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
