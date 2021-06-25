suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("analysis","PAK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/history/%s.rds",
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

#' @examples 
#' TODO: from debugging
#' ggplot(est[compartment == "cases"][, .(asc.value = sum(value), rv = sum(rv)), by=.(sample, date)]) +
#'   aes(date, asc.value, group = sample) +
#'   geom_line(alpha = 0.1) +
#'   geom_line(aes(y=rv), alpha = 0.1, color = "red") +
#'   geom_line(
#'     aes(date, adj),
#'     data = readRDS(file.path(.debug[1], "ins", "adj_data.rds"))[iso3 == .debug[2]],
#'     color = "black", inherit.aes = FALSE
#'   ) +
#'   geom_line(
#'     aes(date, value),
#'     data = readRDS(file.path(.debug[1], "est", "r0", paste0(.debug[2],".rds")))[variable == "infections" & sample < 100],
#'     color = "green", alpha = 0.1
#'   ) +
#'   annotate("rect", xmin=as.Date("2020-09-01"), xmax =as.Date("2020-10-01"), ymin = 0.01, ymax = Inf, alpha = 0.2, fill = "dodgerblue") +
#'   scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
#'   scale_y_log10(labels = scales::label_number_si()) +
#'   coord_cartesian(ylim = c(1, NA)) +
#'   theme_minimal()

parplot <- ggplot(est) +
 aes(date, rv, group = sample) +
  geom_line(aes(color="sim. total"), alpha = 0.05) +
  geom_line(
    aes(color="sim. total", group = NULL),
    data = function(dt) dt[, .(rv=median(rv)), by=.(date)]
  ) +
  geom_line(aes(y=value, color="sim. ascertained"), alpha = 0.05) +
  geom_line(
    aes(y=value, color="sim. ascertained", group = NULL),
    data = function(dt) dt[, .(value=median(value)), by=.(date)]
  ) +
 geom_line(
    aes(date, croll, color = "observed"),
    data = case.dt,
    inherit.aes = FALSE
  ) +
 coord_cartesian(xlim = as.Date(c("2020-03-01",NA)), ylim = c(1, 1e6), expand = F) +
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
