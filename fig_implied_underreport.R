suppressPackageStartupMessages({
  require(data.table)
  require(qs)
  require(ggplot2)
})

.debug <- c("~/Dropbox/covidLMIC", "PAK")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/projections/%s.qs",
  "%s/outputs/introductions/%s.rds",
  "%s/inputs/ecdc_data.rds",
  "%s/outputs/intervention_timing/%s.rds",
  .debug[2],
  sprintf("%s.rds", .debug[2])
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <-tail(.args, 2)[1]

day0 <- readRDS(.args[2])[iso3 == tariso, date[1]]

cases <- melt(
  readRDS(.args[3])[iso3 == tariso],
  id.vars = c("date"),
  measure.vars = c("cases", "deaths"),
  variable.name = "compartment"
)[, c("scen_id", "run") := .(0, 3) ]

lims.dt <- readRDS(.args[4])

dt <- qread(.args[1])[compartment %in% c("cases", "death_o")]
dt[compartment == "death_o", compartment := "deaths" ]
dt[, compartment := factor(compartment, levels = c("cases","deaths"))]
dt[, date := day0 + t ]

allage.dt <- rbind(dt[,
  .(value = sum(value)),
  keyby=.(scen_id, run, compartment, date)
], cases)

reporting_ref <- allage.dt[
  scen_id == 0
][
  order(date), .(date, value=frollmean(value, 7)), by=.(compartment)
][
  !is.na(value)
][date >= dt[, min(date)]]

implied_under_reporting <- allage.dt[
  scen_id > 1,
  .(value = median(value)),
  keyby=.(run, compartment, date)
][
  reporting_ref,
  on = .(compartment, date)
]

rep.dt <- implied_under_reporting[,
  reported_percent := i.value/value
][date >= lims.dt[era == "post", end-7]]

rep.p <- ggplot(rep.dt) + aes(
  date, reported_percent*100, group = compartment
) +
  geom_line(data = function(dt) dt[run == 3]) +
  # TODO: have to do the more complicated time series thing
  # geom_ribbon(
  #   aes(ymin=`4`*100, ymax=`2`*100, y=NULL),
  #   data = function(dt) dcast(dt, date + compartment ~ run, value.var = "reported_percent"),
  #   alpha = 0.1
  # ) +
  theme_minimal() +
  scale_x_date(name = NULL, date_breaks = "months", date_labels = "%b", date_minor_breaks = "weeks") +
  scale_y_continuous(expression("implied reported %"))

saveRDS(rep.p, tail(.args, 1))

#' @examples 
#' p <- rep.p + aes(color = compartment) + 
#'   theme(legend.position = c(.9,.95), legend.justification = c(1,1)) +
#'   scale_color_discrete(NULL)
#' ggsave("reporting.png", p, width = 8, height = 5, units = "in")
#' 