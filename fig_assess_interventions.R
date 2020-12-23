suppressPackageStartupMessages({
  require(data.table)
  require(lubridate)
  require(ggplot2)
})

.debug <- c("~/Dropbox/covidLMIC", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/outputs/interventions.rds",
  "%s/inputs/ecdc_data.rds",
#  "%s/outputs/introductions/%s.rds",
  .debug[2],
  "%s/outputs/intervention_timing/%s.png"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

outcomes <- readRDS(.args[2])[iso3 == tariso]

min.date <- as.Date("2020-02-01")

ref.time <- readRDS(.args[1])[iso3 == tariso]

# intros <- readRDS(.args[3])

if (dim(ref.time)[1]) {
  transition <- ref.time[, .(iso3, start = as.Date(date + intervention), end = as.Date(date + 6 + intervention), era = "transition")]
  post <- copy(transition)[, .(iso3, start = as.Date(end+1), end = as.Date(end + 31), era = "post")]
  censor <- outcomes[which.max(cumsum(cases) > 10)][, .(iso3, start = min.date, end = date + ref.time$censor, era = "censor")]
  pre <- censor[transition[, .(iso3, end = start)], on=.(iso3)][, .(iso3, start = as.Date(end+1), end = as.Date(i.end-1), era = "pre")]
  
#  deaths <- intros[, .(iso3, start = date[1], end = date[.N], era = "intros")]
  
  eras <- rbind(
    censor, pre, transition, post#, deaths
  )
  
} else {
  eras <- data.table(iso3=character(), start = Date(), end = Date(), era = character())
}

p <- ggplot(outcomes) +
  aes(date) +
  geom_line(aes(y=cases)) +
  geom_line(aes(y=deaths), linetype = "dotted") +
  { if (dim(eras)[1]) list(geom_rect(
    aes(
      ymin = 0.1, ymax = Inf,
      xmin=start-0.5, xmax=end+0.5,
      fill = era
    ), data = eras, inherit.aes = FALSE,
    alpha = 0.2
  ), geom_vline(xintercept = eras[era == "pre", end + 3], linetype = "dashed", color = "dodgerblue")) else geom_blank() } +
  theme_minimal() +
  theme(
    legend.position = c(0, 1), legend.justification = c(0,1),
    legend.background = element_rect(fill = "white")
  ) +
  coord_cartesian(ylim = c(1, NA), xlim = as.Date(c("2020-02-01","2020-10-01")), expand = FALSE) +
  scale_y_log10(breaks = 10^(0:5), labels = scales::label_number_si()) +
  scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b")

ggsave(tail(.args, 1), p, height = 3, width = 6, units = "in", dpi = 300)