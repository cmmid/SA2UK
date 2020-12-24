suppressPackageStartupMessages({
  require(data.table)
  require(lubridate)
  require(ggplot2)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/epi_data.rds",
  "%s/outputs/intervention_timing/%s.rds",
  .debug[2],
  "%s/outputs/intervention_timing/%s.png"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

outcomes <- readRDS(.args[1])[iso3 == tariso]

eras <- readRDS(.args[2])

p <- ggplot(outcomes) +
  aes(date) +
  geom_line(aes(y=cases, linetype="cases")) +
  geom_line(aes(y=deaths, linetype="deaths")) +
  geom_rect(
    aes(
      ymin = 0.1, ymax = Inf,
      xmin=start-0.5, xmax=end+0.5,
      fill = era
    ), data = eras, inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_vline(
    xintercept = eras[era == "pre", end + 3],
    linetype = "dashed", color = "dodgerblue") +
  theme_minimal() +
  scale_y_log10("incidence", breaks = 10^(0:5), labels = scales::label_number_si()) +
  scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
  scale_fill_manual(
    name = "fitting era",
    breaks=c("censor","pre","transition","post","modification","variant"),
    values = c(censor="grey", pre="firebrick", transition="grey",post="dodgerblue",modification="goldenrod",variant="red")
  ) +
  scale_linetype_discrete(name="outcome") +
  coord_cartesian(ylim = c(1, NA), xlim = as.Date(c("2020-01-01","2021-01-01")))
  
ggsave(tail(.args, 1), p, height = 3, width = 6, units = "in", dpi = 300)