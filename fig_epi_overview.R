suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/Covid_LMIC/All_Africa_paper", "GHA")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/epi_data.rds",
  .debug[2],
  "%s/outputs/figs/epi/%s.png"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

outcomes <- readRDS(.args[1])[iso3 == tariso]

roll.window <- 7

outcomes[order(date), rollcases := frollmean(cases, roll.window), by=iso3 ]
outcomes[order(date), rolldeaths := frollmean(deaths, roll.window), by=iso3 ]
outcomes$cases <- as.numeric(outcomes$cases)
outcomes$deaths <- as.numeric(outcomes$deaths)

mlt <- melt(
  outcomes, id.vars = "date", measure.vars = c("cases", "deaths", "rollcases", "rolldeaths")
)[!is.na(value)]

mlt[, outcome := fifelse(grepl("case", variable), "cases", "deaths")]
mlt[, measure := fifelse(grepl("roll", variable), "rolling", "raw")]
#mlt[, variant := fifelse(grepl("var", variable), "variant", "all")]

p <- ggplot(mlt) +
  aes(date, y= value, alpha = measure) +
  facet_grid(outcome ~ ., scale = "free_y", switch = "y", labeller = labeller(
    outcome = function(o) sprintf("%s (new per day)", o)
  )) +
  geom_line() +
  geom_blank(aes(alpha = NULL), data = data.table(date=as.Date("2020-01-22"), outcome = c("cases","deaths"), value = c(1600, 30))) +
  theme_minimal(base_size = 18) +
  theme(
    strip.placement = "outside", panel.spacing.y = unit(1, "line"),
    legend.position = c(0.05,.95), legend.justification = c(0,1)
  ) +
  coord_cartesian(
    xlim = as.Date(c("2020-03-01", "2021-04-01")),
    expand = F
  ) +
  scale_x_date(NULL,
    date_breaks = "3 months", date_minor_breaks = "month", date_labels = "%b %Y"
  ) +
  scale_y_continuous(NULL, minor_breaks = NULL) +
  scale_alpha_discrete(NULL)

ggsave(tail(.args, 1), p, width = 16, height = 8, dpi = 900)
ggsave(gsub("\\.(\\w+)$","_log.\\1", tail(.args,1)),
  p + lims(y=c(0.1, NA)) + scale_y_log10(NULL, labels = scales::label_number_si(accuracy = 0.1)),
  width = 16, height = 8, dpi = 900
)  
