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

gm_roll <- function(x, n) prod(x, na.rm = TRUE)^(1/sum(!is.na(x)))

outcomes[order(date), rollcases := frollmean(cases, roll.window), by=iso3 ]
outcomes[order(date), rolldeaths := frollmean(deaths, roll.window), by=iso3 ]
# TODO: think about the appropriate mean here?
# outcomes[order(date), grollcases := { x <- cases; x[x==0] <- NA; frollapply(x, roll.window, gm_roll) }, by=iso3 ]
# outcomes[order(date), grolldeaths := { x <- deaths; x[x==0] <- NA; frollapply(x, roll.window, gm_roll) }, by=iso3 ]

outcomes$cases <- as.numeric(outcomes$cases)
outcomes$deaths <- as.numeric(outcomes$deaths)

mlt <- melt(
  outcomes, id.vars = "date", measure.vars = c("cases", "deaths", "rollcases", "rolldeaths", "grollcases", "grolldeaths")
)[!is.na(value)]

mlt[, outcome := fifelse(grepl("case", variable), "cases", "deaths")]
mlt[, measure := fifelse(grepl("roll", variable), "rolling", "raw")]
# mlt[measure == "rolling", measure := fifelse(grepl("groll", variable), "grolling", "rolling")]
# mlt[, variant := fifelse(grepl("var", variable), "variant", "all")]

mlt[, dow := wday(date) ]

p <- function(dt = mlt, yupper = NA) ggplot(dt) +
  aes(date, y = value, color = outcome) +
  # facet_grid(outcome ~ ., scale = "free_y", switch = "y", labeller = labeller(
  #   outcome = function(o) sprintf("%s (new per day)", o)
  # )) +
  geom_rect(
    aes(ymin = 0, ymax = Inf, xmin = xlo, xmax = xhi),
    data = function(dt) dt[,.(date = unique(date))][, {
      wd <- wday(date)
      i1 <- wd[which.max(wd == 1)]
      xlo <- seq(date[i1], date[.N], by=14)
      xhi <- xlo + 6
      .(xlo = xlo, xhi = xhi)
    }],
    inherit.aes = FALSE,
    fill = "lightgrey", alpha = 0.2
  ) +
  geom_text(
    aes(y = peakv, x = xlo, label = md),
    data = function(dt) dt[,.(date = unique(date), peakv = max(value))][, {
      wd <- wday(date)
      i1 <- wd[which.max(wd == 1)]
      xlo <- seq(date[i1], date[.N], by=14)
      .(xlo = xlo + 3, peakv = 0.3, md = mday(xlo))
    }],
    inherit.aes = FALSE,
    fill = "lightgrey"
  ) +
  geom_line(data = function(dt) dt[measure == "rolling"]) +
#  geom_line(data = function(dt) dt[measure == "grolling"], linetype = "dashed") +
  geom_point(data = function(dt) dt[measure == "raw"], alpha = 0.25) +
#  geom_blank(aes(alpha = NULL), data = data.table(date=as.Date("2020-01-22"), outcome = c("cases","deaths"), value = c(1600, 30))) +
  theme_minimal(base_size = 18) +
  theme(
    strip.placement = "outside", panel.spacing.y = unit(1, "line"),
    legend.position = c(0.05,.95), legend.justification = c(0,1)
  ) +
  coord_cartesian(
    xlim = as.Date(c("2020-03-01", "2021-04-01")),
    ylim = c(NA, yupper),
    expand = F, clip = "off"
  ) +
  scale_x_date(NULL,
    date_breaks = "3 months", date_minor_breaks = "month", date_labels = "%b %Y"
  ) +
  scale_y_log10(NULL, labels = scales::label_number_si(accuracy = 0.1)) +
#  scale_y_continuous(NULL, minor_breaks = NULL) +
  scale_alpha_discrete(NULL) +
  scale_color_manual(name = NULL, values = c(cases="dodgerblue", deaths = "firebrick"))

#ggsave(tail(.args, 1), p(mlt[outcome == "cases" & date >= "2020-03-01"], 1600), width = 16, height = 8, dpi = 300)
ggsave(tail(.args, 1),
  p(mlt[date >= "2020-03-01"]),
  width = 16, height = 8, dpi = 300
)  
