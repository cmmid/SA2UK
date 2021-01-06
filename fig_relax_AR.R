suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/urbanization.rds",
  "%s/inputs/pops/%s.rds",
  "%s/outputs/projections/%s.rds",
  .debug[2],
  "%s/outputs/figs/AR.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

proj.dt <- readRDS(.args[3])[compartment == "R"]

urbfrac <- readRDS(.args[1])[iso3 == tariso, value/100 ]
pop <- round(readRDS(.args[2])$pop[[1]]$size*urbfrac)
capita <- data.table(sentinel = c(TRUE, FALSE), pop = c(sum(pop[5:10]),sum(pop[-c(5:10)])))

AR.dt <- dcast(proj.dt[,
  .(count = sum(value)),
  keyby=.(q, sentinel = between(as.numeric(group), 5, 10), date)
][capita, on=.(sentinel), rate:=count/pop ], sentinel + date ~ q, value.var = "rate")

#' TODO figure out how to force evaluation of this in ggplot
seroposdelay <- 14
serosurvey <- data.table(
  sentinel = TRUE,
  start = as.Date("2020-07-15")-seroposdelay,
  mid = mean(as.Date(c("2020-07-15", "2020-08-07")))-seroposdelay,
  end = as.Date("2020-08-07")-seroposdelay,
  hi = .431,
  md = .411,
  lo = .391
)

AR.p <- force(ggplot(AR.dt) + aes(date, fill = sentinel) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha = 0.4) +
  geom_line(aes(y=md, color = sentinel)) +
  geom_linerange(
    aes(
      y = md, x = NULL,
      xmin = start, xmax = end,
      color = sentinel
    ), data = serosurvey
  ) +
  geom_linerange(
    aes(
      x = mid, y = NULL,
      ymin = lo, ymax = hi,
      color = sentinel
    ), data = serosurvey
  ) +
  scale_color_manual(
    name = NULL,
    breaks = c(TRUE, FALSE),
    labels = c("20-49", "<20 & 50+"),
    values = c(`TRUE`="goldenrod", `FALSE`="grey"),
    aesthetics = c("color", "fill")
  ) +
  scale_y_continuous("Cumulative Attack Fraction") +
  scale_x_date(name = NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
  theme_minimal())


saveRDS(AR.p, tail(.args, 1))

# geom_rect(
#   aes(
#     ymin = 0.1, ymax = Inf,
#     xmin=start-0.5, xmax=end+0.5,
#     fill = era
#   ), data = eras[!(era %in% c("censor","transition"))], inherit.aes = FALSE,
#   alpha = 0.1
# ) +
# scale_fill_manual(
#   name = NULL,
#   breaks=c("pre","post","relaxation","variant"),
#   labels=c(pre="pre-intervention",post="post-intervention",relaxation="relaxation",variant="emergent variant"),
#   values = c(pre="firebrick", post="dodgerblue",relaxation="goldenrod",variant="red")
# ) +
