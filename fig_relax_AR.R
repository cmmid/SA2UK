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

AR.dt <- proj.dt[,
  .(count = sum(value)),
  keyby=.(sample, sentinel = between(as.numeric(group), 5, 10), date)
][capita, on=.(sentinel), rate:=count/pop ][,{
  qs <- quantile(rate, probs = c(0.025, 0.25,0.5,0.75,0.975))
  names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
  as.list(qs)
}, by=.(date, sentinel)
]

serosurvey <- data.table(
  sentinel = c(TRUE, FALSE),
  start = as.Date("2020-07-15"),
  mid = mean(as.Date(c("2020-07-15", "2020-08-07"))),
  end = as.Date("2020-08-07"),
  hi = c(.431, 0.402),
  md = c(.411, .355),
  lo = c(.391, 0.311)
)

AR.p <- force(ggplot(AR.dt) + aes(date, fill = sentinel) +
  geom_ribbon(aes(ymin=lo.lo, ymax=hi.hi), alpha = 0.1) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha = 0.25) +
  geom_line(aes(y=med, color = sentinel), alpha = 0.4) +
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
    name = "Serosurvey",
    breaks = c(TRUE, FALSE),
    labels = c("20-49", "<20 & 50+"),
    values = c(`TRUE`="goldenrod", `FALSE`="grey"),
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
    scale_fill_manual(NULL,
      breaks = c(TRUE, FALSE),
      labels = c("20-49", "<20 & 50+"),
      values = c(`TRUE`="goldenrod", `FALSE`="grey"), guide = "none"
    ) +
  scale_y_continuous("Cumulative Attack Fraction", labels = function(b) sprintf("%2g%%",b*100)) +
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
