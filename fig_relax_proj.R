suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

.debug <- c("~/Dropbox/SA2UK", "ZAF")
.args <- if (interactive()) sprintf(c(
  "%s/inputs/epi_data.rds",
  "%s/outputs/intervention_timing/%s.rds",
  "%s/outputs/phylo.rds",
  "%s/outputs/projections/%s.rds",
  "%s/outputs/params/%s.rds",
  .debug[2],
  "%s/outputs/figs/timeseries.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = TRUE)

tariso <- tail(.args, 2)[1]

outcomes <- readRDS(.args[1])[iso3 == tariso]

eras <- readRDS(.args[2])

phylo.frac <- readRDS(.args[3])

outcomes[, var.frac := 0 ]
outcomes[phylo.frac[!is.na(binop)], on=.(date), var.frac := binop ]
outcomes[!between(date, phylo.frac[, min(date)], phylo.frac[, max(date)]), var.frac := NA ]

death.delay <- 21

outcomes[, del.var.frac := 0 ]
outcomes[phylo.frac[!is.na(binop), .(date = date + 21, binop)], on=.(date), del.var.frac := binop ]
outcomes[!between(date, phylo.frac[, min(date)+21], phylo.frac[, max(date)+21]), del.var.frac := NA ]

roll.window <- 7

outcomes$cases <- as.numeric(outcomes$cases)
outcomes$deaths <- as.numeric(outcomes$deaths)
outcomes[order(date), rollcases := frollmean(cases, roll.window), by=iso3 ]
outcomes[order(date), rolldeaths := frollmean(deaths, roll.window), by=iso3 ]
outcomes[, rollvarcases := rollcases * var.frac ]
outcomes[, rollvardeaths := rolldeaths * del.var.frac ]

mlt <- melt(
  outcomes, id.vars = "date", measure.vars = c("cases", "deaths", "rollcases", "rolldeaths", "rollvarcases", "rollvardeaths")
)[!is.na(value)]

mlt[, outcome := fifelse(grepl("case", variable), "cases", "deaths")]
mlt[, measure := fifelse(grepl("roll", variable), "rolling", "raw")]
mlt[, variant := fifelse(grepl("var", variable), "variant", "all")]

proj.dt <- readRDS(.args[4])[compartment %in% c("cases","death_o")][, .(value = sum(value)), by=.(date, compartment, sample)]
proj.dt[,
  outcome := fifelse(compartment == "cases", "cases", "deaths") 
][, measure := "rolling" ][, variant := "projection" ]

asc.dt <- readRDS(.args[5])[, .(sample, detectp = asc)]
proj.dt[asc.dt, on=.(sample), asc := fifelse(outcome == "cases", value*detectp, value*.5)]

plot.proj <- proj.dt[, {
  qs <- quantile(asc, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
  as.list(qs)
}, by=.(date, outcome, variant, measure)]

p <- force(ggplot() +
  aes(date, y = value, color = variant, alpha = measure, linetype = outcome) +
  geom_ribbon(aes(ymin=lo.lo, ymax=hi.hi,y=NULL,color=NULL,fill=variant), data = plot.proj, alpha = 0.1, show.legend = FALSE) +
  geom_ribbon(aes(ymin=lo, ymax=hi,y=NULL,color=NULL,fill=variant), data = plot.proj, alpha = 0.25, show.legend = FALSE) +
  geom_line(aes(y=med), data = plot.proj) +
  geom_line(data = mlt[measure != "raw" & (value > 0.1)]) +
  scale_y_log10(
    sprintf("Incidence", roll.window),
    breaks = 10^(0:5), labels = scales::label_number_si(),
    minor_breaks = NULL
  ) +
  scale_x_date(NULL, date_breaks = "month", date_minor_breaks = "week", date_labels = "%b") +
  scale_color_manual(
    name = NULL,
    labels = c(
      all = "reported (7 day mean)", variant = "estimated 501Y.V2", projection = "model (no 501Y.V2)"
    ),
    values = c(all="black", variant = "red", projection = "dodgerblue"),
    aesthetics = c("color", "fill")
  ) +
  scale_linetype_manual(name = NULL, values = c(cases="solid", deaths = "longdash")) +
  scale_alpha_manual(name = NULL, values = c(raw=0.4, rolling = 1), guide = "none") +
  coord_cartesian(ylim = c(1, NA), xlim = as.Date(c("2020-04-01", NA)), expand = FALSE) +
  theme_minimal())

saveRDS(p, tail(.args, 1))

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
